#!/usr/bin/env Rscript
# ============================================================================
# Script: 4a_make-short-id-gmts.R
# Purpose: Create short-ID GMT files for MAGMA pathway analysis.
#
# Why this is needed:
#   MAGMA can truncate long pathway names in .gsa.out files, e.g.
#     GOBP_POSITIVE_REGULATION_OF_...670
#   Truncated names are hard to match back to the original GMT files.
#   This script creates GMT files where the pathway name is a short stable ID:
#     GO_000001, GO_000002, ..., REACTOME_000001, ...
#   It also writes a manifest mapping each short ID back to the full pathway name.
#
# Run before:
#   4_run-MAGMA-pathways.sh
#
# Then in script 4 use:
#   pathways/short_ids/c5.go.v2026.1.Hs.entrez.shortids.gmt
#   pathways/short_ids/c2.cp.reactome.v2026.1.Hs.entrez.shortids.gmt
#
# For downstream scoring, script 6a should use the corresponding SYMBOL short-ID
# GMT files produced here.
# ============================================================================

module -f unload compilers mpi gcc-libs
module load r/recommended
R --vanilla

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

pathway_dir <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/pathways"
out_dir <- file.path(pathway_dir, "short_ids")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# GO input files.
# If a combined c5.go file exists, use it. Otherwise use BP/CC/MF files.
go_entrez_combined <- file.path(pathway_dir, "c5.go.v2026.1.Hs.entrez.gmt")
go_entrez_split <- c(
  file.path(pathway_dir, "c5.go.bp.v2026.1.Hs.entrez.gmt"),
  file.path(pathway_dir, "c5.go.cc.v2026.1.Hs.entrez.gmt"),
  file.path(pathway_dir, "c5.go.mf.v2026.1.Hs.entrez.gmt")
)

go_symbol_combined <- file.path(pathway_dir, "c5.go.v2026.1.Hs.symbols.gmt")
go_symbol_split <- c(
  file.path(pathway_dir, "c5.go.bp.v2026.1.Hs.symbols.gmt"),
  file.path(pathway_dir, "c5.go.cc.v2026.1.Hs.symbols.gmt"),
  file.path(pathway_dir, "c5.go.mf.v2026.1.Hs.symbols.gmt")
)

# Reactome input files.
reactome_entrez <- file.path(pathway_dir, "c2.cp.reactome.v2026.1.Hs.entrez.gmt")
reactome_symbol <- file.path(pathway_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt")

# Output files.
go_entrez_out <- file.path(out_dir, "c5.go.v2026.1.Hs.entrez.shortids.gmt")
go_symbol_out <- file.path(out_dir, "c5.go.v2026.1.Hs.symbols.shortids.gmt")
reactome_entrez_out <- file.path(out_dir, "c2.cp.reactome.v2026.1.Hs.entrez.shortids.gmt")
reactome_symbol_out <- file.path(out_dir, "c2.cp.reactome.v2026.1.Hs.symbols.shortids.gmt")
manifest_out <- file.path(out_dir, "pathway_short_id_manifest.tsv")
missing_symbol_out <- file.path(out_dir, "pathway_short_id_missing_symbol_terms.tsv")

# ------------------------------------------------------------
# 2. Basic functions
# ------------------------------------------------------------

choose_gmts <- function(combined, split) {
  if (file.exists(combined)) return(combined)
  split[file.exists(split)]
}

read_gmt <- function(paths) {
  rows <- list()
  k <- 1

  for (path in paths) {
    if (!file.exists(path)) next

    lines <- readLines(path, warn = FALSE)

    for (i in seq_along(lines)) {
      parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
      if (length(parts) < 3) next

      genes <- unique(parts[-c(1, 2)])
      genes <- genes[nzchar(genes)]

      rows[[k]] <- data.frame(
        full_name = parts[1],
        description = parts[2],
        genes = paste(genes, collapse = "\t"),
        n_genes = length(genes),
        source_gmt = basename(path),
        source_line = i,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      full_name = character(), description = character(), genes = character(),
      n_genes = integer(), source_gmt = character(), source_line = integer(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$full_name), ]
  row.names(out) <- NULL
  out
}

write_gmt <- function(tab, path, gene_column = "genes") {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  if (nrow(tab) == 0) return(invisible(NULL))

  for (i in seq_len(nrow(tab))) {
    writeLines(
      paste(tab$short_id[i], tab$short_description[i], tab[[gene_column]][i], sep = "\t"),
      con
    )
  }
}

make_short_ids <- function(collection, id_prefix, entrez_paths, symbol_paths,
                           entrez_out, symbol_out) {

  cat("\n============================================================\n")
  cat("Processing collection:", collection, "\n")
  cat("Entrez GMT files:\n")
  cat(paste("  ", entrez_paths, collapse = "\n"), "\n", sep = "")
  cat("Symbol GMT files:\n")
  cat(paste("  ", symbol_paths, collapse = "\n"), "\n", sep = "")

  entrez <- read_gmt(entrez_paths)
  symbols <- read_gmt(symbol_paths)

  if (nrow(entrez) == 0) {
    cat("WARNING: no Entrez pathways read for ", collection, ".\n", sep = "")
    return(list(manifest = data.frame(), missing = data.frame()))
  }

  if (nrow(symbols) == 0) {
    cat("WARNING: no symbol pathways read for ", collection, ".\n", sep = "")
  }

  entrez$short_id <- sprintf("%s_%06d", id_prefix, seq_len(nrow(entrez)))

  m <- match(entrez$full_name, symbols$full_name)

  short_description <- paste0(
    "full_name=", entrez$full_name,
    ";collection=", collection,
    ";original_desc=", entrez$description,
    ";source_entrez_gmt=", entrez$source_gmt
  )

  entrez_short <- data.frame(
    short_id = entrez$short_id,
    short_description = short_description,
    genes = entrez$genes,
    stringsAsFactors = FALSE
  )

  matched <- !is.na(m)

  symbol_short <- data.frame(
    short_id = entrez$short_id[matched],
    short_description = paste0(
      "full_name=", entrez$full_name[matched],
      ";collection=", collection,
      ";original_desc=", symbols$description[m[matched]],
      ";source_symbol_gmt=", symbols$source_gmt[m[matched]]
    ),
    genes = symbols$genes[m[matched]],
    stringsAsFactors = FALSE
  )

  write_gmt(entrez_short, entrez_out)
  write_gmt(symbol_short, symbol_out)

  manifest <- data.frame(
    collection = collection,
    short_id = entrez$short_id,
    full_name = entrez$full_name,
    entrez_description = entrez$description,
    entrez_source_gmt = entrez$source_gmt,
    entrez_source_line = entrez$source_line,
    n_entrez_genes = entrez$n_genes,
    has_symbol_match = matched,
    symbol_source_gmt = ifelse(matched, symbols$source_gmt[m], NA_character_),
    symbol_source_line = ifelse(matched, symbols$source_line[m], NA_integer_),
    n_symbol_genes = ifelse(matched, symbols$n_genes[m], NA_integer_),
    stringsAsFactors = FALSE
  )

  missing <- manifest[!manifest$has_symbol_match, c("collection", "short_id", "full_name", "entrez_source_gmt")]

  cat("Pathways in Entrez GMT:", nrow(entrez), "\n")
  cat("Pathways with symbol match:", sum(matched), "\n")
  cat("Pathways missing from symbol GMT:", sum(!matched), "\n")
  cat("Wrote Entrez short-ID GMT:", entrez_out, "\n")
  cat("Wrote symbol short-ID GMT:", symbol_out, "\n")

  list(manifest = manifest, missing = missing)
}

# ------------------------------------------------------------
# 3. Run
# ------------------------------------------------------------

go_entrez_paths <- choose_gmts(go_entrez_combined, go_entrez_split)
go_symbol_paths <- choose_gmts(go_symbol_combined, go_symbol_split)

reactome_entrez_paths <- reactome_entrez[file.exists(reactome_entrez)]
reactome_symbol_paths <- reactome_symbol[file.exists(reactome_symbol)]

all_manifest <- list()
all_missing <- list()

res_go <- make_short_ids(
  collection = "GO",
  id_prefix = "GO",
  entrez_paths = go_entrez_paths,
  symbol_paths = go_symbol_paths,
  entrez_out = go_entrez_out,
  symbol_out = go_symbol_out
)
all_manifest[[length(all_manifest) + 1]] <- res_go$manifest
all_missing[[length(all_missing) + 1]] <- res_go$missing

res_reactome <- make_short_ids(
  collection = "Reactome",
  id_prefix = "REACTOME",
  entrez_paths = reactome_entrez_paths,
  symbol_paths = reactome_symbol_paths,
  entrez_out = reactome_entrez_out,
  symbol_out = reactome_symbol_out
)
all_manifest[[length(all_manifest) + 1]] <- res_reactome$manifest
all_missing[[length(all_missing) + 1]] <- res_reactome$missing

manifest <- do.call(rbind, all_manifest)
missing <- do.call(rbind, all_missing)

write.table(manifest, manifest_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(missing, missing_symbol_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n============================================================\n")
cat("Done.\n")
cat("Manifest written:", manifest_out, "\n")
cat("Missing symbol report written:", missing_symbol_out, "\n")
cat("\nUse these files for MAGMA pathway analysis:\n")
cat("  GO:       ", go_entrez_out, "\n", sep = "")
cat("  Reactome: ", reactome_entrez_out, "\n", sep = "")
cat("\nUse these files for PRSice/LDpred2 selected pathway scoring:\n")
cat("  GO:       ", go_symbol_out, "\n", sep = "")
cat("  Reactome: ", reactome_symbol_out, "\n", sep = "")
cat("\nAfter creating these files, rerun script 4 using the short-ID Entrez GMTs.\n")
cat("Then rerun scripts 5 and 6a. Selected pathway IDs should match exactly.\n")
