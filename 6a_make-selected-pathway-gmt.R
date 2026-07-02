#!/usr/bin/env Rscript
# ============================================================================
# Script: 6a_make-selected-pathway-gmt.R
# Purpose: Build a GMT containing only the pathways selected from MAGMA results.
#
# Run after:
#   5_extract-sig-pathways.R
#
# Important:
#   Preferred workflow: run 4a_make-short-id-gmts.R before MAGMA pathway
#   analysis, then run MAGMA using the short-ID Entrez GMTs. This avoids
#   MAGMA truncating long pathway names in .gsa.out files.
#
#   This script then matches selected short IDs such as GO_000001 or
#   REACTOME_000001 back to short-ID SYMBOL GMTs for PRSice/LDpred2 scoring.
#   It also still works with old full-name results where names were not
#   truncated, but that is less robust.
# ============================================================================


module -f unload compilers mpi gcc-libs
module load r/recommended
R --vanilla

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

results_dir <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/results"
pathway_dir <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/pathways"
out_dir     <- file.path(pathway_dir, "selected")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Choose one trait.
#trait_prefix <- "HT_EUR_2022"
trait_prefix <- "BMI_EUR_2018"

# Choose selected table from script 5.
# Common values: "FDR05", "Bonferroni05", "nominal_P05", "top50"
selection <- "Bonferroni05"

run_go <- TRUE
run_reactome <- TRUE

# If the selected table is empty, set fallback_top_n > 0 to use the top N
# pathways from the corresponding top50 table. Use this only as exploratory.
fallback_top_n <- 0

# Symbol GMTs. Short-ID files are tried first.
# These are created by 4a_make-short-id-gmts.R and are preferred.
go_symbol_gmts <- c(
  file.path(pathway_dir, "short_ids/c5.go.v2026.1.Hs.symbols.shortids.gmt"),
  file.path(pathway_dir, "c5.go.v2026.1.Hs.symbols.gmt"),
  file.path(pathway_dir, "c5.go.bp.v2026.1.Hs.symbols.gmt"),
  file.path(pathway_dir, "c5.go.cc.v2026.1.Hs.symbols.gmt"),
  file.path(pathway_dir, "c5.go.mf.v2026.1.Hs.symbols.gmt")
)
go_symbol_gmts <- go_symbol_gmts[file.exists(go_symbol_gmts)]

reactome_symbol_gmts <- c(
  file.path(pathway_dir, "short_ids/c2.cp.reactome.v2026.1.Hs.symbols.shortids.gmt"),
  file.path(pathway_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt")
)
reactome_symbol_gmts <- reactome_symbol_gmts[file.exists(reactome_symbol_gmts)]

short_id_manifest_file <- file.path(pathway_dir, "short_ids/pathway_short_id_manifest.tsv")

# ------------------------------------------------------------
# 2. Small utility code
# ------------------------------------------------------------

empty_gmt_table <- function() {
  data.frame(
    name = character(),
    description = character(),
    genes = character(),
    source_gmt = character(),
    collection = character(),
    stringsAsFactors = FALSE
  )
}

read_table_auto <- function(path) {
  if (!file.exists(path)) return(NULL)

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

first_existing <- function(paths) {
  x <- paths[file.exists(paths)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

read_short_id_manifest <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(short_id = character(), full_name = character(), collection = character(), stringsAsFactors = FALSE))
  }

  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("short_id", "full_name", "collection")
  if (!all(needed %in% names(x))) {
    return(data.frame(short_id = character(), full_name = character(), collection = character(), stringsAsFactors = FALSE))
  }

  x[, needed, drop = FALSE]
}

read_gmt <- function(paths) {
  all_rows <- list()
  k <- 1

  for (path in paths) {
    if (!file.exists(path)) next
    lines <- readLines(path, warn = FALSE)

    for (i in seq_along(lines)) {
      parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
      if (length(parts) < 3) next

      all_rows[[k]] <- data.frame(
        name = parts[1],
        description = parts[2],
        genes = paste(unique(parts[-c(1, 2)]), collapse = "\t"),
        source_gmt = basename(path),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  if (length(all_rows) == 0) {
    return(data.frame(
      name = character(), description = character(), genes = character(),
      source_gmt = character(), stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, all_rows)
  out <- out[!duplicated(out$name), ]
  row.names(out) <- NULL
  out
}

write_gmt <- function(x, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  if (nrow(x) == 0) return(invisible(NULL))

  for (i in seq_len(nrow(x))) {
    writeLines(paste(x$name[i], x$description[i], x$genes[i], sep = "\t"), con)
  }
}

get_selected_names <- function(collection) {
  possible <- c(
    file.path(results_dir, paste0(trait_prefix, "_", collection, "_", selection, ".csv")),
    file.path(results_dir, paste0(trait_prefix, "_", collection, "_", selection, ".tsv"))
  )

  selected_file <- first_existing(possible)

  if (is.na(selected_file)) {
    cat("WARNING: no selected table found for ", collection, ":\n", sep = "")
    cat(paste(possible, collapse = "\n"), "\n")
    return(character())
  }

  tab <- read_table_auto(selected_file)

  if (is.null(tab) || nrow(tab) == 0) {
    cat("WARNING: selected table is empty for ", collection, ": ", selected_file, "\n", sep = "")

    if (fallback_top_n > 0) {
      top_possible <- c(
        file.path(results_dir, paste0(trait_prefix, "_", collection, "_top50.csv")),
        file.path(results_dir, paste0(trait_prefix, "_", collection, "_top50.tsv"))
      )
      top_file <- first_existing(top_possible)
      if (!is.na(top_file)) {
        top_tab <- read_table_auto(top_file)
        if (!is.null(top_tab) && "VARIABLE" %in% names(top_tab)) {
          cat("Using top ", fallback_top_n, " pathways from: ", top_file, "\n", sep = "")
          return(unique(as.character(head(top_tab$VARIABLE, fallback_top_n))))
        }
      }
    }

    return(character())
  }

  if (!"VARIABLE" %in% names(tab)) {
    cat("WARNING: VARIABLE column not found in ", selected_file, "\n", sep = "")
    return(character())
  }

  unique(as.character(tab$VARIABLE))
}

# ------------------------------------------------------------
# 3. Extract selected GO and Reactome terms from symbol GMTs
# ------------------------------------------------------------

cat("Trait:", trait_prefix, "\n")
cat("Selection:", selection, "\n")
cat("Output folder:", out_dir, "\n")
cat("GO symbol GMTs used:\n")
cat(paste("  ", go_symbol_gmts, collapse = "\n"), "\n", sep = "")
cat("Reactome symbol GMTs used:\n")
cat(paste("  ", reactome_symbol_gmts, collapse = "\n"), "\n\n", sep = "")

short_id_manifest <- read_short_id_manifest(short_id_manifest_file)

missing <- data.frame(collection = character(), term = character(), stringsAsFactors = FALSE)
combined <- empty_gmt_table()

n_requested_total <- 0
n_written_total <- 0

if (run_go) {
  go_names <- get_selected_names("GO")
  go_gmt <- read_gmt(go_symbol_gmts)

  if (length(go_symbol_gmts) == 0) {
    cat("WARNING: no GO symbol GMT files found. Downstream scoring needs a symbol GMT.\n")
  }

  selected_go <- go_gmt[go_gmt$name %in% go_names, , drop = FALSE]
  missing_go <- setdiff(go_names, selected_go$name)

  if (nrow(selected_go) > 0) {
    selected_go$collection <- "GO"
    combined <- rbind(combined, selected_go)
  }

  if (length(missing_go) > 0) {
    missing <- rbind(missing, data.frame(collection = "GO", term = missing_go, stringsAsFactors = FALSE))
  }

  go_out <- file.path(out_dir, paste0(trait_prefix, "_selected_", selection, "_GO.symbols.gmt"))
  write_gmt(selected_go, go_out)

  cat("GO selected terms requested:", length(go_names), "\n")
  cat("GO selected terms written:", nrow(selected_go), "\n")
  cat("GO missing terms:", length(missing_go), "\n")
  cat("GO GMT written:", go_out, "\n\n")

  n_requested_total <- n_requested_total + length(go_names)
  n_written_total <- n_written_total + nrow(selected_go)
}

if (run_reactome) {
  reactome_names <- get_selected_names("Reactome")

  if (length(reactome_symbol_gmts) == 0) {
    cat("WARNING: no Reactome symbol GMT files found. Downstream scoring needs a symbol GMT.\n")
    reactome_gmt <- data.frame(
      name = character(), description = character(), genes = character(),
      source_gmt = character(), stringsAsFactors = FALSE
    )
  } else {
    reactome_gmt <- read_gmt(reactome_symbol_gmts)
  }

  selected_reactome <- reactome_gmt[reactome_gmt$name %in% reactome_names, , drop = FALSE]
  missing_reactome <- setdiff(reactome_names, selected_reactome$name)

  if (nrow(selected_reactome) > 0) {
    selected_reactome$collection <- "Reactome"
    combined <- rbind(combined, selected_reactome)
  }

  if (length(missing_reactome) > 0) {
    missing <- rbind(missing, data.frame(collection = "Reactome", term = missing_reactome, stringsAsFactors = FALSE))
  }

  reactome_out <- file.path(out_dir, paste0(trait_prefix, "_selected_", selection, "_Reactome.symbols.gmt"))
  write_gmt(selected_reactome, reactome_out)

  cat("Reactome selected terms requested:", length(reactome_names), "\n")
  cat("Reactome selected terms written:", nrow(selected_reactome), "\n")
  cat("Reactome missing terms:", length(missing_reactome), "\n")
  cat("Reactome GMT written:", reactome_out, "\n\n")

  n_requested_total <- n_requested_total + length(reactome_names)
  n_written_total <- n_written_total + nrow(selected_reactome)
}

# ------------------------------------------------------------
# 4. Write combined selected panel and provenance files
# ------------------------------------------------------------

if (nrow(combined) > 0) {
  combined <- combined[!duplicated(combined$name), , drop = FALSE]
  row.names(combined) <- NULL
}

# If matching is incomplete, use PARTIAL in the file name to avoid accidental use
# as a complete selected pathway panel.
partial <- n_requested_total != n_written_total
suffix <- if (partial) "PARTIAL" else ""
if (nzchar(suffix)) suffix <- paste0("_", suffix)

combined_out <- file.path(out_dir, paste0(trait_prefix, "_selected_", selection, "_GO_Reactome", suffix, ".symbols.gmt"))
manifest_out <- file.path(out_dir, paste0(trait_prefix, "_selected_", selection, "_manifest.tsv"))
missing_out <- file.path(out_dir, paste0(trait_prefix, "_selected_", selection, "_missing_terms.tsv"))

write_gmt(combined, combined_out)

if (nrow(combined) > 0) {
  n_genes <- vapply(strsplit(combined$genes, "\t", fixed = TRUE), length, integer(1))
  full_name <- combined$name
  if (nrow(short_id_manifest) > 0) {
    idx <- match(combined$name, short_id_manifest$short_id)
    full_name[!is.na(idx)] <- short_id_manifest$full_name[idx[!is.na(idx)]]
  }

  manifest <- data.frame(
    trait_prefix = trait_prefix,
    selection = selection,
    collection = combined$collection,
    pathway_id = combined$name,
    full_pathway_name = full_name,
    n_genes = n_genes,
    source_gmt = combined$source_gmt,
    stringsAsFactors = FALSE
  )
} else {
  manifest <- data.frame(
    trait_prefix = character(), selection = character(), collection = character(),
    pathway_id = character(), full_pathway_name = character(),
    n_genes = integer(), source_gmt = character(),
    stringsAsFactors = FALSE
  )
}

write.table(manifest, manifest_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(missing, missing_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Combined selected GMT written:", combined_out, "\n")
cat("Manifest written:", manifest_out, "\n")
cat("Missing term report written:", missing_out, "\n")
cat("Total selected pathways requested:", n_requested_total, "\n")
cat("Total selected pathways matched/written:", n_written_total, "\n")
cat("Total missing selected pathways:", nrow(missing), "\n")

if (partial) {
  cat("\nWARNING: selected pathway matching is incomplete.\n")
  cat("Do NOT use this partial GMT for final pathway/non-pathway scoring.\n")
  cat("The usual cause is either MAGMA truncating long pathway names with '...' or using selected tables created before short-ID GMTs were introduced.\n")
  cat("Recommended fix: run scripts/4a_make-short-id-gmts.R, rerun script 4 using the short-ID Entrez GMTs, then rerun scripts 5 and 6a.\n")
}

if (nrow(combined) == 0) {
  cat("\nWARNING: combined GMT contains zero pathways.\n")
  cat("Check whether selected results exist and whether symbol GMT files are present.\n")
}
