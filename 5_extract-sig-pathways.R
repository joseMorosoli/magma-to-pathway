#!/usr/bin/env Rscript
# ============================================================================
# Script: 5_extract-sig-pathways.R
# Purpose: Add multiple-testing correction and extract significant/top MAGMA
#          pathway results.
#
# Intended use:
#   Run after script 4 has created .gsa.out files.
#
# Required input:
#   One or more MAGMA .gsa.out files with a P column.
#
# Output for each input:
#   <prefix>_all_with_FDR.tsv
#   <prefix>_nominal_P05.tsv
#   <prefix>_FDR05.tsv
#   <prefix>_Bonferroni05.tsv
#   <prefix>_top50.tsv
#
# Notes:
#   - GO and Reactome are corrected separately because they are different
#     pathway collections / multiple-testing families.
#   - Interpret pathway results cautiously: significance depends on pathway
#     definition, SNP-to-gene mapping, LD, ancestry/reference panel and GWAS power.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

RESULTS_DIR <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/results"

INPUTS <- data.table::data.table(
  collection = c("GO", "Reactome"),
  gsa_file = c(
    file.path(RESULTS_DIR, "HT_EUR_2022_GO.gsa.out"),
    file.path(RESULTS_DIR, "HT_EUR_2022_Reactome.gsa.out")
  ),
  out_prefix = c(
    file.path(RESULTS_DIR, "HT_EUR_2022_GO"),
    file.path(RESULTS_DIR, "HT_EUR_2022_Reactome")
  )
)

# ---------------------------------------------------------------------------
# Function
# ---------------------------------------------------------------------------

process_gsa <- function(gsa_file, out_prefix, collection) {
  if (!file.exists(gsa_file)) {
    stop("MAGMA .gsa.out file not found for ", collection, ": ", gsa_file, call. = FALSE)
  }

  gsa <- fread(gsa_file)

  if (!"P" %in% names(gsa)) {
    stop("Column 'P' not found in ", gsa_file,
         ". Observed columns: ", paste(names(gsa), collapse = ", "),
         call. = FALSE)
  }

  gsa <- gsa[!is.na(P) & P >= 0 & P <= 1]
  gsa[, Bonferroni := p.adjust(P, method = "bonferroni")]
  gsa[, FDR := p.adjust(P, method = "BH")]
  gsa[, minus_log10_p := -log10(pmax(P, .Machine$double.xmin))]
  gsa[, minus_log10_fdr := -log10(pmax(FDR, .Machine$double.xmin))]
  setorder(gsa, FDR, P)

  all_file <- paste0(out_prefix, "_all_with_FDR.tsv")
  nom_file <- paste0(out_prefix, "_nominal_P05.tsv")
  fdr_file <- paste0(out_prefix, "_FDR05.tsv")
  bonf_file <- paste0(out_prefix, "_Bonferroni05.tsv")
  top_file <- paste0(out_prefix, "_top50.tsv")

  fwrite(gsa, all_file, sep = "\t")
  fwrite(gsa[P < 0.05], nom_file, sep = "\t")
  fwrite(gsa[FDR < 0.05], fdr_file, sep = "\t")
  fwrite(gsa[Bonferroni < 0.05], bonf_file, sep = "\t")
  fwrite(head(gsa, 50), top_file, sep = "\t")

  cat("\nCollection: ", collection, "\n", sep = "")
  cat("Input: ", gsa_file, "\n", sep = "")
  cat("Total pathways tested: ", nrow(gsa), "\n", sep = "")
  cat("Nominal P < 0.05: ", nrow(gsa[P < 0.05]), "\n", sep = "")
  cat("FDR < 0.05: ", nrow(gsa[FDR < 0.05]), "\n", sep = "")
  cat("Bonferroni < 0.05: ", nrow(gsa[Bonferroni < 0.05]), "\n", sep = "")
  cat("Files written:\n")
  cat("  ", all_file, "\n", sep = "")
  cat("  ", nom_file, "\n", sep = "")
  cat("  ", fdr_file, "\n", sep = "")
  cat("  ", bonf_file, "\n", sep = "")
  cat("  ", top_file, "\n", sep = "")

  cols_to_print <- intersect(c("VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR", "Bonferroni"), names(gsa))
  print(head(gsa[, ..cols_to_print], 20))

  invisible(gsa)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

results <- vector("list", nrow(INPUTS))
for (i in seq_len(nrow(INPUTS))) {
  results[[i]] <- process_gsa(
    gsa_file = INPUTS$gsa_file[i],
    out_prefix = INPUTS$out_prefix[i],
    collection = INPUTS$collection[i]
  )
}
