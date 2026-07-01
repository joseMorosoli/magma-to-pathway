#!/usr/bin/env Rscript
# ============================================================================
# Script: 4_run-MAGMA-pathways.R
# Purpose: Run MAGMA gene-set/pathway analysis from an existing .genes.raw file.
#
# Intended use:
#   Run locally or on Myriad after script 3 has successfully created:
#     <GENE_PREFIX>.genes.raw
#
# Required input:
#   - MAGMA executable
#   - trait-specific MAGMA .genes.raw file
#   - GO and/or Reactome GMT files
#
# Output:
#   <OUT_PREFIX>.gsa.out files for each pathway collection.
#
# Notes:
#   - This script does not rerun gene-level analysis.
#   - GO and Reactome are run as separate gene-set collections.
#   - Multiple-testing correction is handled in script 5.
# ============================================================================

# ---------------------------------------------------------------------------
# USER SETTINGS: edit these lines for each trait/environment.
# ---------------------------------------------------------------------------

PROJECT_DIR <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA"
MAGMA_EXE   <- file.path(PROJECT_DIR, "v1.10/magma")

# Example: height. Change only these output prefixes for each trait.
GENE_RAW    <- file.path(PROJECT_DIR, "results/HT_EUR_2022_genes.genes.raw")
GO_GMT      <- file.path(PROJECT_DIR, "pathways/c5.go.v2026.1.Hs.entrez.gmt")
REACTOME_GMT <- file.path(PROJECT_DIR, "pathways/c2.cp.reactome.v2026.1.Hs.entrez.gmt")

GO_OUT      <- file.path(PROJECT_DIR, "results/HT_EUR_2022_GO")
REACTOME_OUT <- file.path(PROJECT_DIR, "results/HT_EUR_2022_Reactome")

RUN_GO <- TRUE
RUN_REACTOME <- TRUE
OVERWRITE <- TRUE

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

run_magma_set <- function(gene_raw, gmt, out_prefix, collection_name) {
  stop_if_missing(gene_raw, "MAGMA .genes.raw file")
  stop_if_missing(gmt, paste0(collection_name, " GMT file"))

  if (OVERWRITE) {
    unlink(paste0(out_prefix, ".*"))
  }

  message("Running MAGMA pathway analysis: ", collection_name)
  message("Gene results: ", gene_raw)
  message("GMT:          ", gmt)
  message("Out prefix:   ", out_prefix)

  status <- system2(
    MAGMA_EXE,
    args = c(
      "--gene-results", gene_raw,
      "--set-annot", gmt,
      "--out", out_prefix
    )
  )

  if (!identical(status, 0L)) {
    stop("MAGMA failed for ", collection_name, " with exit status ", status, call. = FALSE)
  }

  expected <- paste0(out_prefix, ".gsa.out")
  if (!file.exists(expected) || file.info(expected)$size == 0) {
    stop("Expected output was not created or is empty: ", expected, call. = FALSE)
  }

  message("Created: ", expected)
  invisible(expected)
}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

stop_if_missing(MAGMA_EXE, "MAGMA executable")
stop_if_missing(GENE_RAW, "MAGMA .genes.raw file")

dir.create(dirname(GO_OUT), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REACTOME_OUT), recursive = TRUE, showWarnings = FALSE)

message("MAGMA pathway analysis started: ", Sys.time())

# ---------------------------------------------------------------------------
# Run pathway analyses
# ---------------------------------------------------------------------------

outputs <- character()

if (RUN_GO) {
  outputs <- c(outputs, run_magma_set(GENE_RAW, GO_GMT, GO_OUT, "GO"))
}

if (RUN_REACTOME) {
  outputs <- c(outputs, run_magma_set(GENE_RAW, REACTOME_GMT, REACTOME_OUT, "Reactome"))
}

message("MAGMA pathway analysis completed: ", Sys.time())
message("Outputs:")
print(outputs)
