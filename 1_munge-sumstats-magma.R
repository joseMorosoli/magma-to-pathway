#!/usr/bin/env Rscript
# ============================================================================
# Script: 1_munge-sumstats-magma.R
# Purpose: Create a MAGMA-ready summary-statistics file for one trait.
#
# Intended use:
#   Run locally or on Myriad after creating an LDpred2-ready summary-statistics
#   file with script 0. This script is deliberately one-trait-at-a-time.
#
# Required input:
#   1. LDpred2-ready summary statistics with at least:
#        CHR BP A2 A1 N BETA SE P
#      Optional but recommended:
#        MAF INFO
#   2. MAGMA/1000G reference .bim file whose SNP IDs will be used in the final
#      MAGMA input.
#
# Output:
#   A tab-delimited MAGMA-ready file containing at least:
#        SNP P N
#   Additional columns are retained for auditability.
#
# Reproducibility notes:
#   - The SNP column must match the second column of the reference .bim file.
#   - Matching by CHR:BP assumes the summary statistics and .bim file use the
#     same genome build. For this project, the MAGMA reference is GRCh37.
#   - Do not commit munged summary statistics to GitHub.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# USER SETTINGS: edit these lines for each GWAS/trait.
# ---------------------------------------------------------------------------

# Example local Windows paths. Change these to your actual local or Myriad paths.
SUMSTATS_FILE <- "C:/Users/Jose Morosoli/Documents/UCL/SUMSTATS/ldpred2_ready/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz"
BIM_FILE      <- "C:/MAGMA/g1000_eur/g1000_eur.bim"
OUT_FILE      <- "C:/MAGMA/munged/HT_EUR_2022_MAGMA.txt"

TRAIT_LABEL   <- "height_2022"
GENOME_BUILD  <- "GRCh37"
GWAS_SOURCE   <- "GIANT height 2022; replace with exact citation/download source"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

required_cols <- c("CHR", "BP", "A2", "A1", "N", "BETA", "SE", "P")

# ---------------------------------------------------------------------------
# Checks and read data
# ---------------------------------------------------------------------------

stop_if_missing(SUMSTATS_FILE, "Summary-statistics file")
stop_if_missing(BIM_FILE, "Reference .bim file")
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

cat("Trait label: ", TRAIT_LABEL, "\n", sep = "")
cat("GWAS source: ", GWAS_SOURCE, "\n", sep = "")
cat("Genome build assumed for CHR/BP matching: ", GENOME_BUILD, "\n", sep = "")
cat("Summary statistics: ", SUMSTATS_FILE, "\n", sep = "")
cat("Reference BIM: ", BIM_FILE, "\n", sep = "")
cat("Output: ", OUT_FILE, "\n\n", sep = "")

ss <- fread(SUMSTATS_FILE)
setnames(ss, names(ss), gsub("\\r", "", names(ss)))

missing_cols <- setdiff(required_cols, names(ss))
if (length(missing_cols) > 0) {
  stop("Missing required columns in summary statistics: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

bim <- fread(
  BIM_FILE,
  col.names = c("CHR", "SNP", "CM", "BP", "A1_BIM", "A2_BIM"),
  colClasses = list(character = c("SNP", "A1_BIM", "A2_BIM"))
)

# Keep one BIM SNP per physical position to avoid duplicate CHR:BP ambiguity.
n_bim_before <- nrow(bim)
bim <- bim[!duplicated(paste(CHR, BP))]
cat("Reference BIM variants: ", n_bim_before, "\n", sep = "")
cat("Reference BIM variants after removing duplicated CHR:BP: ", nrow(bim), "\n", sep = "")

# ---------------------------------------------------------------------------
# Basic QC
# ---------------------------------------------------------------------------

ss <- ss[CHR %in% 1:22]
ss <- ss[!is.na(BP) & BP > 0]
ss <- ss[!is.na(P) & P > 0 & P <= 1]
ss <- ss[!is.na(N) & N > 0]
ss <- ss[!is.na(BETA) & is.finite(BETA)]
ss <- ss[!is.na(SE) & is.finite(SE) & SE > 0]
ss[, A1 := toupper(A1)]
ss[, A2 := toupper(A2)]
ss <- ss[A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T")]
ss <- ss[A1 != A2]

n_ss_after_qc <- nrow(ss)
cat("Summary-statistic rows after basic QC: ", n_ss_after_qc, "\n", sep = "")

# ---------------------------------------------------------------------------
# Map to reference .bim SNP IDs by CHR:BP
# ---------------------------------------------------------------------------

ss[, CHR := as.integer(CHR)]
ss[, BP := as.integer(BP)]
bim[, CHR := as.integer(CHR)]
bim[, BP := as.integer(BP)]

setkey(ss, CHR, BP)
setkey(bim, CHR, BP)

magma <- merge(ss, bim[, .(CHR, BP, SNP, A1_BIM, A2_BIM)], by = c("CHR", "BP"), all = FALSE)

cat("Rows matched to reference BIM by CHR:BP: ", nrow(magma), "\n", sep = "")
if (nrow(magma) == 0) {
  stop("No variants matched the reference .bim file. Check genome build and CHR/BP columns.", call. = FALSE)
}

overlap_rate <- nrow(magma) / n_ss_after_qc
cat("Proportion of QC-passing summary-statistic rows matched: ",
    round(overlap_rate, 4), "\n", sep = "")

if (overlap_rate < 0.5) {
  warning("Less than 50% of QC-passing variants matched the reference .bim. Check genome build/SNP map.")
}

# Create audit columns similar to the user's current MAGMA-ready files.
magma[, ID := paste(CHR, BP, A1, A2, sep = ":")]
magma[, RSID := SNP]
if (!"MAF" %in% names(magma)) magma[, MAF := NA_real_]
if (!"INFO" %in% names(magma)) magma[, INFO := NA_real_]

out <- magma[, .(
  SNP,
  CHR,
  BP,
  A1,
  A2,
  ID,
  RSID,
  FRQ = MAF,
  BETA,
  SE,
  P,
  N,
  INFO
)]

# MAGMA can tolerate extra columns, but SNP/P/N are the required ones.
required_magma <- c("SNP", "P", "N")
stopifnot(all(required_magma %in% names(out)))

if (anyDuplicated(out$SNP)) {
  dup_n <- sum(duplicated(out$SNP))
  warning("Duplicated SNP IDs found after mapping: ", dup_n,
          ". Keeping first occurrence per SNP.")
  out <- out[!duplicated(SNP)]
}

fwrite(out, OUT_FILE, sep = "\t", quote = FALSE, na = "NA")

cat("\nCreated MAGMA-ready file: ", OUT_FILE, "\n", sep = "")
cat("Rows written: ", nrow(out), "\n", sep = "")
cat("Header:\n")
print(names(out))
cat("Preview:\n")
print(head(out))
