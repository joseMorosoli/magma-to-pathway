#!/usr/bin/env Rscript

# ============================================================
# Script: 5b_check-selected-pathways.R
# Purpose:
#   Check which MAGMA pathways are selected for each trait at
#   different multiple-testing thresholds.
#
# Run after:
#   5_extract-sig-pathways.R
#
# Input expected for each trait/collection:
#   <trait>_<collection>_all_with_FDR.csv
#
# Outputs:
#   selected_pathways_counts_by_trait_level.csv
#   selected_pathways_by_trait_level.csv
#
# Notes:
#   This is a diagnostic script only. It does not create GMT files.
#   It uses base R only.
# ============================================================

module -f unload compilers mpi gcc-libs
module load r/recommended
R --vanilla

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

results_dir <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/results"

# Traits to check.
traits <- c(
  "HT_EUR_2022",
  "BMI_EUR_2018",
  "F4_2025"
)

# Pathway collections to check.
collections <- c(
  "GO",
  "Reactome"
)

# Thresholds.
nominal_alpha <- 0.05
fdr_alpha <- 0.05
bonferroni_alpha <- 0.05

top_n <- 50

# Output files.
out_counts <- file.path(results_dir, "selected_pathways_counts_by_trait_level.csv")
out_long   <- file.path(results_dir, "selected_pathways_by_trait_level.csv")


# ------------------------------------------------------------
# 2. Main loop
# ------------------------------------------------------------

all_counts <- data.frame()
all_selected <- data.frame()

cat("Checking selected MAGMA pathways\n")
cat("Results folder:", results_dir, "\n\n")

for (trait in traits) {
  for (collection in collections) {

    file_csv <- file.path(results_dir, paste0(trait, "_", collection, "_all_with_FDR.csv"))
    file_tsv <- file.path(results_dir, paste0(trait, "_", collection, "_all_with_FDR.tsv"))

    if (file.exists(file_csv)) {
      gsa_file <- file_csv
      gsa <- read.csv(gsa_file, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (file.exists(file_tsv)) {
      gsa_file <- file_tsv
      gsa <- read.delim(gsa_file, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      cat("WARNING: no all_with_FDR file found for", trait, collection, "\n")
      next
    }

    cat("============================================================\n")
    cat("Trait:", trait, "\n")
    cat("Collection:", collection, "\n")
    cat("Input:", gsa_file, "\n")

    if (!"VARIABLE" %in% names(gsa)) {
      cat("WARNING: VARIABLE column not found. Skipping.\n\n")
      next
    }

    if (!"P" %in% names(gsa)) {
      cat("WARNING: P column not found. Skipping.\n\n")
      next
    }

    if (!"FDR" %in% names(gsa)) {
      cat("WARNING: FDR column not found. Recalculating FDR from P.\n")
      gsa$FDR <- p.adjust(gsa$P, method = "BH")
    }

    if (!"Bonferroni" %in% names(gsa)) {
      if ("BF" %in% names(gsa)) {
        gsa$Bonferroni <- gsa$BF
      } else {
        cat("WARNING: Bonferroni column not found. Recalculating Bonferroni from P.\n")
        gsa$Bonferroni <- p.adjust(gsa$P, method = "bonferroni")
      }
    }

    gsa <- gsa[!is.na(gsa$P), ]
    gsa <- gsa[order(gsa$P), ]

    selected_nominal <- gsa[gsa$P < nominal_alpha, ]
    selected_fdr <- gsa[gsa$FDR < fdr_alpha, ]
    selected_bonferroni <- gsa[gsa$Bonferroni < bonferroni_alpha, ]
    selected_top <- head(gsa, top_n)

    cat("Total pathways tested:", nrow(gsa), "\n")
    cat("Nominal P <", nominal_alpha, ":", nrow(selected_nominal), "\n")
    cat("FDR <", fdr_alpha, ":", nrow(selected_fdr), "\n")
    cat("Bonferroni <", bonferroni_alpha, ":", nrow(selected_bonferroni), "\n")
    cat("Top", top_n, ":", nrow(selected_top), "\n")

    counts <- data.frame(
      trait = trait,
      collection = collection,
      n_tested = nrow(gsa),
      n_nominal_P05 = nrow(selected_nominal),
      n_FDR05 = nrow(selected_fdr),
      n_Bonferroni05 = nrow(selected_bonferroni),
      n_top50 = nrow(selected_top),
      stringsAsFactors = FALSE
    )

    all_counts <- rbind(all_counts, counts)

    keep_cols <- c("VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR", "Bonferroni")
    keep_cols <- keep_cols[keep_cols %in% names(gsa)]

    if (nrow(selected_nominal) > 0) {
      tmp <- selected_nominal[, keep_cols, drop = FALSE]
      tmp$trait <- trait
      tmp$collection <- collection
      tmp$selection_level <- "nominal_P05"
      tmp <- tmp[, c("trait", "collection", "selection_level", keep_cols), drop = FALSE]
      all_selected <- rbind(all_selected, tmp)
    }

    if (nrow(selected_fdr) > 0) {
      tmp <- selected_fdr[, keep_cols, drop = FALSE]
      tmp$trait <- trait
      tmp$collection <- collection
      tmp$selection_level <- "FDR05"
      tmp <- tmp[, c("trait", "collection", "selection_level", keep_cols), drop = FALSE]
      all_selected <- rbind(all_selected, tmp)
    }

    if (nrow(selected_bonferroni) > 0) {
      tmp <- selected_bonferroni[, keep_cols, drop = FALSE]
      tmp$trait <- trait
      tmp$collection <- collection
      tmp$selection_level <- "Bonferroni05"
      tmp <- tmp[, c("trait", "collection", "selection_level", keep_cols), drop = FALSE]
      all_selected <- rbind(all_selected, tmp)
    }

    if (nrow(selected_top) > 0) {
      tmp <- selected_top[, keep_cols, drop = FALSE]
      tmp$trait <- trait
      tmp$collection <- collection
      tmp$selection_level <- paste0("top", top_n)
      tmp <- tmp[, c("trait", "collection", "selection_level", keep_cols), drop = FALSE]
      all_selected <- rbind(all_selected, tmp)
    }

    cat("\nTop 20 by P-value:\n")
    print(head(gsa[, c("VARIABLE", "NGENES", "P", "FDR", "Bonferroni"), drop = FALSE], 20))
    cat("\n")
  }
}


# ------------------------------------------------------------
# 3. Save outputs
# ------------------------------------------------------------

write.csv(all_counts, out_counts, row.names = FALSE)
write.csv(all_selected, out_long, row.names = FALSE)

cat("============================================================\n")
cat("Summary counts written to:\n")
cat(out_counts, "\n")
cat("\nSelected pathway table written to:\n")
cat(out_long, "\n")
cat("\nDone.\n")
