
#!/usr/bin/env Rscript

# ============================================================
# Script: 5_extract-sig-pathways.R
# Purpose:
#   Add Bonferroni and FDR correction to MAGMA pathway results.
#
# Run after:
#   4_run-MAGMA-pathways.sh
#
# Required inputs:
#   <TRAIT>_GO.gsa.out
#   <TRAIT>_Reactome.gsa.out
#
# Outputs:
#   <TRAIT>_<collection>_all_with_FDR.csv
#   <TRAIT>_<collection>_nominal_P05.csv
#   <TRAIT>_<collection>_FDR05.csv
#   <TRAIT>_<collection>_Bonferroni05.csv
#   <TRAIT>_<collection>_top50.csv
#
# Notes:
#   GO and Reactome are corrected separately.
#   This script uses base R only.
# ============================================================


module -f unload compilers mpi gcc-libs
module load r/recommended
R --vanilla

# ------------------------------------------------------------
# 1. Paths and trait settings
# ------------------------------------------------------------

results_dir <- "/myriadfs/home/ucju659/SOFTWARE/MAGMA/results"

# Choose trait by uncommenting one block.

# Height
#trait_prefix <- "HT_EUR_2022"

# BMI
#trait_prefix <- "BMI_EUR_2018"

# F4
trait_prefix <- "F4_2025"


# Choose which pathway collections to process.

run_go <- TRUE
run_reactome <- TRUE


# ------------------------------------------------------------
# 2. GO analysis
# ------------------------------------------------------------

if (run_go) {
  
  collection <- "GO"
  
  gsa_file <- file.path(results_dir, paste0(trait_prefix, "_GO.gsa.out"))
  
  out_all <- file.path(results_dir, paste0(trait_prefix, "_GO_all_with_FDR.csv"))
  out_nom <- file.path(results_dir, paste0(trait_prefix, "_GO_nominal_P05.csv"))
  out_fdr <- file.path(results_dir, paste0(trait_prefix, "_GO_FDR05.csv"))
  out_bf  <- file.path(results_dir, paste0(trait_prefix, "_GO_Bonferroni05.csv"))
  out_top <- file.path(results_dir, paste0(trait_prefix, "_GO_top50.csv"))
  
  cat("\n============================================================\n")
  cat("Processing:", collection, "\n")
  cat("Input:", gsa_file, "\n")
  
  if (!file.exists(gsa_file)) {
    cat("WARNING: file not found. Skipping GO:\n")
    cat(gsa_file, "\n")
  } else {
    
    gsa <- read.table(
      gsa_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("Columns found:\n")
    print(names(gsa))
    
    if (!"P" %in% names(gsa)) {
      cat("WARNING: column P not found. Skipping GO.\n")
    } else {
      
      gsa <- gsa[!is.na(gsa$P) & gsa$P >= 0 & gsa$P <= 1, ]
      
      gsa$Bonferroni <- p.adjust(gsa$P, method = "bonferroni")
      gsa$FDR <- p.adjust(gsa$P, method = "BH")
      gsa$minus_log10_p <- -log10(pmax(gsa$P, .Machine$double.xmin))
      gsa$minus_log10_fdr <- -log10(pmax(gsa$FDR, .Machine$double.xmin))
      
      gsa <- gsa[order(gsa$FDR, gsa$P), ]
      
      sig_nom <- gsa[gsa$P < 0.05, ]
      sig_fdr <- gsa[gsa$FDR < 0.05, ]
      sig_bf  <- gsa[gsa$Bonferroni < 0.05, ]
      top50   <- head(gsa, 50)
      
      write.csv(gsa, out_all, row.names = FALSE)
      write.csv(sig_nom, out_nom, row.names = FALSE)
      write.csv(sig_fdr, out_fdr, row.names = FALSE)
      write.csv(sig_bf, out_bf, row.names = FALSE)
      write.csv(top50, out_top, row.names = FALSE)
      
      cat("Total pathways tested:", nrow(gsa), "\n")
      cat("Nominal P < 0.05:", nrow(sig_nom), "\n")
      cat("FDR < 0.05:", nrow(sig_fdr), "\n")
      cat("Bonferroni < 0.05:", nrow(sig_bf), "\n")
      
      cat("Files written:\n")
      cat(out_all, "\n")
      cat(out_nom, "\n")
      cat(out_fdr, "\n")
      cat(out_bf, "\n")
      cat(out_top, "\n")
      
      cols_to_print <- c("VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR", "Bonferroni")
      cols_to_print <- cols_to_print[cols_to_print %in% names(gsa)]
      
      cat("\nTop pathways:\n")
      print(head(gsa[, cols_to_print, drop = FALSE], 20))
    }
  }
}


# ------------------------------------------------------------
# 3. Reactome analysis
# ------------------------------------------------------------

if (run_reactome) {
  
  collection <- "Reactome"
  
  gsa_file <- file.path(results_dir, paste0(trait_prefix, "_Reactome.gsa.out"))
  
  out_all <- file.path(results_dir, paste0(trait_prefix, "_Reactome_all_with_FDR.csv"))
  out_nom <- file.path(results_dir, paste0(trait_prefix, "_Reactome_nominal_P05.csv"))
  out_fdr <- file.path(results_dir, paste0(trait_prefix, "_Reactome_FDR05.csv"))
  out_bf  <- file.path(results_dir, paste0(trait_prefix, "_Reactome_Bonferroni05.csv"))
  out_top <- file.path(results_dir, paste0(trait_prefix, "_Reactome_top50.csv"))
  
  cat("\n============================================================\n")
  cat("Processing:", collection, "\n")
  cat("Input:", gsa_file, "\n")
  
  if (!file.exists(gsa_file)) {
    cat("WARNING: file not found. Skipping Reactome:\n")
    cat(gsa_file, "\n")
  } else {
    
    gsa <- read.table(
      gsa_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    cat("Columns found:\n")
    print(names(gsa))
    
    if (!"P" %in% names(gsa)) {
      cat("WARNING: column P not found. Skipping Reactome.\n")
    } else {
      
      gsa <- gsa[!is.na(gsa$P) & gsa$P >= 0 & gsa$P <= 1, ]
      
      gsa$Bonferroni <- p.adjust(gsa$P, method = "bonferroni")
      gsa$FDR <- p.adjust(gsa$P, method = "BH")
      gsa$minus_log10_p <- -log10(pmax(gsa$P, .Machine$double.xmin))
      gsa$minus_log10_fdr <- -log10(pmax(gsa$FDR, .Machine$double.xmin))
      
      gsa <- gsa[order(gsa$FDR, gsa$P), ]
      
      sig_nom <- gsa[gsa$P < 0.05, ]
      sig_fdr <- gsa[gsa$FDR < 0.05, ]
      sig_bf  <- gsa[gsa$Bonferroni < 0.05, ]
      top50   <- head(gsa, 50)
      
      write.csv(gsa, out_all, row.names = FALSE)
      write.csv(sig_nom, out_nom, row.names = FALSE)
      write.csv(sig_fdr, out_fdr, row.names = FALSE)
      write.csv(sig_bf, out_bf, row.names = FALSE)
      write.csv(top50, out_top, row.names = FALSE)
      
      cat("Total pathways tested:", nrow(gsa), "\n")
      cat("Nominal P < 0.05:", nrow(sig_nom), "\n")
      cat("FDR < 0.05:", nrow(sig_fdr), "\n")
      cat("Bonferroni < 0.05:", nrow(sig_bf), "\n")
      
      cat("Files written:\n")
      cat(out_all, "\n")
      cat(out_nom, "\n")
      cat(out_fdr, "\n")
      cat(out_bf, "\n")
      cat(out_top, "\n")
      
      cols_to_print <- c("VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "FDR", "Bonferroni")
      cols_to_print <- cols_to_print[cols_to_print %in% names(gsa)]
      
      cat("\nTop pathways:\n")
      print(head(gsa[, cols_to_print, drop = FALSE], 20))
    }
  }
}

cat("\nDone.\n")
