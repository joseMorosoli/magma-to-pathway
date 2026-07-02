#!/bin/bash -l
# ============================================================================
# Script: 6b_create-nonpathway-snp-list.sh
# Purpose: Create SNP lists for selected pathway SNPs and their non-pathway
#          complement in the target genotype data.
#
# Run after:
#   6a_make-selected-pathway-gmt.R
#
# Required input:
#   - selected pathway GMT using gene symbols
#   - target PLINK .bim file
#   - GRCh37 GTF with gene_name values
#
# Output:
#   pathway_union_snps.txt
#   nonpathway_snps.extract.txt
#   nonpathway_snp_definition_metadata.tsv
#
# Notes:
#   - This uses the same 35 kb upstream / 10 kb downstream convention.
#   - Run one trait at a time by editing TRAIT_PREFIX and SELECTION.
# ============================================================================

set +e
set -o pipefail

module -f unload compilers mpi gcc-libs
module load r/recommended

export R_LIBS="/myriadfs/home/ucju659/MyRLibs/"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# Choose one trait.
# TRAIT_PREFIX="HT_EUR_2022"
TRAIT_PREFIX="BMI_EUR_2018"

SELECTION="FDR05"

SELECTED_GMT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_GO_Reactome.symbols.gmt"
OUTDIR="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_nonpathway"

TARGET_BIM="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only.bim"
GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

UPSTREAM_BP="35000"
DOWNSTREAM_BP="10000"

mkdir -p "${OUTDIR}"

CAN_RUN=true

if [[ ! -f "${SELECTED_GMT}" ]]; then
  echo "WARNING: selected GMT not found: ${SELECTED_GMT}" >&2
  CAN_RUN=false
fi

if [[ ! -f "${TARGET_BIM}" ]]; then
  echo "WARNING: target BIM not found: ${TARGET_BIM}" >&2
  CAN_RUN=false
fi

if [[ ! -f "${GTF}" ]]; then
  echo "WARNING: GTF not found: ${GTF}" >&2
  CAN_RUN=false
fi

if [[ "${CAN_RUN}" != true ]]; then
  echo "WARNING: skipping non-pathway SNP list creation because an input is missing." >&2
else
  export SELECTED_GMT OUTDIR TARGET_BIM GTF UPSTREAM_BP DOWNSTREAM_BP TRAIT_PREFIX SELECTION

  Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
})

selected_gmt <- Sys.getenv("SELECTED_GMT")
outdir <- Sys.getenv("OUTDIR")
target_bim <- Sys.getenv("TARGET_BIM")
gtf_path <- Sys.getenv("GTF")
upstream <- as.integer(Sys.getenv("UPSTREAM_BP"))
downstream <- as.integer(Sys.getenv("DOWNSTREAM_BP"))
trait_prefix <- Sys.getenv("TRAIT_PREFIX")
selection <- Sys.getenv("SELECTION")

cat("Selected GMT:", selected_gmt, "\n")
cat("Target BIM:", target_bim, "\n")
cat("GTF:", gtf_path, "\n")
cat("Output folder:", outdir, "\n")

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    x <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(x) >= 3) {
      genes <- unique(toupper(x[-c(1, 2)]))
      genes <- genes[nzchar(genes)]
      out[[x[1]]] <- genes
    }
  }
  out
}

pathways <- read_gmt(selected_gmt)
if (length(pathways) == 0) stop("Selected GMT contained zero pathways: ", selected_gmt)

pathway_genes <- unique(unlist(pathways, use.names = FALSE))
cat("Pathways in GMT:", length(pathways), "\n")
cat("Unique pathway genes in GMT:", length(pathway_genes), "\n")

bim <- fread(target_bim, header = FALSE)
setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
bim <- bim[as.character(CHR) %in% as.character(1:22)]
bim[, CHR := as.character(CHR)]
bim[, BP := as.integer(BP)]
bim[, SNP := as.character(SNP)]

cat("Autosomal SNPs in target BIM:", nrow(bim), "\n")

snpGR <- GRanges(seqnames = paste0("chr", bim$CHR), ranges = IRanges(bim$BP, bim$BP))
snpGR$SNP <- bim$SNP

gtf <- rtracklayer::import(gtf_path)
genes <- gtf[gtf$type == "gene"]

if (length(genes) == 0) stop("No gene features found in GTF.")

if ("gene_name" %in% names(mcols(genes))) {
  gene_names <- genes$gene_name
} else if ("gene_id" %in% names(mcols(genes))) {
  gene_names <- genes$gene_id
} else {
  stop("No gene_name or gene_id found in GTF metadata.")
}

genes$gene_clean <- toupper(as.character(gene_names))

cat("Unique genes in GTF:", length(unique(genes$gene_clean)), "\n")

# Strand-aware upstream/downstream expansion.
geneWin <- genes
st <- as.character(strand(genes))
new_start <- start(genes)
new_end <- end(genes)

plus <- st == "+"
minus <- st == "-"
other <- !(plus | minus)

new_start[plus] <- pmax(1L, start(genes)[plus] - upstream)
new_end[plus] <- end(genes)[plus] + downstream

new_start[minus] <- pmax(1L, start(genes)[minus] - downstream)
new_end[minus] <- end(genes)[minus] + upstream

new_start[other] <- pmax(1L, start(genes)[other] - upstream)
new_end[other] <- end(genes)[other] + downstream

ranges(geneWin) <- IRanges(new_start, new_end)
names(geneWin) <- genes$gene_clean

geneWin <- keepStandardChromosomes(geneWin, pruning.mode = "coarse")
snpGR <- keepStandardChromosomes(snpGR, pruning.mode = "coarse")
seqlevelsStyle(snpGR) <- seqlevelsStyle(geneWin)[1]

selected_gene_windows <- geneWin[names(geneWin) %in% pathway_genes]
cat("Pathway genes found in GTF:", length(unique(names(selected_gene_windows))), "\n")

if (length(selected_gene_windows) == 0) {
  stop("No selected GMT genes matched the GTF. Check gene identifier type: symbols vs Entrez/Ensembl.")
}

ov <- findOverlaps(snpGR, selected_gene_windows, ignore.strand = TRUE)
pathway_snps <- unique(snpGR$SNP[queryHits(ov)])
nonpathway_snps <- setdiff(bim$SNP, pathway_snps)

cat("SNPs mapped to selected pathway genes:", length(pathway_snps), "\n")
cat("SNPs outside selected pathway genes:", length(nonpathway_snps), "\n")

pathway_file <- file.path(outdir, "pathway_union_snps.txt")
nonpath_file <- file.path(outdir, "nonpathway_snps.extract.txt")
meta_file <- file.path(outdir, "nonpathway_snp_definition_metadata.tsv")
genes_file <- file.path(outdir, "selected_pathway_genes_from_gmt.txt")
genes_in_gtf_file <- file.path(outdir, "selected_pathway_genes_found_in_gtf.txt")

fwrite(data.table(SNP = pathway_snps), pathway_file, col.names = FALSE)
fwrite(data.table(SNP = nonpathway_snps), nonpath_file, col.names = FALSE)
fwrite(data.table(gene = sort(pathway_genes)), genes_file, col.names = FALSE)
fwrite(data.table(gene = sort(unique(names(selected_gene_windows)))), genes_in_gtf_file, col.names = FALSE)

meta <- data.table(
  trait_prefix = trait_prefix,
  selection = selection,
  selected_gmt = selected_gmt,
  target_bim = target_bim,
  gtf = gtf_path,
  upstream_bp = upstream,
  downstream_bp = downstream,
  n_pathways = length(pathways),
  n_unique_pathway_genes_gmt = length(pathway_genes),
  n_pathway_genes_found_in_gtf = length(unique(names(selected_gene_windows))),
  n_target_autosomal_snps = nrow(bim),
  n_pathway_union_snps = length(pathway_snps),
  n_nonpathway_snps = length(nonpathway_snps)
)

fwrite(meta, meta_file, sep = "\t")

cat("Written:\n")
cat(pathway_file, "\n")
cat(nonpath_file, "\n")
cat(meta_file, "\n")
cat(genes_file, "\n")
cat(genes_in_gtf_file, "\n")
RSCRIPT
fi
