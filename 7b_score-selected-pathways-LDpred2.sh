#!/bin/bash -l
# ============================================================================
# Script: 7b_score-selected-pathways-LDpred2.sh
# Purpose: Calculate one LDpred2 posterior-beta score per selected pathway.
#
# Run after:
#   6a_make-selected-pathway-gmt.R
#
# Required input:
#   - selected pathway GMT using gene symbols
#   - target bigsnpr .rds
#   - LDpred2 posterior beta file with columns chr pos a0 a1 beta
#   - GRCh37 GTF file
#
# Output:
#   <OUT_PREFIX_NAME>.tsv.gz
#   <OUT_PREFIX_NAME>.metadata.tsv
#   <OUT_PREFIX_NAME>.snps_used.tsv.gz
#   <OUT_PREFIX_NAME>.rds
# ============================================================================

set +e
set -o pipefail

module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0
#module load r/recommended

unset R_LIBS
export R_LIBS_USER="/myriadfs/home/ucju659/MyRLibs/R-4.5.1"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export THREADS="${NSLOTS:-2}"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# Choose one trait.
#TRAIT_PREFIX="HT_EUR_2022"
#TRAIT_PREFIX="BMI_EUR_2018"
TRAIT_PREFIX="F4_2025"

# Common values: "FDR05", "Bonferroni05", "nominal_P05", "top50"
SELECTION="Bonferroni05"

export SELECTED_GMT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_GO_Reactome.symbols.gmt"
export BIGSNP_RDS="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only.rds"
export GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

# Edit this for height/BMI.
# Expected columns: chr pos a0 a1 beta, with or without a header.
#export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz_17June2026/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz_final_beta_auto.txt"
#export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz_17June2026/GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz_final_beta_auto.txt"
export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/F4_Internalizing_2025.ldpred2.tsv.gz_24April2026/F4_Internalizing_2025.ldpred2.tsv.gz_final_beta_auto.txt"

export OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/${TRAIT_PREFIX}/selected_pathways_ldpred2"
export OUT_PREFIX_NAME="${TRAIT_PREFIX}_selected_${SELECTION}_LDpred2_pathway_scores"

export UPSTREAM_BP="35000"
export DOWNSTREAM_BP="10000"

# Optional FID/IID keep file. Leave empty for all individuals.
export IDS_KEEP=""

mkdir -p "${OUTDIR}"

CAN_RUN=true
for file in "${SELECTED_GMT}" "${BIGSNP_RDS}" "${GTF}" "${RAW_BETAS}"; do
  if [[ ! -f "${file}" ]]; then
    echo "WARNING: required file not found: ${file}" >&2
    CAN_RUN=false
  fi
done

if [[ "${CAN_RUN}" != true ]]; then
  echo "WARNING: skipping LDpred2 selected pathway scores because an input is missing." >&2
fi

R --vanilla

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(bigsnpr)
  library(bigstatsr)
})

selected_gmt <- Sys.getenv("SELECTED_GMT")
bigsnp <- Sys.getenv("BIGSNP_RDS")
gtf_path <- Sys.getenv("GTF")
raw_betas <- Sys.getenv("RAW_BETAS")
outdir <- Sys.getenv("OUTDIR")
out_prefix_name <- Sys.getenv("OUT_PREFIX_NAME")
upstream <- as.integer(Sys.getenv("UPSTREAM_BP"))
downstream <- as.integer(Sys.getenv("DOWNSTREAM_BP"))
threads <- as.integer(Sys.getenv("THREADS"))
ids_keep <- Sys.getenv("IDS_KEEP")

normalise_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT", "M")] <- "25"
  x
}

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

cat("Selected GMT:", selected_gmt, "\n")
cat("bigSNP:", bigsnp, "\n")
cat("GTF:", gtf_path, "\n")
cat("LDpred2 betas:", raw_betas, "\n")
cat("Output folder:", outdir, "\n")

pathways <- read_gmt(selected_gmt)
if (length(pathways) == 0) stop("Selected GMT contained zero pathways.")

cat("Selected pathways:", length(pathways), "\n")

obj.bigSNP <- snp_attach(bigsnp)
G <- obj.bigSNP$genotypes
map <- as.data.table(obj.bigSNP$map)

names(map)[1:6] <- c("chr", "rsid", "dist", "pos", "a1", "a0")
map[, chr := as.integer(normalise_chr(chr))]
map[, pos := as.integer(pos)]
map[, rsid := as.character(rsid)]
map[, a1 := toupper(as.character(a1))]
map[, a0 := toupper(as.character(a0))]

map_small <- map[, .(snp_ind = .I, chr, pos, rsid, map_a1 = a1, map_a0 = a0)]
cat("Target SNPs in bigSNP:", nrow(map_small), "\n")

snpGR <- GRanges(seqnames = paste0("chr", map_small$chr), ranges = IRanges(map_small$pos, map_small$pos))
snpGR$SNP <- map_small$rsid

gtf <- rtracklayer::import(gtf_path)
genes <- gtf[gtf$type == "gene"]

if ("gene_name" %in% names(mcols(genes))) {
  gene_names <- genes$gene_name
} else if ("gene_id" %in% names(mcols(genes))) {
  gene_names <- genes$gene_id
} else {
  stop("No gene_name or gene_id found in GTF metadata.")
}

genes$gene_clean <- toupper(as.character(gene_names))

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

cat("Loading LDpred2 beta weights.\n")
first_line <- system2("bash", c("-lc", paste("zcat -f", shQuote(raw_betas), "| head -n 1")), stdout = TRUE)
first_token <- strsplit(trimws(first_line[1]), "\\s+")[[1]][1]
has_header <- grepl("[A-Za-z]", first_token)

betas <- fread(cmd = paste("zcat -f", shQuote(raw_betas)), header = has_header)
betas <- betas[, 1:5]
setnames(betas, c("chr", "pos", "a0", "a1", "weight"))
betas[, chr := as.integer(normalise_chr(chr))]
betas[, pos := as.integer(pos)]
betas[, a0 := toupper(as.character(a0))]
betas[, a1 := toupper(as.character(a1))]
betas[, weight := as.numeric(weight)]
betas <- betas[!is.na(chr) & !is.na(pos) & !is.na(weight)]

map_coord <- map_small[!duplicated(paste(chr, pos)) & !duplicated(paste(chr, pos), fromLast = TRUE)]
beta_coord <- betas[!duplicated(paste(chr, pos)) & !duplicated(paste(chr, pos), fromLast = TRUE)]

keep <- merge(map_coord, beta_coord, by = c("chr", "pos"), all = FALSE, allow.cartesian = FALSE)
keep <- keep[!duplicated(snp_ind)]

cat("Matched LDpred2 beta SNPs:", nrow(keep), "\n")
if (nrow(keep) == 0) stop("No LDpred2 betas matched target bigSNP by chr:pos.")

ind.keep <- rows_along(G)
if (nzchar(ids_keep)) {
  ids <- fread(ids_keep, header = FALSE)
  fam <- obj.bigSNP$fam
  sel <- match(paste(fam$family.ID, fam$sample.ID), paste(ids[[1]], ids[[2]]))
  ind.keep <- which(!is.na(sel))
}
cat("Individuals scored:", length(ind.keep), "\n")

fam <- obj.bigSNP$fam
scores <- data.table(FID = fam$family.ID[ind.keep], IID = fam$sample.ID[ind.keep])
metadata_rows <- list()
snps_rows <- list()

for (p in names(pathways)) {
  p_genes <- pathways[[p]]
  p_geneWin <- geneWin[names(geneWin) %in% p_genes]

  if (length(p_geneWin) == 0) {
    pathway_snps <- character()
  } else {
    ov <- findOverlaps(snpGR, p_geneWin, ignore.strand = TRUE)
    pathway_snps <- unique(snpGR$SNP[queryHits(ov)])
  }

  p_keep <- keep[rsid %in% pathway_snps]

  score_col <- make.names(p, unique = TRUE)

  if (nrow(p_keep) == 0) {
    scores[, (score_col) := NA_real_]
  } else {
    score <- big_prodVec(
      G,
      p_keep$weight,
      ind.row = ind.keep,
      ind.col = p_keep$snp_ind,
      ncores = threads
    )
    scores[, (score_col) := score]
  }

  metadata_rows[[length(metadata_rows) + 1]] <- data.table(
    pathway = p,
    score_column = score_col,
    n_genes_in_gmt = length(p_genes),
    n_genes_found_in_gtf = length(unique(names(p_geneWin))),
    n_target_snps_mapped_to_pathway = length(pathway_snps),
    n_ldpred2_snps_used = nrow(p_keep)
  )

  if (nrow(p_keep) > 0) {
    tmp <- p_keep[, .(rsid, chr, pos, map_a1, map_a0, beta_a1 = a1, beta_a0 = a0, LDpred2_beta = weight)]
    tmp[, pathway := p]
    snps_rows[[length(snps_rows) + 1]] <- tmp
  }

  cat("Pathway:", p, "| genes:", length(p_genes), "| SNPs used:", nrow(p_keep), "\n")
}

metadata <- rbindlist(metadata_rows, fill = TRUE)
snps_used <- if (length(snps_rows) > 0) rbindlist(snps_rows, fill = TRUE) else data.table()

score_file <- file.path(outdir, paste0(out_prefix_name, ".tsv.gz"))
meta_file <- file.path(outdir, paste0(out_prefix_name, ".metadata.tsv"))
snps_file <- file.path(outdir, paste0(out_prefix_name, ".snps_used.tsv.gz"))
rds_file <- file.path(outdir, paste0(out_prefix_name, ".rds"))

fwrite(scores, score_file, sep = "\t")
fwrite(metadata, meta_file, sep = "\t")
fwrite(snps_used, snps_file, sep = "\t")
saveRDS(list(scores = scores, metadata = metadata, snps_used = snps_used), rds_file)

cat("Written:\n")
cat(score_file, "\n")
cat(meta_file, "\n")
cat(snps_file, "\n")
cat(rds_file, "\n")


