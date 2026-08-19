#!/bin/bash -l
# ============================================================================
# Script: 7b_score-FTO-LDpred2.sh
# Purpose: Calculate an LDpred2 posterior-beta score for the FTO gene as a
#          positive control for BMI.
#
# The FTO gene is processed through the same pipeline as selected pathways:
#   GMT -> GRCh37 GTF -> 35 kb / 10 kb gene window -> SNPs -> LDpred2 score
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

unset R_LIBS
export R_LIBS_USER="/myriadfs/home/ucju659/MyRLibs/R-4.5.1"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export THREADS="${NSLOTS:-2}"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# ============================================================================
# Inputs
# ============================================================================

TRAIT_PREFIX="BMI_EUR_2018"

export FTO_GMT="${PROJECT_DIR}/pathways/positive_controls/FTO.symbols.gmt"

export BIGSNP_RDS="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only_GRCh37.rds"

export GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/GRCh37_GIANT_UKBB_BMI_2018_ALL_SITES.cleaned.tsv.gz_15August2026/GIANT_UKBB_BMI_2018_ALL_SITES.cleaned.tsv.gz_final_beta_auto.txt"

export OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/BMI_EUR_2018/FTO_ldpred2"

export OUT_PREFIX_NAME="${TRAIT_PREFIX}_FTO_LDpred2_GRCh37_35kb_10kb"

export UPSTREAM_BP="35000"
export DOWNSTREAM_BP="10000"

# Optional FID/IID keep file. Leave empty for all individuals.
export IDS_KEEP=""

mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${FTO_GMT}")"

# FTO GMT uses the same three-column structure expected by read_gmt():
# set name, description, gene symbol.
printf "FTO_PINGAULT\tNA\tFTO\n" > "${FTO_GMT}"

echo "FTO gene set:"
cat "${FTO_GMT}"

# Check required inputs.
for file in \
  "${FTO_GMT}" \
  "${BIGSNP_RDS}" \
  "${GTF}" \
  "${RAW_BETAS}"
do
  if [[ ! -f "${file}" ]]; then
    echo "ERROR: required file not found: ${file}" >&2
  fi
done

# ============================================================================
# Calculate FTO LDpred2 score
# ============================================================================

R --vanilla

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(bigsnpr)
  library(bigstatsr)
})

# ============================================================================
# Inputs
# ============================================================================

selected_gmt <- Sys.getenv("FTO_GMT")
bigsnp <- Sys.getenv("BIGSNP_RDS")
gtf_path <- Sys.getenv("GTF")
raw_betas <- Sys.getenv("RAW_BETAS")
outdir <- Sys.getenv("OUTDIR")
out_prefix_name <- Sys.getenv("OUT_PREFIX_NAME")
ids_keep <- Sys.getenv("IDS_KEEP")

upstream <- as.integer(Sys.getenv("UPSTREAM_BP"))
downstream <- as.integer(Sys.getenv("DOWNSTREAM_BP"))
threads <- as.integer(Sys.getenv("THREADS", unset = "1"))

normalise_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x[toupper(x) == "X"] <- "23"
  x[toupper(x) == "Y"] <- "24"
  x[toupper(x) %in% c("MT", "M")] <- "25"
  x
}

read_gmt <- function(path) {

  lines <- readLines(path, warn = FALSE)

  out <- lapply(lines, function(x) {
    x <- strsplit(x, "\t", fixed = TRUE)[[1]]
    unique(toupper(x[-c(1, 2)]))
  })

  names(out) <- vapply(
    lines,
    function(x) strsplit(x, "\t", fixed = TRUE)[[1]][1],
    character(1)
  )

  out
}

# ============================================================================
# Read FTO gene set
# ============================================================================

pathways <- read_gmt(selected_gmt)

if (length(pathways) != 1L) {
  stop("Expected exactly one gene set in FTO GMT.")
}

cat("Gene set:", names(pathways), "\n")
cat("Genes:", paste(pathways[[1]], collapse = ", "), "\n")

# ============================================================================
# Load target bigSNP
# ============================================================================

obj.bigSNP <- snp_attach(bigsnp)
G <- obj.bigSNP$genotypes
map <- as.data.table(obj.bigSNP$map)

names(map)[1:6] <- c(
  "chr",
  "rsid",
  "dist",
  "pos",
  "a1",
  "a0"
)

map[, chr := as.integer(normalise_chr(chr))]
map[, pos := as.integer(pos)]
map[, rsid := as.character(rsid)]
map[, a0 := toupper(as.character(a0))]
map[, a1 := toupper(as.character(a1))]

map_small <- map[, .(
  chr,
  pos,
  rsid,
  a0,
  a1
)]

cat("Target SNPs:", nrow(map_small), "\n")

# Keep genotype-column index attached to each SNP coordinate.
snpGR <- GRanges(
  seqnames = paste0("chr", map_small$chr),
  ranges = IRanges(
    start = map_small$pos,
    end = map_small$pos
  )
)

snpGR$snp_ind <- seq_len(nrow(map_small))

# ============================================================================
# Load GTF and define gene windows
# ============================================================================

gtf <- rtracklayer::import(gtf_path)
genes <- gtf[gtf$type == "gene"]

if (!"gene_name" %in% names(mcols(genes))) {
  stop("gene_name not found in GTF.")
}

genes$gene_clean <- toupper(
  as.character(genes$gene_name)
)

genes <- genes[
  !is.na(genes$gene_clean) &
  nzchar(genes$gene_clean)
]

geneWin <- genes

st <- as.character(strand(genes))
new_start <- start(genes)
new_end <- end(genes)

plus <- st == "+"
minus <- st == "-"

# 35 kb at the 5' end and 10 kb at the 3' end.
new_start[plus] <- pmax(
  1L,
  start(genes)[plus] - upstream
)

new_end[plus] <- (
  end(genes)[plus] + downstream
)

new_start[minus] <- pmax(
  1L,
  start(genes)[minus] - downstream
)

new_end[minus] <- (
  end(genes)[minus] + upstream
)

ranges(geneWin) <- IRanges(
  new_start,
  new_end
)

names(geneWin) <- genes$gene_clean

geneWin <- keepStandardChromosomes(
  geneWin,
  pruning.mode = "coarse"
)

snpGR <- keepStandardChromosomes(
  snpGR,
  pruning.mode = "coarse"
)

seqlevelsStyle(snpGR) <- seqlevelsStyle(geneWin)[1]

# ============================================================================
# Load genome-wide BMI LDpred2 posterior betas
#
# Expected columns:
# chr pos a0 a1 final_beta_auto
# ============================================================================

betas <- fread(
  raw_betas,
  select = c(
    "chr",
    "pos",
    "a0",
    "a1",
    "final_beta_auto"
  )
)

setnames(
  betas,
  "final_beta_auto",
  "beta"
)

betas[, chr := as.integer(normalise_chr(chr))]
betas[, pos := as.integer(pos)]
betas[, a0 := toupper(as.character(a0))]
betas[, a1 := toupper(as.character(a1))]
betas[, beta := as.numeric(beta)]

betas <- betas[
  !is.na(chr) &
  !is.na(pos) &
  !is.na(a0) &
  !is.na(a1) &
  !is.na(beta)
]

cat("LDpred2 beta rows:", nrow(betas), "\n")

# ============================================================================
# Allele-aware matching
#
# snp_match() aligns posterior betas to the target genotype orientation.
# ============================================================================

keep <- as.data.table(
  bigsnpr::snp_match(
    betas,
    map_small[, .(
      chr,
      pos,
      a0,
      a1
    )],
    strand_flip = FALSE,
    join_by_pos = TRUE,
    remove_dups = TRUE,
    return_flip_and_rev = TRUE
  )
)

if (nrow(keep) == 0L) {
  stop("No LDpred2 betas matched the target bigSNP.")
}

# _NUM_ID_ indexes the second snp_match input, i.e. the target map.
keep[, snp_ind := `_NUM_ID_`]
keep[, rsid := map_small$rsid[snp_ind]]

cat("Matched LDpred2 SNPs:", nrow(keep), "\n")
cat(
  "Alleles reversed:",
  sum(keep$`_REV_`, na.rm = TRUE),
  "\n"
)

# ============================================================================
# Individuals to score
# ============================================================================

fam <- obj.bigSNP$fam
ind.keep <- rows_along(G)

if (nzchar(ids_keep)) {

  ids <- fread(
    ids_keep,
    header = FALSE
  )

  wanted <- paste(
    ids[[1]],
    ids[[2]]
  )

  target <- paste(
    fam$family.ID,
    fam$sample.ID
  )

  ind.keep <- which(
    target %in% wanted
  )
}

cat(
  "Individuals scored:",
  length(ind.keep),
  "\n"
)

# ============================================================================
# Calculate FTO-restricted LDpred2 score
#
# This is intentionally the same pathway loop used in the selected-pathway
# pipeline. Here it runs only once because the GMT contains only FTO.
# ============================================================================

scores <- data.table(
  FID = fam$family.ID[ind.keep],
  IID = fam$sample.ID[ind.keep]
)

metadata_rows <- list()
snps_rows <- list()

for (p in names(pathways)) {

  p_genes <- pathways[[p]]

  p_geneWin <- geneWin[
    names(geneWin) %in% p_genes
  ]

  if (length(p_geneWin) == 0L) {

    pathway_ind <- integer()

  } else {

    ov <- findOverlaps(
      snpGR,
      p_geneWin,
      ignore.strand = TRUE
    )

    pathway_ind <- unique(
      snpGR$snp_ind[queryHits(ov)]
    )
  }

  # Restrict matched genome-wide LDpred2 weights to FTO SNPs.
  p_keep <- keep[
    snp_ind %in% pathway_ind
  ]

  score_col <- make.names(p)

  if (nrow(p_keep) == 0L) {

    stop(
      "No matched LDpred2 SNPs found for ",
      p
    )

  } else {

    scores[, (score_col) := big_prodVec(
      G,
      p_keep$beta,
      ind.row = ind.keep,
      ind.col = p_keep$snp_ind,
      ncores = threads
    )]
  }

  metadata_rows[[p]] <- data.table(
    pathway = p,
    score_column = score_col,
    n_genes_in_gmt = length(p_genes),
    n_genes_found_in_gtf = length(
      unique(names(p_geneWin))
    ),
    n_target_snps_mapped = length(pathway_ind),
    n_ldpred2_snps_used = nrow(p_keep)
  )

  snps_rows[[p]] <- p_keep[, .(
    pathway = p,
    rsid,
    chr,
    pos,
    a0,
    a1,
    allele_reversed = `_REV_`,
    LDpred2_beta = beta
  )]

  cat(
    p,
    "| genes:",
    length(p_genes),
    "| genes found:",
    length(unique(names(p_geneWin))),
    "| target SNPs:",
    length(pathway_ind),
    "| LDpred2 SNPs:",
    nrow(p_keep),
    "\n"
  )
}

metadata <- rbindlist(metadata_rows)
snps_used <- rbindlist(snps_rows)

# ============================================================================
# Save outputs
# ============================================================================

score_file <- file.path(
  outdir,
  paste0(
    out_prefix_name,
    ".tsv.gz"
  )
)

meta_file <- file.path(
  outdir,
  paste0(
    out_prefix_name,
    ".metadata.tsv"
  )
)

snps_file <- file.path(
  outdir,
  paste0(
    out_prefix_name,
    ".snps_used.tsv.gz"
  )
)

rds_file <- file.path(
  outdir,
  paste0(
    out_prefix_name,
    ".rds"
  )
)

fwrite(
  scores,
  score_file,
  sep = "\t"
)

fwrite(
  metadata,
  meta_file,
  sep = "\t"
)

fwrite(
  snps_used,
  snps_file,
  sep = "\t"
)

saveRDS(
  list(
    scores = scores,
    metadata = metadata,
    snps_used = snps_used
  ),
  rds_file
)

cat("\nCompleted FTO LDpred2 positive control.\n")
cat("Scores:   ", score_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("SNPs:     ", snps_file, "\n", sep = "")
cat("RDS:      ", rds_file, "\n", sep = "")

