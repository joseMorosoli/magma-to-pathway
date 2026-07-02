#!/bin/bash -l
# ============================================================================
# Script: 6d_score-nonpathway-LDpred2.sh
# Purpose: Calculate a non-pathway LDpred2 score using posterior LDpred2 betas
#          restricted to SNPs outside the selected pathway panel.
#
# Run after:
#   6b_create-nonpathway-snp-list.sh
#
# Required input:
#   - target bigsnpr .rds
#   - LDpred2 posterior beta file with columns chr pos a0 a1 beta
#   - nonpathway_snps.extract.txt
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
TRAIT_PREFIX="HT_EUR_2022"
#TRAIT_PREFIX="BMI_EUR_2018"

# Common values: "FDR05", "Bonferroni05", "nominal_P05", "top50"
SELECTION="Bonferroni05"

export BIGSNP_RDS="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only.rds"

# Edit this for height/BMI.
# Expected columns: chr pos a0 a1 beta, with or without a header.
export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz_17June2026/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz_final_beta_auto.txt"
#export RAW_BETAS="/myriadfs/home/ucju659/uclhg-mcs-pgs/GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz_17June2026/GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz_final_beta_auto.txt"

export NONPATHWAY_EXTRACT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_nonpathway/nonpathway_snps.extract.txt"

export OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/${TRAIT_PREFIX}/nonpathway_ldpred2"
export OUT_PREFIX_NAME="${TRAIT_PREFIX}_selected_${SELECTION}_nonpathway_LDpred2"

# Optional FID/IID keep file. Leave empty for all individuals.
export IDS_KEEP=""

mkdir -p "${OUTDIR}"

CAN_RUN=true
for file in "${BIGSNP_RDS}" "${RAW_BETAS}" "${NONPATHWAY_EXTRACT}"; do
  if [[ ! -f "${file}" ]]; then
    echo "WARNING: required file not found: ${file}" >&2
    CAN_RUN=false
  fi
done

if [[ "${CAN_RUN}" != true ]]; then
  echo "WARNING: skipping LDpred2 non-pathway score because an input is missing." >&2
fi

R --vanilla

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
})

bigsnp <- Sys.getenv("BIGSNP_RDS")
raw_betas <- Sys.getenv("RAW_BETAS")
extract <- Sys.getenv("NONPATHWAY_EXTRACT")
outdir <- Sys.getenv("OUTDIR")
out_prefix_name <- Sys.getenv("OUT_PREFIX_NAME")
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

cat("Loading extract list:\n", extract, "\n", sep = "")
extract_ids <- fread(extract, header = FALSE)[[1]]
extract_ids <- unique(as.character(extract_ids))
cat("SNPs in extract:", length(extract_ids), "\n")

cat("Loading bigSNP:\n", bigsnp, "\n", sep = "")
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
cat("Extract-bigSNP rsID overlap:", length(intersect(extract_ids, map_small$rsid)), "\n")

cat("Loading LDpred2 betas:\n", raw_betas, "\n", sep = "")
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

cat("LDpred2 beta rows:", nrow(betas), "\n")

# LDpred2 posterior beta files usually match target genotype by chr:pos.
map_coord <- map_small[!duplicated(paste(chr, pos)) & !duplicated(paste(chr, pos), fromLast = TRUE)]
beta_coord <- betas[!duplicated(paste(chr, pos)) & !duplicated(paste(chr, pos), fromLast = TRUE)]

keep <- merge(map_coord, beta_coord, by = c("chr", "pos"), all = FALSE, allow.cartesian = FALSE)
keep <- keep[!duplicated(snp_ind)]

cat("Matched LDpred2 beta SNPs:", nrow(keep), "\n")

nonpathway <- keep[rsid %in% extract_ids]
cat("Matched non-pathway LDpred2 SNPs:", nrow(nonpathway), "\n")

if (nrow(nonpathway) == 0) stop("No matched non-pathway LDpred2 SNPs.")

if (all(c("a0", "a1", "map_a0", "map_a1") %in% names(nonpathway))) {
  same <- nonpathway$a0 == nonpathway$map_a0 & nonpathway$a1 == nonpathway$map_a1
  swap <- nonpathway$a0 == nonpathway$map_a1 & nonpathway$a1 == nonpathway$map_a0
  mismatch <- !(same | swap)
  cat("Allele same orientation:", sum(same, na.rm = TRUE), "\n")
  cat("Allele swapped orientation:", sum(swap, na.rm = TRUE), "\n")
  cat("Allele mismatches:", sum(mismatch, na.rm = TRUE), "\n")
}

ind.keep <- rows_along(G)

if (nzchar(ids_keep)) {
  ids <- fread(ids_keep, header = FALSE)
  fam <- obj.bigSNP$fam
  sel <- match(paste(fam$family.ID, fam$sample.ID), paste(ids[[1]], ids[[2]]))
  ind.keep <- which(!is.na(sel))
}

cat("Individuals scored:", length(ind.keep), "\n")

score <- big_prodVec(
  G,
  nonpathway$weight,
  ind.row = ind.keep,
  ind.col = nonpathway$snp_ind,
  ncores = threads
)

fam <- obj.bigSNP$fam
scores <- data.table(FID = fam$family.ID[ind.keep], IID = fam$sample.ID[ind.keep])
scores[, (out_prefix_name) := score]

score_file <- file.path(outdir, paste0(out_prefix_name, ".tsv.gz"))
meta_file <- file.path(outdir, paste0(out_prefix_name, ".metadata.tsv"))
snps_file <- file.path(outdir, paste0(out_prefix_name, ".snps_used.tsv.gz"))
rds_file <- file.path(outdir, paste0(out_prefix_name, ".rds"))

metadata <- data.table(
  bigsnp = bigsnp,
  raw_betas = raw_betas,
  nonpathway_extract = extract,
  out_prefix_name = out_prefix_name,
  n_extract_snps = length(extract_ids),
  n_matched_ldpred2_betas = nrow(keep),
  n_nonpathway_snps_used = nrow(nonpathway),
  n_individuals_scored = length(ind.keep)
)

snps_used <- nonpathway[, .(rsid, chr, pos, map_a1, map_a0, beta_a1 = a1, beta_a0 = a0, LDpred2_beta = weight)]

fwrite(scores, score_file, sep = "\t")
fwrite(metadata, meta_file, sep = "\t")
fwrite(snps_used, snps_file, sep = "\t")
saveRDS(list(scores = scores, metadata = metadata, snps_used = snps_used), rds_file)

cat("Written:\n")
cat(score_file, "\n")
cat(meta_file, "\n")
cat(snps_file, "\n")
cat(rds_file, "\n")
