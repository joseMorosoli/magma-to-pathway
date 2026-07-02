#!/bin/bash -l
# ============================================================================
# Script: 6c_score-nonpathway-PRSice.sh
# Purpose: Calculate a PRSice C+T score using SNPs outside the selected pathway
#          panel.
#
# Run after:
#   6b_create-nonpathway-snp-list.sh
#
# Required input:
#   - PRSice-ready base GWAS file
#   - target PLINK files
#   - nonpathway_snps.extract.txt
#
# Output:
#   PRSice score files under OUTDIR.
# ============================================================================

set +e
set -o pipefail

module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0
#module load r/recommended

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

THREADS="${NSLOTS:-2}"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# Choose one trait.
TRAIT_PREFIX="HT_EUR_2022"
#TRAIT_PREFIX="BMI_EUR_2018"

# Common values: "FDR05", "Bonferroni05", "nominal_P05", "top50"
SELECTION="Bonferroni05"

# Edit these base GWAS paths for height/BMI.
# The base file must contain: MarkerName CHR POS A1 A2 BETA P
# BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_HEIGHT_YENGO_2022_EUR.rsID.prsice.tsv.gz"
BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_UKBB_BMI_2018_ALL_SITES.rsID.prsice.tsv.gz"

SOFTWARE="/myriadfs/home/ucju659/SOFTWARE/PRSice2"
PRSICE_R="${SOFTWARE}/PRSice.R"
PRSICE_BIN="${SOFTWARE}/PRSice_linux"

TARGET_PLINK="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only"

NONPATHWAY_EXTRACT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_nonpathway/nonpathway_snps.extract.txt"

OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/${TRAIT_PREFIX}/nonpathway_prsice"
OUT_PREFIX="${OUTDIR}/${TRAIT_PREFIX}_selected_${SELECTION}_nonpathway_PRSice"

# Use "1" for p <= 1 only, or e.g. "0.001,0.01,0.05,0.1,0.5,1".
P_THRESHOLDS="1"

mkdir -p "${OUTDIR}"

CAN_RUN=true

for file in "${PRSICE_R}" "${PRSICE_BIN}" "${BASE}" "${NONPATHWAY_EXTRACT}" "${TARGET_PLINK}.bed" "${TARGET_PLINK}.bim" "${TARGET_PLINK}.fam"; do
  if [[ ! -e "${file}" ]]; then
    echo "WARNING: required file not found: ${file}" >&2
    CAN_RUN=false
  fi
done

if [[ -e "${PRSICE_BIN}" && ! -x "${PRSICE_BIN}" ]]; then
  echo "WARNING: PRSice binary is not executable: ${PRSICE_BIN}" >&2
  echo "WARNING: try chmod +x ${PRSICE_BIN}" >&2
  CAN_RUN=false
fi

if [[ "${CAN_RUN}" != true ]]; then
  echo "WARNING: skipping PRSice non-pathway score because an input is missing." >&2
else
  echo "Running PRSice non-pathway score"
  echo "Trait: ${TRAIT_PREFIX}"
  echo "Base: ${BASE}"
  echo "Extract: ${NONPATHWAY_EXTRACT}"
  echo "Output prefix: ${OUT_PREFIX}"
  echo "P thresholds: ${P_THRESHOLDS}"

  Rscript --vanilla "${PRSICE_R}" \
    --prsice "${PRSICE_BIN}" \
    --base "${BASE}" \
    --target "${TARGET_PLINK}" \
    --extract "${NONPATHWAY_EXTRACT}" \
    --snp MarkerName \
    --chr CHR \
    --bp POS \
    --A1 A1 \
    --A2 A2 \
    --stat BETA \
    --pvalue P \
    --beta \
    --clump-kb 250 \
    --clump-r2 0.1 \
    --clump-p 1 \
    --bar-levels "${P_THRESHOLDS}" \
    --fastscore \
    --no-regress \
    --all-score \
    --print-snp \
    --nonfounders \
    --thread "${THREADS}" \
    --out "${OUT_PREFIX}"

  STATUS=$?
  if [[ "${STATUS}" -ne 0 ]]; then
    echo "WARNING: PRSice returned non-zero status: ${STATUS}" >&2
  fi
fi
