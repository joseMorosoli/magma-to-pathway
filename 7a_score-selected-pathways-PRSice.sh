#!/bin/bash -l
# ============================================================================
# Script: 7a_score-selected-pathways-PRSice.sh
# Purpose: Calculate PRSice/PRSet C+T pathway scores for each selected pathway.
#
# Run after:
#   6a_make-selected-pathway-gmt.R
#
# Required input:
#   - selected pathway GMT using gene symbols
#   - PRSice-ready base GWAS file
#   - target PLINK genotype files
#   - GRCh37 GTF file
#
# Output:
#   PRSet/PRSice pathway score files under OUTDIR.
# ============================================================================

set +e
set -o pipefail

module -f unload compilers mpi gcc-libs
module load r/recommended

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

THREADS="${NSLOTS:-4}"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# Choose one trait.
# TRAIT_PREFIX="HT_EUR_2022"
TRAIT_PREFIX="BMI_EUR_2018"

SELECTION="FDR05"

SELECTED_GMT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_GO_Reactome.symbols.gmt"

# Edit these base GWAS paths for height/BMI.
# The base file must contain: MarkerName CHR POS A1 A2 BETA P
# BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/HEIGHT.prsice.tsv.gz"
BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/BMI.prsice.tsv.gz"

SOFTWARE="/myriadfs/home/ucju659/SOFTWARE/PRSice2"
PRSICE_R="${SOFTWARE}/PRSice.R"
PRSICE_BIN="${SOFTWARE}/PRSice_linux"

TARGET_PLINK="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only"
GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

WIND5="35kb"
WIND3="10kb"

OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSice-pathway/${TRAIT_PREFIX}/selected_pathways_prsice"
OUT_PREFIX="${OUTDIR}/${TRAIT_PREFIX}_selected_${SELECTION}_PRSet_CT_GRCh37_35kb_10kb"

# Use "1" for p <= 1 only, or e.g. "0.001,0.01,0.05,0.1,0.5,1".
P_THRESHOLDS="1"

mkdir -p "${OUTDIR}"

CAN_RUN=true

for file in "${PRSICE_R}" "${PRSICE_BIN}" "${BASE}" "${SELECTED_GMT}" "${GTF}" "${TARGET_PLINK}.bed" "${TARGET_PLINK}.bim" "${TARGET_PLINK}.fam"; do
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
  echo "WARNING: skipping PRSice selected-pathway scores because an input is missing." >&2
else
  echo "Running PRSice/PRSet selected pathway scores"
  echo "Trait: ${TRAIT_PREFIX}"
  echo "Selected GMT: ${SELECTED_GMT}"
  echo "Base: ${BASE}"
  echo "Output prefix: ${OUT_PREFIX}"
  echo "P thresholds: ${P_THRESHOLDS}"

  Rscript --vanilla "${PRSICE_R}" \
    --prsice "${PRSICE_BIN}" \
    --base "${BASE}" \
    --target "${TARGET_PLINK}" \
    --snp MarkerName \
    --chr CHR \
    --bp POS \
    --A1 A1 \
    --A2 A2 \
    --stat BETA \
    --pvalue P \
    --beta \
    --gtf "${GTF}" \
    --msigdb "${SELECTED_GMT}" \
    --wind-5 "${WIND5}" \
    --wind-3 "${WIND3}" \
    --clump-kb 250 \
    --clump-r2 0.1 \
    --clump-p 1 \
    --bar-levels "${P_THRESHOLDS}" \
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
