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
#   A slim p1-only .all_score file:
#     <OUT_PREFIX>.p1_only.all_score
# ============================================================================

set +e
set -o pipefail

module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0
# module load r/recommended

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

SELECTED_GMT="${PROJECT_DIR}/pathways/selected/${TRAIT_PREFIX}_selected_${SELECTION}_GO_Reactome.symbols.gmt"

# Edit these base GWAS paths for height/BMI.
# The base file must contain: MarkerName CHR POS A1 A2 BETA P
# BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_HEIGHT_YENGO_2022_EUR.rsID.prsice.tsv.gz"
BASE="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_UKBB_BMI_2018_ALL_SITES.rsID.prsice.tsv.gz"

SOFTWARE="/myriadfs/home/ucju659/SOFTWARE/PRSice2"
PRSICE_R="${SOFTWARE}/PRSice.R"
PRSICE_BIN="${SOFTWARE}/PRSice_linux"

TARGET_PLINK="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only"
GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

WIND5="35kb"
WIND3="10kb"

OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/${TRAIT_PREFIX}/selected_pathways_prsice"
OUT_PREFIX="${OUTDIR}/${TRAIT_PREFIX}_selected_${SELECTION}_PRSet_CT_GRCh37_35kb_10kb"

# Use "1" for p <= 1 only, or e.g. "0.001,0.01,0.05,0.1,0.5,1".
P_THRESHOLDS="1"

mkdir -p "${OUTDIR}"

CAN_RUN=true

for file in \
  "${PRSICE_R}" \
  "${PRSICE_BIN}" \
  "${BASE}" \
  "${SELECTED_GMT}" \
  "${GTF}" \
  "${TARGET_PLINK}.bed" \
  "${TARGET_PLINK}.bim" \
  "${TARGET_PLINK}.fam"
do
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
  echo "Trait:        ${TRAIT_PREFIX}"
  echo "Selection:    ${SELECTION}"
  echo "Selected GMT: ${SELECTED_GMT}"
  echo "Base:         ${BASE}"
  echo "Output dir:   ${OUTDIR}"
  echo "Output prefix:${OUT_PREFIX}"
  echo "P thresholds: ${P_THRESHOLDS}"
  echo "Threads:      ${THREADS}"

  echo "Removing previous output files for this exact OUT_PREFIX:"
  echo "${OUT_PREFIX}.*"
  rm -f "${OUT_PREFIX}".*

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

# ============================================================================
# Final check and p1-only slimming for the trait currently being run
# ============================================================================

FILE="${OUT_PREFIX}.all_score"
OUT="${OUT_PREFIX}.p1_only.all_score"

echo "Checking PRSice output for current trait"
echo "Trait:        ${TRAIT_PREFIX}"
echo "Selection:    ${SELECTION}"
echo "Output dir:   ${OUTDIR}"
echo "All-score:    ${FILE}"
echo "Slim output:  ${OUT}"

echo "Files in output directory:"
ls -lh "${OUTDIR}"

if [[ ! -s "${FILE}" ]]; then
  echo
  echo "WARNING: .all_score file was not created or is empty:"
  echo "${FILE}"
  echo "Skipping p1-only slimming."
  exit 0
fi

echo "Original .all_score size:"
ls -lh "${FILE}"

echo "Original .all_score dimensions:"
awk 'NR == 1 {print "columns:", NF} END {print "rows:", NR}' "${FILE}"

echo "First 80 header fields:"
awk '
NR == 1 {
  for (i = 1; i <= NF; i++) print i, $i
}
' "${FILE}" | head -n 80

N_P1_COLS=$(
  awk '
  NR == 1 {
    n = 0
    for (i = 1; i <= NF; i++) {
      if ($i ~ /_1$/) n++
    }
    print n
  }
  ' "${FILE}"
)

N_PATHWAYS=$(
  wc -l < "${SELECTED_GMT}"
)

echo "Selected pathways in GMT: ${N_PATHWAYS}"
echo "Detected *_1 score columns in .all_score: ${N_P1_COLS}"

if [[ "${N_P1_COLS}" -eq 0 ]]; then
  echo
  echo "ERROR: no *_1 score columns detected in:"
  echo "${FILE}"
  echo "Cannot create p1-only slim file."
  echo "Inspect the header to confirm how PRSice named the p = 1 columns."
  exit 0
fi

echo "Creating p1-only slim .all_score file:"
echo "${OUT}"

awk '
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "FID" || $i == "IID" || $i ~ /_1$/) {
      keep[i] = 1
    }
  }
}
{
  out = ""
  for (i = 1; i <= NF; i++) {
    if (keep[i]) {
      out = out (out == "" ? "" : OFS) $i
    }
  }
  print out
}
' OFS="\t" "${FILE}" > "${OUT}"

echo "Slim .all_score size:"
ls -lh "${OUT}"

echo "Slim .all_score dimensions:"
awk 'NR == 1 {print "columns:", NF} END {print "rows:", NR}' "${OUT}"

echo "First 80 slim header fields:"
awk '
NR == 1 {
  for (i = 1; i <= NF; i++) print i, $i
}
' "${OUT}" | head -n 80

EXPECTED_SLIM_COLS=$((N_PATHWAYS + 2))

echo "Expected slim columns if all pathways are present:"
echo "${EXPECTED_SLIM_COLS} = 2 ID columns + ${N_PATHWAYS} pathway score columns"

OBSERVED_SLIM_COLS=$(
  awk 'NR == 1 {print NF}' "${OUT}"
)

if [[ "${OBSERVED_SLIM_COLS}" -ne "${EXPECTED_SLIM_COLS}" ]]; then
  echo
  echo "WARNING: observed slim column count does not equal expected count."
  echo "Observed: ${OBSERVED_SLIM_COLS}"
  echo "Expected: ${EXPECTED_SLIM_COLS}"
  echo "This may be okay only if some pathways failed or PRSice used unexpected column names."
else
  echo
  echo "Slim file has expected number of columns."
fi

echo "Completed selected-pathway PRSice scoring and p1-only slimming."
echo "Use this file for trio merging:"
echo "${OUT}"

