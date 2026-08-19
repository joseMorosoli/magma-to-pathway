#!/bin/bash -l
# ============================================================================
# Script: 7c_score-FTO-PRSice.sh
# Purpose: Positive-control BMI PRSet C+T score restricted to FTO.
#
# SNP-to-gene mapping:
#   FTO GRCh37 GTF boundaries
#   +35 kb at the 5' end
#   +10 kb at the 3' end
#
# PRSet performs set-specific LD clumping.
#
# Output:
#   Native PRSet output under OUTDIR
#   A slim FTO-only score file:
#     BMI_FTO_PRSet_CT.p1_only.all_score
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

THREADS="${NSLOTS:-2}"

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# ============================================================
# FTO gene set
# ============================================================

FTO_GMT="${PROJECT_DIR}/pathways/FTO.symbols.gmt"

mkdir -p "$(dirname "${FTO_GMT}")"

printf "FTO_PINGAULT\tFTO\n" > "${FTO_GMT}"

echo "FTO gene set:"
cat "${FTO_GMT}"

# ============================================================
# BMI base GWAS
# ============================================================

BASE="/myriadfs/home/ucju659/SUMSTATS/software-ready/GIANT_UKBB_BMI_2018_ALL_SITES.cleaned.tsv.gz"

SOFTWARE="/myriadfs/home/ucju659/SOFTWARE/PRSice2"
PRSICE_R="${SOFTWARE}/PRSice.R"
PRSICE_BIN="${SOFTWARE}/PRSice_linux"

TARGET_PLINK="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only_GRCh37"

GTF="/myriadfs/home/ucju659/misc/ANNOTATIONS/gtf/Homo_sapiens.GRCh37.87.gtf"

WIND5="35kb"
WIND3="10kb"

OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/PRSet/BMI_EUR_2018/FTO_prsice"

OUT_PREFIX="${OUTDIR}/BMI_FTO_PRSet_CT_GRCh37_35kb_10kb"

mkdir -p "${OUTDIR}"

# ============================================================
# Check inputs
# ============================================================

for file in \
  "${PRSICE_R}" \
  "${PRSICE_BIN}" \
  "${BASE}" \
  "${FTO_GMT}" \
  "${GTF}" \
  "${TARGET_PLINK}.bed" \
  "${TARGET_PLINK}.bim" \
  "${TARGET_PLINK}.fam"
do
  if [[ ! -e "${file}" ]]; then
    echo "ERROR: required file not found: ${file}" >&2
  fi
done

if [[ ! -x "${PRSICE_BIN}" ]]; then
  echo "ERROR: PRSice binary is not executable:" >&2
  echo "${PRSICE_BIN}" >&2
fi

# ============================================================
# Run PRSet
#
# No p-value threshold parameters are supplied:
# PRSet therefore calculates the FTO score at P = 1.
# ============================================================

echo "Running FTO PRSet positive control"
echo "Base:       ${BASE}"
echo "Gene set:   ${FTO_GMT}"
echo "Target:     ${TARGET_PLINK}"
echo "GTF:        ${GTF}"
echo "Output:     ${OUT_PREFIX}"

rm -f "${OUT_PREFIX}".*

Rscript --vanilla "${PRSICE_R}" \
  --prsice "${PRSICE_BIN}" \
  --base "${BASE}" \
  --target "${TARGET_PLINK}" \
  --snp RSID \
  --chr CHR \
  --A1 A1 \
  --A2 A2 \
  --stat BETA \
  --pvalue P \
  --base-info INFO:0.8 \
  --beta \
  --gtf "${GTF}" \
  --msigdb "${FTO_GMT}" \
  --wind-5 "${WIND5}" \
  --wind-3 "${WIND3}" \
  --clump-kb 250 \
  --clump-r2 0.1 \
  --clump-p 1 \
  --no-regress \
  --all-score \
  --print-snp \
  --nonfounders \
  --thread "${THREADS}" \
  --out "${OUT_PREFIX}"

# ============================================================
# Extract the single FTO P=1 score column
# ============================================================

ALL_SCORE="${OUT_PREFIX}.all_score"
SLIM_SCORE="${OUTDIR}/BMI_FTO_PRSet_CT.p1_only.all_score"

if [[ ! -s "${ALL_SCORE}" ]]; then
  echo "ERROR: PRSet .all_score file was not created:" >&2
  echo "${ALL_SCORE}" >&2
fi

awk '
NR == 1 {

  fid = 0
  iid = 0
  score = 0

  for (i = 1; i <= NF; i++) {

    if ($i == "FID") fid = i
    if ($i == "IID") iid = i

    if ($i ~ /^FTO_POSITIVE_CONTROL/) {
      score = i
    }
  }

  if (fid == 0 || iid == 0 || score == 0) {
    print "ERROR: could not identify FID, IID and FTO score columns." > "/dev/stderr"
    exit 2
  }

  print $fid, $iid, $score
  next
}

{
  print $fid, $iid, $score
}
' OFS="\t" "${ALL_SCORE}" > "${SLIM_SCORE}"

echo "Completed FTO PRSet positive control."
echo "FTO score file:"
echo "${SLIM_SCORE}"

echo "Dimensions:"
awk '
NR == 1 {print "columns:", NF}
END {print "rows:", NR}
' "${SLIM_SCORE}"

echo
echo "Header:"
head -n 1 "${SLIM_SCORE}"

echo
echo "PRSet SNP-membership output:"
ls -lh "${OUT_PREFIX}"*snp* 2>/dev/null || true