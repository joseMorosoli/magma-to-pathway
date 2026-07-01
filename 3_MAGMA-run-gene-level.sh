#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=12:00:00
#$ -l mem=8G
#$ -pe smp 1
#$ -N MAGMA_gene
#$ -wd /myriadfs/home/ucju659/SOFTWARE/MAGMA
#$ -j y
#$ -o /myriadfs/home/ucju659/SOFTWARE/MAGMA/logs/MAGMA_gene.log
# ============================================================================
# Script: 3_MAGMA-run-gene-level.sh
# Purpose: Run MAGMA gene-level analysis for one trait at a time.
#
# Intended use:
#   Submit to Myriad with:
#     qsub scripts/3_MAGMA-run-gene-level.sh
#
# Required input:
#   - MAGMA executable
#   - Reference PLINK files: .bed/.bim/.fam
#   - MAGMA-ready summary statistics with SNP, P, N columns
#   - MAGMA .genes.annot file from script 2
#
# Output:
#   <GENE_OUT>.genes.raw
#   <GENE_OUT>.genes.out
#   <GENE_OUT>.log
#
# To switch trait, edit only MUNGED_FILE and GENE_OUT in USER SETTINGS.
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# USER SETTINGS: edit these two lines for each trait.
# ---------------------------------------------------------------------------

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

# Example: height
MUNGED_FILE="${PROJECT_DIR}/munged/HT_EUR_2022_MAGMA.txt"
GENE_OUT="${PROJECT_DIR}/results/HT_EUR_2022_genes"

# Example: BMI. Uncomment these two lines and comment height lines above.
# MUNGED_FILE="${PROJECT_DIR}/munged/BMI_all.sites_2022_MAGMA.txt"
# GENE_OUT="${PROJECT_DIR}/results/BMI_EUR_2018_genes"

# ---------------------------------------------------------------------------
# Shared project paths
# ---------------------------------------------------------------------------

MAGMA_EXE="${PROJECT_DIR}/v1.10/magma"
REFERENCE_PREFIX="${PROJECT_DIR}/g1000_eur/g1000_eur"
GENE_ANNOT="${PROJECT_DIR}/gene_locations/GRCh37_35kb_10kb.genes.annot"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

# If TRUE, existing files with this output prefix are deleted before running.
OVERWRITE="TRUE"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

echo "============================================================"
echo "MAGMA gene-level analysis"
echo "Job ID:      ${JOB_ID:-not_available}"
echo "Host:        $(hostname)"
echo "Start time:  $(date)"
echo "Munged file: $MUNGED_FILE"
echo "Output:      $GENE_OUT"
echo "============================================================"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

for file in \
  "$MAGMA_EXE" \
  "$MUNGED_FILE" \
  "${REFERENCE_PREFIX}.bed" \
  "${REFERENCE_PREFIX}.bim" \
  "${REFERENCE_PREFIX}.fam" \
  "$GENE_ANNOT"
do
  if [[ ! -e "$file" ]]; then
    echo "ERROR: required file not found: $file" >&2
    exit 1
  fi
done

if [[ ! -x "$MAGMA_EXE" ]]; then
  echo "ERROR: MAGMA is not executable: $MAGMA_EXE" >&2
  echo "Fix with: chmod +x \"$MAGMA_EXE\"" >&2
  exit 1
fi

HEADER_FIELDS=$(
  head -n 1 "$MUNGED_FILE" |
    tr -d '\r' |
    tr '\t ,' '\n' |
    sed '/^$/d'
)

for col in SNP P N; do
  if ! grep -Fxq "$col" <<< "$HEADER_FIELDS"; then
    echo "ERROR: required column '$col' not found in $MUNGED_FILE" >&2
    echo "Observed header:" >&2
    head -n 1 "$MUNGED_FILE" | cat -A >&2
    exit 1
  fi
done

if [[ "$OVERWRITE" == "TRUE" ]]; then
  rm -f "${GENE_OUT}".*
elif compgen -G "${GENE_OUT}.*" > /dev/null; then
  echo "ERROR: output files already exist for prefix: $GENE_OUT" >&2
  echo "Set OVERWRITE=TRUE or choose a new GENE_OUT prefix." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Run MAGMA gene-level analysis
# ---------------------------------------------------------------------------

"$MAGMA_EXE" \
  --bfile "$REFERENCE_PREFIX" \
  --pval "$MUNGED_FILE" ncol=N \
  --gene-annot "$GENE_ANNOT" \
  --out "$GENE_OUT"

# ---------------------------------------------------------------------------
# Final checks
# ---------------------------------------------------------------------------

if [[ ! -s "${GENE_OUT}.genes.raw" ]]; then
  echo "ERROR: MAGMA did not create ${GENE_OUT}.genes.raw" >&2
  exit 1
fi

if [[ ! -s "${GENE_OUT}.genes.out" ]]; then
  echo "ERROR: MAGMA did not create ${GENE_OUT}.genes.out" >&2
  exit 1
fi

echo "MAGMA gene-level analysis completed successfully."
echo "Created:"
echo "  ${GENE_OUT}.genes.raw"
echo "  ${GENE_OUT}.genes.out"
echo "  ${GENE_OUT}.log"
echo "End time: $(date)"
