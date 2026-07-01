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
#   It can also be run interactively after editing the USER SETTINGS block.
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
#
# No-exit mode:
#   Checks print WARNING and skip the MAGMA command rather than terminating the shell.
#   This prevents accidental logout if the script is sourced interactively.
#   For submitted jobs, read the log carefully because warnings may not mark the
#   job as failed to the scheduler.
# ============================================================================

set +e
set -o pipefail

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

CAN_RUN=true

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

echo "============================================================"
echo "MAGMA gene-level analysis"
echo "Job ID:      ${JOB_ID:-not_available}"
echo "Host:        $(hostname)"
echo "Start time:  $(date)"
echo "Munged file: $MUNGED_FILE"
echo "Output:      $GENE_OUT"
echo "============================================================"

for file in \
  "$MAGMA_EXE" \
  "$MUNGED_FILE" \
  "${REFERENCE_PREFIX}.bed" \
  "${REFERENCE_PREFIX}.bim" \
  "${REFERENCE_PREFIX}.fam" \
  "$GENE_ANNOT"
do
  if [[ ! -e "$file" ]]; then
    echo "WARNING: required file not found: $file" >&2
    CAN_RUN=false
  fi
done

if [[ -e "$MAGMA_EXE" && ! -x "$MAGMA_EXE" ]]; then
  echo "WARNING: MAGMA is not executable: $MAGMA_EXE" >&2
  echo "WARNING: Fix with: chmod +x \"$MAGMA_EXE\"" >&2
  CAN_RUN=false
fi

if [[ -s "$MUNGED_FILE" ]]; then
  HEADER_FIELDS=$(head -n 1 "$MUNGED_FILE" | tr -d '\r' | tr '\t ,' '\n' | sed '/^$/d')

  for col in SNP P N; do
    if ! grep -Fxq "$col" <<< "$HEADER_FIELDS"; then
      echo "WARNING: required column '$col' not found in $MUNGED_FILE" >&2
      echo "WARNING: Observed header: $(head -n 1 "$MUNGED_FILE" | cat -A)" >&2
      CAN_RUN=false
    fi
  done
else
  echo "WARNING: cannot check header because MUNGED_FILE is missing or empty: $MUNGED_FILE" >&2
  CAN_RUN=false
fi

if [[ "$OVERWRITE" == "TRUE" ]]; then
  if [[ "$CAN_RUN" == true ]]; then
    rm -f "${GENE_OUT}".*
  fi
elif compgen -G "${GENE_OUT}.*" > /dev/null; then
  echo "WARNING: output files already exist for prefix: $GENE_OUT" >&2
  echo "WARNING: Set OVERWRITE=TRUE or choose a new GENE_OUT prefix. Skipping MAGMA." >&2
  CAN_RUN=false
fi

# ---------------------------------------------------------------------------
# Run MAGMA gene-level analysis
# ---------------------------------------------------------------------------

if [[ "$CAN_RUN" == true ]]; then
  "$MAGMA_EXE" \
    --bfile "$REFERENCE_PREFIX" \
    --pval "$MUNGED_FILE" ncol=N \
    --gene-annot "$GENE_ANNOT" \
    --out "$GENE_OUT"

  MAGMA_STATUS=$?

  if [[ "$MAGMA_STATUS" -ne 0 ]]; then
    echo "WARNING: MAGMA gene-level command returned non-zero status: $MAGMA_STATUS" >&2
  fi

  # -------------------------------------------------------------------------
  # Final checks
  # -------------------------------------------------------------------------

  if [[ ! -s "${GENE_OUT}.genes.raw" ]]; then
    echo "WARNING: MAGMA did not create a non-empty file: ${GENE_OUT}.genes.raw" >&2
  fi

  if [[ ! -s "${GENE_OUT}.genes.out" ]]; then
    echo "WARNING: MAGMA did not create a non-empty file: ${GENE_OUT}.genes.out" >&2
  fi

  if [[ -s "${GENE_OUT}.genes.raw" && -s "${GENE_OUT}.genes.out" ]]; then
    echo "MAGMA gene-level analysis completed successfully."
    echo "Created:"
    echo "  ${GENE_OUT}.genes.raw"
    echo "  ${GENE_OUT}.genes.out"
    echo "  ${GENE_OUT}.log"
  else
    echo "WARNING: MAGMA gene-level analysis did not produce all expected outputs." >&2
  fi
else
  echo "WARNING: Skipping MAGMA gene-level analysis because one or more checks failed." >&2
fi

echo "End time: $(date)"
