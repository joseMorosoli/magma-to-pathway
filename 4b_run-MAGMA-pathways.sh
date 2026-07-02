#!/bin/bash -l

#$ -S /bin/bash
#$ -l h_rt=02:00:00
#$ -l mem=8G
#$ -pe smp 1
#$ -N MAGMA_pathways
#$ -wd /myriadfs/home/ucju659/SOFTWARE/MAGMA
#$ -j y
#$ -o /myriadfs/home/ucju659/SOFTWARE/MAGMA/logs/MAGMA_pathways.log

# =============================================================================
# Script: 4_run-MAGMA-pathways.sh
# Purpose: Run MAGMA gene-set/pathway analysis from an existing .genes.raw file.
#
# Intended use:
#   Submit as a Myriad job after script 3 has successfully created:
#     <GENE_PREFIX>.genes.raw
#
# Required input:
#   - MAGMA executable
#   - trait-specific MAGMA .genes.raw file
#   - GO and/or Reactome GMT files
#
# Output:
#   - <OUT_PREFIX>.gsa.out files for each pathway collection
#   - <OUT_PREFIX>.log files written by MAGMA
#
# Notes:
#   - This script does not rerun gene-level analysis.
#   - GO and Reactome are run as separate gene-set collections.
#   - Multiple-testing correction is handled in script 5.
#   - Edit only the trait-specific block to switch traits.
#
# No-exit mode:
#   Checks print WARNING and skip unsafe commands rather than terminating the shell.
# =============================================================================

set +e
set -o pipefail

# -----------------------------------------------------------------------------
# USER SETTINGS: edit these lines for each trait/environment.
# -----------------------------------------------------------------------------

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"

MAGMA_EXE="${PROJECT_DIR}/v1.10/magma"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

GO_GMT="${PROJECT_DIR}/pathways/short_ids/c5.go.v2026.1.Hs.entrez.shortids.gmt"
REACTOME_GMT="${PROJECT_DIR}/pathways/short_ids/c2.cp.reactome.v2026.1.Hs.entrez.shortids.gmt"

# -----------------------------------------------------------------------------
# TRAIT-SPECIFIC LINES: edit these only.
# -----------------------------------------------------------------------------

# Height
#GENE_RAW="${RESULTS_DIR}/HT_EUR_2022_genes.genes.raw"
#GO_OUT="${RESULTS_DIR}/HT_EUR_2022_GO"
#REACTOME_OUT="${RESULTS_DIR}/HT_EUR_2022_Reactome"

# BMI
GENE_RAW="${RESULTS_DIR}/BMI_EUR_2018_genes.genes.raw"
GO_OUT="${RESULTS_DIR}/BMI_EUR_2018_GO"
REACTOME_OUT="${RESULTS_DIR}/BMI_EUR_2018_Reactome"

RUN_GO=true
RUN_REACTOME=true
OVERWRITE=true

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

check_file() {
  local file="$1"
  local label="$2"

  if [[ ! -e "$file" ]]; then
    echo "WARNING: ${label} not found: $file" >&2
    return 1
  fi
  return 0
}

check_nonempty_file() {
  local file="$1"
  local label="$2"

  if [[ ! -s "$file" ]]; then
    echo "WARNING: ${label} is missing or empty: $file" >&2
    return 1
  fi
  return 0
}

run_magma_set() {
  local gene_raw="$1"
  local gmt="$2"
  local out_prefix="$3"
  local collection_name="$4"
  local can_run_collection=true

  check_nonempty_file "$gene_raw" "MAGMA .genes.raw file" || can_run_collection=false
  check_nonempty_file "$gmt" "${collection_name} GMT file" || can_run_collection=false

  if [[ "$can_run_collection" != true ]]; then
    echo "WARNING: Skipping ${collection_name} pathway analysis because required input is missing." >&2
    return 0
  fi

  if [[ "$OVERWRITE" == true ]]; then
    rm -f "${out_prefix}".*
  fi

  echo
  echo "------------------------------------------------------------"
  echo "Running MAGMA pathway analysis: ${collection_name}"
  echo "Gene results: ${gene_raw}"
  echo "GMT:          ${gmt}"
  echo "Out prefix:   ${out_prefix}"
  echo "Start time:   $(date)"
  echo "------------------------------------------------------------"

  "$MAGMA_EXE" \
    --gene-results "$gene_raw" \
    --set-annot "$gmt" \
    --out "$out_prefix"

  MAGMA_STATUS=$?

  if [[ "$MAGMA_STATUS" -ne 0 ]]; then
    echo "WARNING: MAGMA ${collection_name} pathway command returned non-zero status: $MAGMA_STATUS" >&2
  fi

  if check_nonempty_file "${out_prefix}.gsa.out" "MAGMA ${collection_name} .gsa.out output"; then
    echo "Created: ${out_prefix}.gsa.out"
  else
    echo "WARNING: Expected ${collection_name} output was not created." >&2
  fi

  echo "Finished ${collection_name}: $(date)"
}

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

CAN_RUN=true
check_file "$MAGMA_EXE" "MAGMA executable" || CAN_RUN=false
check_nonempty_file "$GENE_RAW" "MAGMA .genes.raw file" || CAN_RUN=false

if [[ -e "$MAGMA_EXE" && ! -x "$MAGMA_EXE" ]]; then
  echo "WARNING: MAGMA is not executable: $MAGMA_EXE" >&2
  echo "WARNING: Run: chmod +x \"$MAGMA_EXE\"" >&2
  CAN_RUN=false
fi

if [[ "$RUN_GO" != true && "$RUN_REACTOME" != true ]]; then
  echo "WARNING: RUN_GO and RUN_REACTOME are both false. Nothing to run." >&2
  CAN_RUN=false
fi

# -----------------------------------------------------------------------------
# Run pathway analyses
# -----------------------------------------------------------------------------

echo "MAGMA pathway analysis started: $(date)"
echo "Working directory: $(pwd)"
echo "Host: $(hostname)"
echo "Job ID: ${JOB_ID:-not_available}"
echo "Gene results: $GENE_RAW"

if [[ "$CAN_RUN" == true ]]; then
  if [[ "$RUN_GO" == true ]]; then
    run_magma_set "$GENE_RAW" "$GO_GMT" "$GO_OUT" "GO"
  fi

  if [[ "$RUN_REACTOME" == true ]]; then
    run_magma_set "$GENE_RAW" "$REACTOME_GMT" "$REACTOME_OUT" "Reactome"
  fi
else
  echo "WARNING: Skipping MAGMA pathway analysis because one or more checks failed." >&2
fi

echo
echo "MAGMA pathway analysis completed or skipped: $(date)"
echo "Expected outputs:"
if [[ "$RUN_GO" == true ]]; then
  echo "  ${GO_OUT}.gsa.out"
fi
if [[ "$RUN_REACTOME" == true ]]; then
  echo "  ${REACTOME_OUT}.gsa.out"
fi
