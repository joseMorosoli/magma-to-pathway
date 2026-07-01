#!/bin/bash -l
# ============================================================================
# Script: 2_create-annotation-file.sh
# Purpose: Create the reusable MAGMA SNP-to-gene annotation file.
#
# Intended use:
#   Run once on Myriad or locally, provided MAGMA and reference files exist.
#   This step is common to all traits and should only be rerun if any of these
#   change: reference .bim, gene-location file, genome build, or annotation window.
#
# Required input:
#   - MAGMA executable
#   - Reference .bim file
#   - NCBI gene-location file for the matching genome build
#
# Output:
#   <OUT_PREFIX>.genes.annot
#
# Current project convention:
#   Genome build: GRCh37
#   Window: 35 kb upstream, 10 kb downstream
#
# No-exit mode:
#   Checks print WARNING and skip the MAGMA command rather than terminating the shell.
#   This prevents accidental logout if the script is sourced interactively.
# ============================================================================

set +e
set -o pipefail

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

PROJECT_DIR="/myriadfs/home/ucju659/SOFTWARE/MAGMA"
MAGMA_EXE="${PROJECT_DIR}/v1.10/magma"
REFERENCE_BIM="${PROJECT_DIR}/g1000_eur/g1000_eur.bim"
GENE_LOC="${PROJECT_DIR}/gene_locations/NCBI37.3.gene.loc"
OUT_PREFIX="${PROJECT_DIR}/gene_locations/GRCh37_35kb_10kb"

UPSTREAM_KB=35
DOWNSTREAM_KB=10

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

CAN_RUN=true

for file in "$MAGMA_EXE" "$REFERENCE_BIM" "$GENE_LOC"; do
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

mkdir -p "$(dirname "$OUT_PREFIX")"

if [[ "$CAN_RUN" == true ]]; then
  echo "Creating MAGMA gene annotation"
  echo "MAGMA:        $MAGMA_EXE"
  echo "Reference BIM:$REFERENCE_BIM"
  echo "Gene loc:     $GENE_LOC"
  echo "Output prefix:$OUT_PREFIX"
  echo "Window:       ${UPSTREAM_KB}kb upstream, ${DOWNSTREAM_KB}kb downstream"
  echo "Start time:   $(date)"

  "$MAGMA_EXE" \
    --annotate window=${UPSTREAM_KB},${DOWNSTREAM_KB} \
    --snp-loc "$REFERENCE_BIM" \
    --gene-loc "$GENE_LOC" \
    --out "$OUT_PREFIX"

  MAGMA_STATUS=$?

  if [[ "$MAGMA_STATUS" -ne 0 ]]; then
    echo "WARNING: MAGMA annotation command returned non-zero status: $MAGMA_STATUS" >&2
  fi

  if [[ ! -s "${OUT_PREFIX}.genes.annot" ]]; then
    echo "WARNING: MAGMA did not create a non-empty annotation file: ${OUT_PREFIX}.genes.annot" >&2
  else
    echo "Created: ${OUT_PREFIX}.genes.annot"
  fi

  echo "End time: $(date)"
else
  echo "WARNING: Skipping MAGMA annotation because one or more required checks failed." >&2
fi
