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
# ============================================================================

set -euo pipefail

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

for file in "$MAGMA_EXE" "$REFERENCE_BIM" "$GENE_LOC"; do
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

mkdir -p "$(dirname "$OUT_PREFIX")"

echo "Creating MAGMA gene annotation"
echo "MAGMA:        $MAGMA_EXE"
echo "Reference BIM:$REFERENCE_BIM"
echo "Gene loc:     $GENE_LOC"
echo "Output prefix:$OUT_PREFIX"
echo "Window:       ${UPSTREAM_KB}kb upstream, ${DOWNSTREAM_KB}kb downstream"
echo "Start time:   $(date)"

# ---------------------------------------------------------------------------
# Run MAGMA annotation
# ---------------------------------------------------------------------------

"$MAGMA_EXE" \
  --annotate window=${UPSTREAM_KB},${DOWNSTREAM_KB} \
  --snp-loc "$REFERENCE_BIM" \
  --gene-loc "$GENE_LOC" \
  --out "$OUT_PREFIX"

if [[ ! -s "${OUT_PREFIX}.genes.annot" ]]; then
  echo "ERROR: MAGMA did not create ${OUT_PREFIX}.genes.annot" >&2
  exit 1
fi

echo "Created: ${OUT_PREFIX}.genes.annot"
echo "End time: $(date)"
