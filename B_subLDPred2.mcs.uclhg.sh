#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=12:00:00
#$ -l mem=12G
#$ -l tmpfs=50G
#$ -pe smp 10
#$ -t 1
#$ -N LDpred2_PGS
#$ -wd /myriadfs/home/ucju659/uclhg-mcs-pgs
#$ -j y
#$ -o /myriadfs/home/ucju659/uclhg-mcs-pgs/logs/LDpred2_PGS.log
#$ -m be
# ============================================================================
# Script: B_subLDPred2.mcs.uclhg.sh
# Purpose: Submit one or more LDpred2 PGS jobs on UCL Myriad.
#
# Intended use:
#   qsub scripts/B_subLDPred2.mcs.uclhg.sh
#
# Required input:
#   - sumstats_list.csv in WORKDIR, with no header by default.
#     Format, one trait per row:
#       <sumstats_filename>,<is_binary_trait>
#     Example:
#       GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz,FALSE
#       some_case_control_gwas.ldpred2.gz,TRUE
#
#   - LDpred2-ready summary statistics in SUMSTATDIR.
#   - A target bigsnpr .rds genotype object.
#   - HapMap3+ LD reference directory.
#   - A_ldpred2_auto_inf_qc_lift2_custom.R available at LDPRED2_SCRIPT.
#
# Output:
#   One output directory per sumstats file under OUTDIR, including LDpred2-inf
#   and LDpred2-auto scores, SNP weights, logs and diagnostics.
#
# Notes:
#   - The SGE array range (#$ -t) must match the number of rows in sumstats_list.csv.
#   - The R script is custom and includes stability settings for LDpred2-auto.
#   - Do not commit generated PGS outputs, logs, .sbk files or large RDS objects.
#   - No-exit mode: failed checks print WARNING and skip Rscript rather than
#     terminating the shell. Read the log carefully.
# ============================================================================

set +e
set -o pipefail

# ---------------------------------------------------------------------------
# Threading: keep BLAS libraries from oversubscribing cores.
# ---------------------------------------------------------------------------

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

WORKDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs"
SUMSTATS_LIST="${WORKDIR}/sumstats_list.csv"

# Recommended: keep the custom R script inside the GitHub repository/scripts folder.
LDPRED2_SCRIPT="${WORKDIR}/scripts/A_ldpred2_auto_inf_qc_lift2_custom.R"

LDREF="/myriadfs/home/ucju659/misc/hapmap3plus/"
GENOFILE="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only.rds"
OUTDIR="/myriadfs/home/ucju659/uclhg-mcs-pgs/"
SUMSTATDIR="/myriadfs/home/ucju659/SUMSTATS/ldpred2_ready"
TARGET_BUILD="hg38"
NCORES=10

# R libraries from personal Myriad folder.
export R_LIBS="/myriadfs/home/ucju659/MyRLibs/R-4.2.0/"

# ---------------------------------------------------------------------------
# Load R
# ---------------------------------------------------------------------------

module -f unload compilers mpi gcc-libs
module load r/recommended

mkdir -p "${WORKDIR}/logs" "$OUTDIR"
cd "$WORKDIR" || echo "WARNING: could not cd to WORKDIR: $WORKDIR" >&2

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

CAN_RUN=true

for path in "$SUMSTATS_LIST" "$LDPRED2_SCRIPT" "$GENOFILE"; do
  if [[ ! -s "$path" ]]; then
    echo "WARNING: required file not found or empty: $path" >&2
    CAN_RUN=false
  fi
done

if [[ ! -d "$LDREF" ]]; then
  echo "WARNING: LD reference directory not found: $LDREF" >&2
  CAN_RUN=false
fi

if [[ ! -d "$SUMSTATDIR" ]]; then
  echo "WARNING: summary-statistics directory not found: $SUMSTATDIR" >&2
  CAN_RUN=false
fi

if [[ -z "${SGE_TASK_ID:-}" ]]; then
  echo "WARNING: SGE_TASK_ID is not set. Assuming task 1 for interactive testing." >&2
  TASK_ID=1
else
  TASK_ID="$SGE_TASK_ID"
fi

# ---------------------------------------------------------------------------
# Select GWAS from CSV
# ---------------------------------------------------------------------------

LINE=""
if [[ -s "$SUMSTATS_LIST" ]]; then
  LINE=$(awk -v task="$TASK_ID" 'NR == task {print}' "$SUMSTATS_LIST" | tr -d '\r')
fi

if [[ -z "$LINE" ]]; then
  echo "WARNING: no row found in $SUMSTATS_LIST for task ID=$TASK_ID" >&2
  CAN_RUN=false
fi

SUMSTAT=""
TYPE=""
if [[ -n "$LINE" ]]; then
  SUMSTAT=$(printf '%s' "$LINE" | awk -F',' '{gsub(/^ +| +$/, "", $1); print $1}')
  TYPE=$(printf '%s' "$LINE" | awk -F',' '{gsub(/^ +| +$/, "", $2); print $2}')
fi

if [[ -z "$SUMSTAT" || -z "$TYPE" ]]; then
  echo "WARNING: malformed row in $SUMSTATS_LIST for task $TASK_ID: $LINE" >&2
  CAN_RUN=false
fi

if [[ -n "$TYPE" && "$TYPE" != "TRUE" && "$TYPE" != "FALSE" ]]; then
  echo "WARNING: trait type should be TRUE or FALSE in second CSV column. Got: $TYPE" >&2
  CAN_RUN=false
fi

if [[ -n "$SUMSTAT" && ! -s "${SUMSTATDIR}/${SUMSTAT}" ]]; then
  echo "WARNING: selected summary-statistics file not found: ${SUMSTATDIR}/${SUMSTAT}" >&2
  CAN_RUN=false
fi

# ---------------------------------------------------------------------------
# Run LDpred2
# ---------------------------------------------------------------------------

echo "============================================================"
echo "LDpred2 PGS job"
echo "Job ID:       ${JOB_ID:-not_available}"
echo "Task ID:      ${TASK_ID:-not_available}"
echo "Host:         $(hostname)"
echo "Workdir:      $WORKDIR"
echo "Sumstats:     ${SUMSTAT:-not_available}"
echo "Binary trait: ${TYPE:-not_available}"
echo "Target build: $TARGET_BUILD"
echo "Start time:   $(date)"
echo "============================================================"

if [[ "$CAN_RUN" == true ]]; then
  Rscript --vanilla "$LDPRED2_SCRIPT" \
    -s "$SUMSTAT" \
    -t "$TYPE" \
    -m "$LDREF" \
    -l "$TARGET_BUILD" \
    -g "$GENOFILE" \
    -o "$OUTDIR" \
    -c "$NCORES" \
    -d "$SUMSTATDIR"

  RSCRIPT_STATUS=$?

  if [[ "$RSCRIPT_STATUS" -ne 0 ]]; then
    echo "WARNING: LDpred2 Rscript returned non-zero status: $RSCRIPT_STATUS" >&2
  else
    echo "LDpred2 job completed: $(date)"
  fi
else
  echo "WARNING: Skipping LDpred2 Rscript because one or more checks failed." >&2
fi
