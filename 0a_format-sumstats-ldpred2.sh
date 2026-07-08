# ============================================================================
# Script: 0_format-sumstats-ldpred2.sh
# Purpose: Convert raw GWAS summary statistics to a simple LDpred2-ready format.
#
# Intended use:
#   Run interactively on UCL Myriad, line by line or as a normal shell script.
#   This script is deliberately one-trait-at-a-time. To switch trait, edit only
#   the USER SETTINGS section.
#
# Required input:
#   A raw GWAS summary-statistics file, optionally gzipped.
#   Required conceptual columns: chromosome, base-pair position, effect allele,
#   other allele, beta, standard error, p-value, and either per-SNP N or fallback N.
#
# Output:
#   A gzipped tab-delimited file with columns:
#     CHR BP A2 A1 N BETA SE MAF P INFO
#
# Notes:
#   - A1 is the effect allele.
#   - A2 is the other/non-effect allele.
#   - Ambiguous A/T and C/G SNPs are removed.
#   - Duplicate physical positions are removed, keeping the first occurrence.
#   - If the raw file has no N column, set N_GWAS to a real defensible value.
#   - Do not commit raw or formatted summary statistics to GitHub.
#   - No-exit mode: checks print WARNING and skip unsafe commands rather than
#     terminating the shell. This is safer if you accidentally source the script in an
#     interactive Myriad session, but it also means you must read the warnings.
# ============================================================================

# ---------------------------------------------------------------------------
# Optional Myriad interactive session command. Run this in the terminal before
# the script if the login node should not do the processing.
# ---------------------------------------------------------------------------
# qrsh -pe smp 1 -l mem=16G,h_rt=0:30:00 -now no

# ---------------------------------------------------------------------------
# USER SETTINGS: edit these lines for each GWAS/trait.
# ---------------------------------------------------------------------------

IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/BMI_GIANT_UKB_2018_all_sites.txt.gz"
OUTDIR="/myriadfs/home/ucju659/SUMSTATS/ldpred2_ready"
OUT="${OUTDIR}/GIANT_UKBB_BMI_2018_ALL_SITES.test.ldpred2.gz"

# Used only when no N column exists in the raw summary statistics.
# Set to a real value, not a placeholder, before running such files.
N_GWAS="PUT_REAL_N_HERE"

# Optional filters.
MIN_MAF="0.01"
MIN_INFO="0.8"     # e.g. 0.6, or NA to skip INFO filtering here.

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

mkdir -p "$OUTDIR"
CAN_RUN=true

if [[ ! -s "$IN" ]]; then
  echo "WARNING: input summary statistics not found or empty: $IN" >&2
  echo "WARNING: Skipping LDpred2 formatting. Edit IN or check the file path." >&2
  CAN_RUN=false
fi

if [[ "$IN" =~ \.gz$ ]]; then
  READCMD="zcat"
else
  READCMD="cat"
fi

if [[ "$CAN_RUN" == true ]]; then
  echo "Input:          $IN"
  echo "Output:         $OUT"
  echo "N_GWAS fallback:$N_GWAS"
  echo "MIN_MAF:        $MIN_MAF"
  echo "MIN_INFO:       $MIN_INFO"

  # -------------------------------------------------------------------------
  # Format summary statistics
  # -------------------------------------------------------------------------

  TMP_OUT="${OUT}.tmp.$$"

  $READCMD "$IN" | awk \
    -v FS='[,\t ]+' \
    -v OFS='\t' \
    -v N_FALLBACK="$N_GWAS" \
    -v MIN_MAF="$MIN_MAF" \
    -v MIN_INFO="$MIN_INFO" '
  BEGIN {
    IGNORECASE = 1
  }

  function clean(x) {
    gsub(/\r/, "", x)
    gsub(/^ +| +$/, "", x)
    return x
  }

  function findcol(alias_string,    n,i,a,key) {
    n = split(alias_string, a, "|")
    for (i = 1; i <= n; i++) {
      key = toupper(a[i])
      if (key in h) return h[key]
    }
    return 0
  }

  NR == 1 {
    for (i = 1; i <= NF; i++) {
      col = toupper(clean($i))
      h[col] = i
    }

    cCHR  = findcol("CHR|CHROM|CHROMOSOME")
    cBP   = findcol("BP|POS|POSITION|BASE_PAIR_LOCATION")
    cA1   = findcol("A1|EFFECT_ALLELE|EA|ALLELE1|TESTED_ALLELE")
    cA2   = findcol("A2|OTHER_ALLELE|NON_EFFECT_ALLELE|NEA|ALLELE2")
    cBETA = findcol("BETA|EFFECT|B")
    cSE   = findcol("SE|STDERR|STANDARD_ERROR")
    cP    = findcol("P|PVAL|PVALUE|P_VALUE")

    cN    = findcol("N|TOTAL_N|NMISS|N_TOTAL")
    cMAF  = findcol("MAF|MINOR_ALLELE_FREQ|MINOR_ALLELE_FREQUENCY")
    cEAF  = findcol("EAF|EAF_HRC|EFFECT_ALLELE_FREQ|EFFECT_ALLELE_FREQUENCY|AF1|FREQ_TESTED_ALLELE")
    cINFO = findcol("INFO|INFO_SCORE|IMPINFO|RSQ")

    ok = 1
    if (!cCHR)  { print "WARNING: could not identify chromosome column" > "/dev/stderr"; ok = 0 }
    if (!cBP)   { print "WARNING: could not identify base-pair position column" > "/dev/stderr"; ok = 0 }
    if (!cA1)   { print "WARNING: could not identify effect allele column" > "/dev/stderr"; ok = 0 }
    if (!cA2)   { print "WARNING: could not identify other allele column" > "/dev/stderr"; ok = 0 }
    if (!cBETA) { print "WARNING: could not identify beta/effect column" > "/dev/stderr"; ok = 0 }
    if (!cSE)   { print "WARNING: could not identify standard error column" > "/dev/stderr"; ok = 0 }
    if (!cP)    { print "WARNING: could not identify p-value column" > "/dev/stderr"; ok = 0 }

    if (!cN && (N_FALLBACK == "" || N_FALLBACK == "NA" || N_FALLBACK == "PUT_REAL_N_HERE")) {
      print "WARNING: no N column found and N_GWAS fallback is not set" > "/dev/stderr"
      ok = 0
    }

    if (!ok) {
      print "WARNING: not creating a formatted LDpred2 file because required columns are missing" > "/dev/stderr"
      skip_all = 1
      next
    }

    print "CHR", "BP", "A2", "A1", "N", "BETA", "SE", "MAF", "P", "INFO"
    next
  }

  NR > 1 {
    if (skip_all == 1) next
    chr  = clean($cCHR)
    bp   = clean($cBP)
    a1   = toupper(clean($cA1))
    a2   = toupper(clean($cA2))
    beta = clean($cBETA)
    se   = clean($cSE)
    p    = clean($cP)

    sub(/^chr/, "", chr)
    sub(/^CHR/, "", chr)

    if (cN) n = clean($cN); else n = N_FALLBACK

    maf = "NA"
    if (cMAF && clean($cMAF) != "" && clean($cMAF) != "NA" && clean($cMAF) != "NaN") {
      maf = clean($cMAF) + 0
      if (maf > 0.5) maf = 1 - maf
    } else if (cEAF && clean($cEAF) != "" && clean($cEAF) != "NA" && clean($cEAF) != "NaN") {
      maf = clean($cEAF) + 0
      if (maf > 0.5) maf = 1 - maf
    }

    if (cINFO && clean($cINFO) != "" && clean($cINFO) != "NA" && clean($cINFO) != "NaN") info = clean($cINFO)
    else info = "NA"

    if (chr !~ /^[0-9]+$/ || chr < 1 || chr > 22) next
    if (bp !~ /^[0-9]+$/) next
    if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/) next
    if (a1 == a2) next
    if ((a1 == "A" && a2 == "T") || (a1 == "T" && a2 == "A") ||
        (a1 == "C" && a2 == "G") || (a1 == "G" && a2 == "C")) next

    if (n == "" || n == "NA" || n == "NaN" || n <= 0) next
    if (beta == "" || beta == "NA" || beta == "NaN") next
    if (se == "" || se == "NA" || se == "NaN" || se <= 0) next
    if (p == "" || p == "NA" || p == "NaN") next
    if (p <= 0) p = 1e-300
    if (p > 1) next

    if (maf != "NA" && maf < MIN_MAF) next
    if (MIN_INFO != "NA" && info != "NA" && info < MIN_INFO) next

    key = chr ":" bp
    if (seen[key]++) next

    print chr, bp, a2, a1, n, beta, se, maf, p, info
  }' | gzip -c > "$TMP_OUT"

  PIPE_STATUSES=("${PIPESTATUS[@]}")
  READ_STATUS=${PIPE_STATUSES[0]:-999}
  AWK_STATUS=${PIPE_STATUSES[1]:-999}
  GZIP_STATUS=${PIPE_STATUSES[2]:-999}

  if [[ ! -s "$TMP_OUT" ]]; then
    echo "WARNING: temporary output was not created or is empty: $TMP_OUT" >&2
    echo "WARNING: No final output was written. Check warnings above." >&2
    rm -f "$TMP_OUT"
  else
    ROWS=$(zcat "$TMP_OUT" 2>/dev/null | wc -l)
    if [[ "$ROWS" -lt 2 ]]; then
      echo "WARNING: formatted file has fewer than 2 rows ($ROWS). Not moving to final output." >&2
      echo "WARNING: Check column names and N_GWAS. Temporary file kept at: $TMP_OUT" >&2
    else
      mv "$TMP_OUT" "$OUT"
      echo "Created: $OUT"
      echo "Rows including header:"
      zcat "$OUT" | wc -l
      echo "Preview:"
      zcat "$OUT" | head
    fi
  fi

  if [[ "$READ_STATUS" -ne 0 || "$AWK_STATUS" -ne 0 || "$GZIP_STATUS" -ne 0 ]]; then
    echo "WARNING: formatting command returned non-zero status: read=$READ_STATUS awk=$AWK_STATUS gzip=$GZIP_STATUS" >&2
  fi
fi
