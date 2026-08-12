#!/bin/bash

# ============================================================================
# Format GWAS summary statistics for LDpred2 and PRSice
#
# LDpred2 output:
#   CHR BP A2 A1 N BETA SE MAF P INFO RSID
#
# PRSice output:
#   RSID CHR POS A1 A2 N BETA SE MAF P INFO
#
# Features:
#   - Recognises multiple common names for the same column.
#   - Automatically detects tab, comma, or whitespace delimiters.
#   - Handles both:
#         rs12345
#         rs12345:A:G
#   - Converts EAF to MAF when necessary.
#   - Uses per-SNP N when available; otherwise N_GWAS.
#   - Removes ambiguous A/T and C/G SNPs.
#   - Removes duplicate rsIDs and duplicate positions.
# ============================================================================


# ============================================================================
# USER SETTINGS
# Keep only ONE trait active
# ============================================================================


# ---- BMI --------------------------------------------------------------------
# IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/BMI_GIANT_UKB_2018_all_sites.txt.gz"
# OUT_NAME="GIANT_UKBB_BMI_2018_ALL_SITES"
# N_GWAS="NA"       # per-SNP N exists


# ---- Height -----------------------------------------------------------------

IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
OUT_NAME="GIANT_HEIGHT_YENGO_2022_EUR"
N_GWAS="NA"         # per-SNP N exists


# ---- F4 ---------------------------------------------------------------------
# IN="/myriadfs/home/ucju659/SUMSTATS/PFACTOR_2025/F4_Internalizing_2025.tsv"
# OUT_NAME="F4_Internalizing_2025"
# N_GWAS="1637337"


# ============================================================================
# General settings
# ============================================================================

OUTDIR="/myriadfs/home/ucju659/SUMSTATS/software-ready"

OUT_LDPRED2="${OUTDIR}/${OUT_NAME}.ldpred2.tsv.gz"
OUT_PRSICE="${OUTDIR}/${OUT_NAME}.prsice.tsv.gz"

MIN_MAF="0.01"
MIN_INFO="0.8"

mkdir -p "$OUTDIR"


# ============================================================================
# Input checks
# ============================================================================

CAN_RUN=true
LDPRED2_OK=false

if [[ ! -s "$IN" ]]; then
    echo "WARNING: input file not found or empty:"
    echo "$IN"
    CAN_RUN=false
fi


# gzip or plain text
if [[ "$IN" == *.gz ]]; then
    READCMD="zcat"
else
    READCMD="cat"
fi


# ============================================================================
# Detect delimiter
# ============================================================================

if [[ "$CAN_RUN" == true ]]; then

    HEADER=$($READCMD "$IN" | head -1)

    if [[ "$HEADER" == *$'\t'* ]]; then
        AWK_FS='\t'
        echo "Detected delimiter: TAB"

    elif [[ "$HEADER" == *,* ]]; then
        AWK_FS=','
        echo "Detected delimiter: COMMA"

    else
        AWK_FS='[[:space:]]+'
        echo "Detected delimiter: WHITESPACE"
    fi
fi


# ============================================================================
# 1. Raw summary statistics -> LDpred2
# ============================================================================

if [[ "$CAN_RUN" == true ]]; then

    echo
    echo "Input:          $IN"
    echo "LDpred2 output: $OUT_LDPRED2"
    echo "PRSice output:  $OUT_PRSICE"

    TMP_LDPRED2="${OUT_LDPRED2}.tmp.$$"


    $READCMD "$IN" | awk \
        -v FS="$AWK_FS" \
        -v OFS="\t" \
        -v N_FALLBACK="$N_GWAS" \
        -v MIN_MAF="$MIN_MAF" \
        -v MIN_INFO="$MIN_INFO" '

    # ------------------------------------------------------------------------
    # Functions
    # ------------------------------------------------------------------------

    function clean(x) {
        gsub(/\r/, "", x)
        gsub(/^ +| +$/, "", x)
        return x
    }


    function findcol(names, n,i,a,key) {

        n = split(names, a, "|")

        for (i = 1; i <= n; i++) {

            key = toupper(a[i])

            if (key in h)
                return h[key]
        }

        return 0
    }


    # ------------------------------------------------------------------------
    # Header
    # ------------------------------------------------------------------------

    NR == 1 {

        for (i = 1; i <= NF; i++)
            h[toupper(clean($i))] = i


        # Variant identifier
        cRSID = findcol("RSID|RS_ID|RS_NUMBER|SNP|MARKERNAME|SNPID")

        # Chromosome / position
        cCHR = findcol("CHR|CHROM|CHROMOSOME")
        cBP  = findcol("BP|POS|POSITION|BASE_PAIR_LOCATION")

        # Alleles
        cA1 = findcol("A1|EFFECT_ALLELE|EA|ALLELE1|TESTED_ALLELE")
        cA2 = findcol("A2|OTHER_ALLELE|NON_EFFECT_ALLELE|NEA|ALLELE2")

        # Association statistics
        cBETA = findcol("BETA|EFFECT|B")
        cSE   = findcol("SE|STDERR|STANDARD_ERROR")
        cP    = findcol("P|PVAL|PVALUE|P_VALUE")

        # Sample size
        cN = findcol("N|TOTAL_N|NMISS|N_TOTAL")

        # Allele frequency
        cMAF = findcol("MAF|MINOR_ALLELE_FREQ|MINOR_ALLELE_FREQUENCY")
        cEAF = findcol("EAF|EAF_HRC|EFFECT_ALLELE_FREQ|EFFECT_ALLELE_FREQUENCY|AF1|FREQ_TESTED_ALLELE")

        # Imputation quality
        cINFO = findcol("INFO|INFO_SCORE|IMPINFO|RSQ")


        # --------------------------------------------------------------------
        # Required columns
        # --------------------------------------------------------------------

        missing = ""

        if (!cRSID) missing = missing " rsID/SNP"
        if (!cCHR)  missing = missing " CHR"
        if (!cBP)   missing = missing " BP/POS"
        if (!cA1)   missing = missing " A1/effect allele"
        if (!cA2)   missing = missing " A2/other allele"
        if (!cBETA) missing = missing " BETA"
        if (!cSE)   missing = missing " SE"
        if (!cP)    missing = missing " P"

        if (!cN && (N_FALLBACK == "" || N_FALLBACK == "NA" || N_FALLBACK == "PUT_REAL_F4_N_HERE"))
            missing = missing " N"


        if (missing != "") {

            print "WARNING: missing required column(s):" missing > "/dev/stderr"

            skip_all = 1
            next
        }


        print "CHR", "BP", "A2", "A1", "N", \
              "BETA", "SE", "MAF", "P", "INFO", "RSID"

        next
    }


    # ------------------------------------------------------------------------
    # Variants
    # ------------------------------------------------------------------------

    NR > 1 {

        if (skip_all)
            next


        rsid = clean($cRSID)

        chr = clean($cCHR)
        bp  = clean($cBP)

        a1 = toupper(clean($cA1))
        a2 = toupper(clean($cA2))

        beta = clean($cBETA)
        se   = clean($cSE)
        p    = clean($cP)


        # Remove chr prefix if present
        sub(/^[Cc][Hh][Rr]/, "", chr)


        # --------------------------------------------------------------------
        # rsID
        #
        # Handles:
        #   rs12345
        #   rs12345:A:G
        # --------------------------------------------------------------------

        if (rsid ~ /^rs[0-9]+:/)
            sub(/:.*/, "", rsid)

        if (rsid !~ /^rs[0-9]+$/)
            next


        # --------------------------------------------------------------------
        # Sample size
        # ------------------------------------------------------------------------

        if (cN)
            n = clean($cN)
        else
            n = N_FALLBACK


        # --------------------------------------------------------------------
        # MAF
        #
        # Prefer MAF if supplied.
        # Otherwise convert EAF -> MAF.
        # ------------------------------------------------------------------------

        maf = "NA"

        if (cMAF && $cMAF != "" && $cMAF != "NA" && $cMAF != "NaN") {

            maf = $cMAF + 0

        } else if (cEAF && $cEAF != "" && $cEAF != "NA" && $cEAF != "NaN") {

            maf = $cEAF + 0
        }


        if (maf != "NA" && maf > 0.5)
            maf = 1 - maf


        # --------------------------------------------------------------------
        # INFO
        # ------------------------------------------------------------------------

        if (cINFO && $cINFO != "" && $cINFO != "NA" && $cINFO != "NaN")
            info = $cINFO
        else
            info = "NA"


        # --------------------------------------------------------------------
        # Basic QC
        # ------------------------------------------------------------------------

        if (chr !~ /^[0-9]+$/ || chr < 1 || chr > 22)
            next

        if (bp !~ /^[0-9]+$/)
            next

        if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/ || a1 == a2)
            next


        # Remove strand-ambiguous SNPs
        if ((a1 == "A" && a2 == "T") ||
            (a1 == "T" && a2 == "A") ||
            (a1 == "C" && a2 == "G") ||
            (a1 == "G" && a2 == "C"))
            next


        # --------------------------------------------------------------------
        # Numeric QC
        # ------------------------------------------------------------------------

        if (n == "" || n == "NA" || n == "NaN" || n <= 0)
            next

        if (beta == "" || beta == "NA" || beta == "NaN")
            next

        if (se == "" || se == "NA" || se == "NaN" || se <= 0)
            next

        if (p == "" || p == "NA" || p == "NaN")
            next

        if (p <= 0)
            p = 1e-300

        if (p > 1)
            next


        # MAF filter
        if (maf != "NA" && maf < MIN_MAF)
            next


        # INFO filter only when INFO exists
        if (MIN_INFO != "NA" && info != "NA" && info < MIN_INFO)
            next


        # --------------------------------------------------------------------
        # Remove duplicate rsIDs and positions
        # ------------------------------------------------------------------------

        if (seen_rsid[rsid]++)
            next

        poskey = chr ":" bp

        if (seen_pos[poskey]++)
            next


        # --------------------------------------------------------------------
        # LDpred2 output
        # ------------------------------------------------------------------------

        print chr, bp, a2, a1, n, \
              beta, se, maf, p, info, rsid
    }

    ' | gzip -c > "$TMP_LDPRED2"


    STATUS=("${PIPESTATUS[@]}")


    # ------------------------------------------------------------------------
    # Check output
    # ------------------------------------------------------------------------

    if [[ "${STATUS[0]}" -eq 0 &&
          "${STATUS[1]}" -eq 0 &&
          "${STATUS[2]}" -eq 0 &&
          -s "$TMP_LDPRED2" &&
          "$(zcat "$TMP_LDPRED2" | head -2 | wc -l)" -ge 2 ]]; then

        mv "$TMP_LDPRED2" "$OUT_LDPRED2"

        LDPRED2_OK=true

        echo
        echo "Created LDpred2 file:"
        echo "$OUT_LDPRED2"

        echo "Preview:"
        zcat "$OUT_LDPRED2" | head

        echo "Rows:"
        zcat "$OUT_LDPRED2" | wc -l

    else

        echo "WARNING: LDpred2 formatting failed."
        rm -f "$TMP_LDPRED2"
    fi
fi


# ============================================================================
# 2. LDpred2 -> PRSice
# ============================================================================

if [[ "$LDPRED2_OK" == true ]]; then

    TMP_PRSICE="${OUT_PRSICE}.tmp.$$"


    zcat "$OUT_LDPRED2" | awk '

    BEGIN {
        OFS="\t"
    }

    NR == 1 {
        print "RSID","CHR","POS","A1","A2","N","BETA","SE","MAF","P","INFO"
        next
    }

    {
        print $11,$1,$2,$4,$3,$5,$6,$7,$8,$9,$10
    }

    ' | gzip -c > "$TMP_PRSICE"


    STATUS=("${PIPESTATUS[@]}")


    if [[ "${STATUS[0]}" -eq 0 &&
          "${STATUS[1]}" -eq 0 &&
          "${STATUS[2]}" -eq 0 &&
          -s "$TMP_PRSICE" &&
          "$(zcat "$TMP_PRSICE" | head -2 | wc -l)" -ge 2 ]]; then

        mv "$TMP_PRSICE" "$OUT_PRSICE"

        echo
        echo "Created PRSice file:"
        echo "$OUT_PRSICE"

        echo "Preview:"
        zcat "$OUT_PRSICE" | head

        echo "Rows:"
        zcat "$OUT_PRSICE" | wc -l

    else

        echo "WARNING: PRSice conversion failed."
        rm -f "$TMP_PRSICE"
    fi

else

    echo "WARNING: PRSice conversion skipped because LDpred2 formatting failed."
fi

# ============================================================================
# 3. CHECKS
# ============================================================================


# ---- BMI --------------------------------------------------------------------
# IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/BMI_GIANT_UKB_2018_all_sites.txt.gz"
# OUT_NAME="GIANT_UKBB_BMI_2018_ALL_SITES"


# ---- Height -----------------------------------------------------------------

IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
OUT_NAME="GIANT_HEIGHT_YENGO_2022_EUR"


# ---- F4 ---------------------------------------------------------------------
# IN="/myriadfs/home/ucju659/SUMSTATS/PFACTOR_2025/F4_Internalizing_2025.tsv"
# OUT_NAME="F4_Internalizing_2025"


OUTDIR="/myriadfs/home/ucju659/SUMSTATS/software-ready"

LDPRED2="${OUTDIR}/${OUT_NAME}.ldpred2.tsv.gz"
PRSICE="${OUTDIR}/${OUT_NAME}.prsice.tsv.gz"


# gzip or plain-text input
if [[ "$IN" == *.gz ]]; then
    READCMD="zcat"
else
    READCMD="cat"
fi


echo "=== RAW INPUT ==="
$READCMD "$IN" | head

echo "=== LDPRED2 ==="
zcat "$LDPRED2" | head

echo "=== PRSICE ==="
zcat "$PRSICE" | head

echo "=== ROW COUNTS ==="
$READCMD "$IN" | wc -l
zcat "$LDPRED2" | wc -l
zcat "$PRSICE" | wc -l

echo "=== BAD rsIDs ==="
zcat "$LDPRED2" | awk 'NR>1 && $11 !~ /^rs[0-9]+$/ {n++} END {print "LDpred2:", n+0}'
zcat "$PRSICE"  | awk 'NR>1 && $1  !~ /^rs[0-9]+$/ {n++} END {print "PRSice:", n+0}'

echo "=== DUPLICATE rsIDs ==="
zcat "$LDPRED2" | awk 'NR>1 {print $11}' | sort | uniq -d | wc -l
zcat "$PRSICE"  | awk 'NR>1 {print $1}'  | sort | uniq -d | wc -l





