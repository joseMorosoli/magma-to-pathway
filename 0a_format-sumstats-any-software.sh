#!/bin/bash

# ============================================================================
# Format GWAS summary statistics for LDpred2 + PRSice/PRSet
#
# Output:
#   RSID CHR BP A1 A2 N BETA SE MAF P INFO
#
# Handles the current project formats:
#
# Height:
#   SNPID RSID CHR POS EFFECT_ALLELE OTHER_ALLELE
#   EFFECT_ALLELE_FREQ BETA SE P N
#
# BMI:
#   CHR POS SNP Tested_Allele Other_Allele Freq_Tested_Allele
#   BETA SE P N INFO
#   where SNP may look like rs140052487:C:A
#
# F4:
#   SNP CHR BP MAF A1 A2 BETA SE P Q_P
#   Requires N_GWAS because no per-SNP N is supplied.
#
# This formatter:
#   - recognises alternative column names automatically
#   - detects tab/comma/whitespace delimiter
#   - retains/extracts genuine rsIDs
#   - converts EAF to MAF when necessary
#   - applies MAF and INFO QC when available
#   - removes A/T and C/G variants
#   - removes duplicate rsIDs and duplicate positions
#
# Intended for continuous-trait GWAS.
# ============================================================================


# ============================================================================
# 1. USER SETTINGS
# Keep only ONE trait active
# ============================================================================


# ---- BMI --------------------------------------------------------------------
#IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/BMI_GIANT_UKB_2018_all_sites.txt.gz"
#OUT_NAME="GIANT_UKBB_BMI_2018_ALL_SITES"
#N_GWAS="NA"       # per-SNP N exists

# ---- Height -----------------------------------------------------------------
#IN="/myriadfs/home/ucju659/SUMSTATS/GIANT/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
#OUT_NAME="GIANT_HEIGHT_YENGO_2022_EUR"
#N_GWAS="NA"         # per-SNP N exists

# ---- F4 ---------------------------------------------------------------------
IN="/myriadfs/home/ucju659/SUMSTATS/PFACTOR_2025/F4_Internalizing_2025.tsv"
OUT_NAME="F4_Internalizing_2025"
N_GWAS="1637337"


# ============================================================================
# General settings
# ============================================================================

OUTDIR="/myriadfs/home/ucju659/SUMSTATS/software-ready"

OUT="${OUTDIR}/${OUT_NAME}.cleaned.tsv.gz"

MIN_MAF="0.01"
MIN_INFO="0.8"

mkdir -p "$OUTDIR"


# ============================================================================
# Check input and determine file type
# ============================================================================

CAN_RUN=true

if [[ ! -s "$IN" ]]; then
    echo "WARNING: input file not found or empty:"
    echo "$IN"
    CAN_RUN=false
fi


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
# 2. Format summary statistics
# ============================================================================

if [[ "$CAN_RUN" == true ]]; then

    echo
    echo "Input:  $IN"
    echo "Output: $OUT"

    TMP="${OUT}.tmp.$$"


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
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", x)
        return x
    }


    function findcol(names, n, i, a, key) {

        n = split(names, a, "|")

        for (i = 1; i <= n; i++) {

            key = toupper(a[i])

            if (key in h)
                return h[key]
        }

        return 0
    }


    function missing_value(x) {
    return (x == "" || x == "NA" || x == "NaN" || x == "NAN" || x == ".")
}


    # ------------------------------------------------------------------------
    # Header
    # ------------------------------------------------------------------------

    NR == 1 {

        for (i = 1; i <= NF; i++)
            h[toupper(clean($i))] = i


        # Variant identifier
        # RSID is preferred when explicitly present.
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


        # Canonical output
        print "RSID","CHR","BP","A1","A2","N","BETA","SE","MAF","P","INFO"

        next
    }


    # ------------------------------------------------------------------------
    # Data rows
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


        # --------------------------------------------------------------------
        # Chromosome
        # ------------------------------------------------------------------------

        sub(/^[Cc][Hh][Rr]/, "", chr)


        # --------------------------------------------------------------------
        # rsID
        #
        # Handles both:
        #   rs12345
        #   rs12345:A:G
        # ------------------------------------------------------------------------

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
        # Prefer supplied MAF.
        # Otherwise convert EAF -> MAF.
        # ------------------------------------------------------------------------

        maf = "NA"


        if (cMAF && !missing_value(clean($cMAF))) {

            maf = clean($cMAF) + 0

        } else if (cEAF && !missing_value(clean($cEAF))) {

            maf = clean($cEAF) + 0
        }


        if (maf != "NA" && maf > 0.5)
            maf = 1 - maf


        # --------------------------------------------------------------------
        # INFO
        # ------------------------------------------------------------------------

        if (cINFO && !missing_value(clean($cINFO)))
            info = clean($cINFO)
        else
            info = "NA"


        # --------------------------------------------------------------------
        # Structural QC
        # ------------------------------------------------------------------------

        if (chr !~ /^[0-9]+$/ || chr < 1 || chr > 22)
            next


        if (bp !~ /^[0-9]+$/ || bp <= 0)
            next


        if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/)
            next


        if (a1 == a2)
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

        if (missing_value(n) || n <= 0)
            next


        if (missing_value(beta))
            next


        if (missing_value(se) || se <= 0)
            next


        if (missing_value(p))
            next


        if (p <= 0)
            p = 1e-300


        if (p > 1)
            next


        # MAF QC when frequency exists
        if (maf != "NA" && (maf < MIN_MAF || maf > 0.5))
            next


        # INFO QC only when INFO exists
        if (MIN_INFO != "NA" && info != "NA" && info < MIN_INFO)
            next


        # --------------------------------------------------------------------
        # Remove duplicates
        # ------------------------------------------------------------------------

        if (seen_rsid[rsid]++)
            next


        poskey = chr ":" bp

        if (seen_pos[poskey]++)
            next


        # --------------------------------------------------------------------
        # Output
        # ------------------------------------------------------------------------

        print rsid,chr,bp,a1,a2,n,beta,se,maf,p,info
    }

    ' | gzip -c > "$TMP"


    STATUS=("${PIPESTATUS[@]}")


    # ------------------------------------------------------------------------
    # Validate output
    # ------------------------------------------------------------------------

    if [[ "${STATUS[0]}" -eq 0 &&
          "${STATUS[1]}" -eq 0 &&
          "${STATUS[2]}" -eq 0 &&
          -s "$TMP" &&
          "$(zcat "$TMP" | head -2 | wc -l)" -ge 2 ]]; then

        mv "$TMP" "$OUT"

        echo
        echo "Created cleaned summary statistics:"
        echo "$OUT"

        echo
        echo "Preview:"
        zcat "$OUT" | head

        echo
        echo "Rows:"
        zcat "$OUT" | wc -l

    else

        echo
        echo "WARNING: summary-statistic formatting failed."

        rm -f "$TMP"
    fi
fi


# ============================================================================
# 3. BASIC CHECKS
# ============================================================================

if [[ -s "$OUT" ]]; then

    echo
    echo "=== BASIC QC ==="

    echo "Raw rows:"
    $READCMD "$IN" | wc -l

    echo "Cleaned rows:"
    zcat "$OUT" | wc -l


    echo "Bad rsIDs:"
    zcat "$OUT" | awk '
        NR > 1 && $1 !~ /^rs[0-9]+$/ {n++}
        END {print n+0}
    '


    echo "Duplicate rsIDs:"
    zcat "$OUT" | awk '
        NR > 1 {print $1}
    ' | sort | uniq -d | wc -l

fi

zcat "$IN" | wc -l
zcat "$OUT" | wc -l

# How many raw variants have valid rsIDs?
zcat "$IN" | awk -F'\t' 'NR>1 && $2 ~ /^rs[0-9]+$/ {n++} END{print n+0}'

# How many cleaned variants?
zcat "$OUT" | awk 'NR>1 {n++} END{print n+0}'

# Final check for the ambiguity filter:

zcat "$IN" | awk -F'\t' '
NR>1 {
    a1=toupper($5); a2=toupper($6)
    if ((a1=="A"&&a2=="T")||(a1=="T"&&a2=="A")||
        (a1=="C"&&a2=="G")||(a1=="G"&&a2=="C")) n++
}
END{print "Ambiguous SNPs:",n+0}'




# ============================================================================
# 4. OPTIONAL:  Remove INFO column if it exists and contains no non-missing values
# ============================================================================

cd SUMSTATS/software-ready

zcat GIANT_UKBB_BMI_2018_ALL_SITES.cleaned.tsv.gz \
  | awk 'BEGIN{OFS="\t"} {NF--; print}' \
  | gzip -c > GIANT_UKBB_BMI_2018_ALL_SITES.cleaned.noINFO.tsv.gz


