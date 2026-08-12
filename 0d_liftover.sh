#!/bin/bash

# ============================================================
# Convert MCS PLINK SNP coordinates: GRCh38 -> GRCh37/hg19
# Uses rsID mapping from UCSC dbSNP155
# Original GRCh38 files are NOT modified
# ============================================================

IN="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only"

OUTDIR="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS"
OUT="${OUTDIR}/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only_GRCh37"

TOOLS="/myriadfs/home/ucju659/SOFTWARE/ucsc"
TMP="${OUTDIR}/liftover_tmp"

BBNI="${TOOLS}/bigBedNamedItems"
DBSNP="${TOOLS}/dbSnp155Common.hg19.bb"

CONTAINER_DIR="/myriadfs/home/ucju659/SOFTWARE/containers"
CONTAINER="${CONTAINER_DIR}/ubuntu24.sif"

mkdir -p "${OUTDIR}" "${TOOLS}" "${TMP}" "${CONTAINER_DIR}"


# ------------------------------------------------------------
# 1. Check PLINK
# ------------------------------------------------------------

module load plink/1.90b3.40

if ! command -v plink >/dev/null 2>&1; then
    echo "ERROR: PLINK is not available."
    return 1
fi


# ------------------------------------------------------------
# 2. Check required UCSC files
# ------------------------------------------------------------

if [[ ! -x "${BBNI}" ]]; then
    echo "ERROR: bigBedNamedItems not found or not executable:"
    echo "${BBNI}"
    return 1
fi

if [[ ! -s "${DBSNP}" ]]; then
    echo "ERROR: dbSNP hg19 file not found:"
    echo "${DBSNP}"
    return 1
fi


# ------------------------------------------------------------
# 3. Extract rsIDs from MCS BIM
# ------------------------------------------------------------

awk '$2 ~ /^rs[0-9]+$/ {print $2}' "${IN}.bim" \
    | sort -u \
    > "${TMP}/mcs_rsids.txt"

echo "Variants in original BIM:"
wc -l "${IN}.bim"

echo "Unique rsIDs:"
wc -l "${TMP}/mcs_rsids.txt"

echo "Preview of rsIDs:"
head "${TMP}/mcs_rsids.txt"


# ------------------------------------------------------------
# 4. Load Apptainer and obtain container if required
# ------------------------------------------------------------

module load apptainer

if [[ ! -s "${CONTAINER}" ]]; then
    echo "Ubuntu container not found; downloading..."
    apptainer pull "${CONTAINER}" docker://ubuntu:24.04

    if [[ ! -s "${CONTAINER}" ]]; then
        echo "ERROR: Failed to create Apptainer container."
        return 1
    fi
fi


# ------------------------------------------------------------
# 5. Check that paths are visible inside container
# ------------------------------------------------------------

echo "Checking files from inside Apptainer..."

apptainer exec "${CONTAINER}" \
    ls -lh "${TMP}/mcs_rsids.txt"

apptainer exec "${CONTAINER}" \
    ls -lh "${DBSNP}"


# ------------------------------------------------------------
# 6. Look up each rsID directly in hg19 dbSNP
# ------------------------------------------------------------

apptainer exec "${CONTAINER}" \
    "${BBNI}" -nameFile \
    "${DBSNP}" \
    "${TMP}/mcs_rsids.txt" \
    "${TMP}/dbsnp_hg19_raw.bed"


# ------------------------------------------------------------
# 7. Check output
# ------------------------------------------------------------

if [[ ! -s "${TMP}/dbsnp_hg19_raw.bed" ]]; then
    echo "ERROR: dbSNP lookup produced no output."
    return 1
fi

echo
echo "dbSNP lookup successful:"
wc -l "${TMP}/dbsnp_hg19_raw.bed"

echo "Preview:"
head "${TMP}/dbsnp_hg19_raw.bed"

# ------------------------------------------------------------
# 8. Reduce dbSNP output to rsID, chromosome and GRCh37 position
# ------------------------------------------------------------

awk 'BEGIN {OFS="\t"}
{
    chr=$1
    sub(/^chr/, "", chr)

    # Keep autosomal, single-base variants only.
    # UCSC BED start is 0-based; PLINK position is 1-based.
    if (chr ~ /^([1-9]|1[0-9]|2[0-2])$/ && ($3-$2)==1)
        print $4, chr, $2+1
}' "${TMP}/dbsnp_hg19_raw.bed" \
| sort -u \
> "${TMP}/hg19_candidates.txt"

echo
echo "Mapping rows:"
wc -l "${TMP}/hg19_candidates.txt"

echo "Unique rsIDs represented:"
cut -f1 "${TMP}/hg19_candidates.txt" | sort -u | wc -l

echo "Original MCS rsIDs:"
wc -l "${TMP}/mcs_rsids.txt"


# ------------------------------------------------------------
# 9. Keep only rsIDs with exactly one GRCh37 mapping
# ------------------------------------------------------------

awk '
{
    n[$1]++
    row[$1]=$0
}
END {
    for (id in n)
        if (n[id] == 1)
            print row[id]
}' "${TMP}/hg19_candidates.txt" \
| sort -k1,1 \
> "${TMP}/hg19_unique_mapping.txt"

echo
echo "Unique unambiguous GRCh37 mappings:"
wc -l "${TMP}/hg19_unique_mapping.txt"

echo "Preview:"
head "${TMP}/hg19_unique_mapping.txt"

echo
echo "Check example rs12562034:"
grep -w "rs12562034" "${TMP}/hg19_unique_mapping.txt" || true


# ------------------------------------------------------------
# 10. Create PLINK update files
# ------------------------------------------------------------

# Format:
# rsID    new chromosome
awk 'BEGIN {OFS="\t"} {print $1,$2}' \
    "${TMP}/hg19_unique_mapping.txt" \
    > "${TMP}/update_chr.txt"

# Format:
# rsID    new GRCh37 base-pair position
awk 'BEGIN {OFS="\t"} {print $1,$3}' \
    "${TMP}/hg19_unique_mapping.txt" \
    > "${TMP}/update_bp.txt"

# rsIDs that can safely be retained
cut -f1 "${TMP}/hg19_unique_mapping.txt" \
    > "${TMP}/mapped_rsids.txt"


# ------------------------------------------------------------
# 11. Record rsIDs that could not be uniquely mapped
# ------------------------------------------------------------

comm -23 \
    <(sort "${TMP}/mcs_rsids.txt") \
    <(sort "${TMP}/mapped_rsids.txt") \
    > "${TMP}/unmapped_rsids.txt"

echo "Mapped rsIDs:"
wc -l "${TMP}/mapped_rsids.txt"

echo "Unmapped or ambiguous rsIDs:"
wc -l "${TMP}/unmapped_rsids.txt"


# ------------------------------------------------------------
# 12. Create the new GRCh37 PLINK dataset
# ------------------------------------------------------------

plink \
    --bfile "${IN}" \
    --extract "${TMP}/mapped_rsids.txt" \
    --update-chr "${TMP}/update_chr.txt" \
    --update-map "${TMP}/update_bp.txt" \
    --keep-allele-order \
    --make-bed \
    --out "${OUT}"

if [[ ! -s "${OUT}.bim" ]]; then
    echo "ERROR: GRCh37 PLINK dataset was not created."
    return 1
fi


# ------------------------------------------------------------
# 13. Final QC
# ------------------------------------------------------------

echo "============================================"
echo "GRCh38 -> GRCh37 conversion complete"
echo "============================================"

echo "Original variants:"
wc -l "${IN}.bim"

echo "GRCh37 variants:"
wc -l "${OUT}.bim"

echo "Original first five variants:"
head -5 "${IN}.bim"

echo "GRCh37 first five variants:"
head -5 "${OUT}.bim"

echo "Check rs12562034 before and after:"
grep -w "rs12562034" "${IN}.bim" || true
grep -w "rs12562034" "${OUT}.bim" || true

echo "Output files:"
echo "${OUT}.bed"
echo "${OUT}.bim"
echo "${OUT}.fam"

# Checks
IN="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only"
OUT="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only_GRCh37"

wc -l "${IN}.bim" "${OUT}.bim"
wc -l "${IN}.fam" "${OUT}.fam"

head "${OUT}.bim"

grep -w "rs12562034" "${IN}.bim"
grep -w "rs12562034" "${OUT}.bim"
