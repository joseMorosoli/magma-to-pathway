TARGET_BIM="/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS/MCS_topmed_EUR_KING_QCd_rsID_PCs_SD4-hapmap-only.bim"

#IN="/myriadfs/home/ucju659/SUMSTATS/ldpred2_ready/GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz"
IN="/myriadfs/home/ucju659/SUMSTATS/ldpred2_ready/GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz"

#OUT="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_HEIGHT_YENGO_2022_EUR.rsID.prsice.tsv.gz"
OUT="/myriadfs/home/ucju659/SUMSTATS/prsice_ready/GIANT_UKBB_BMI_2018_ALL_SITES.rsID.prsice.tsv.gz"

mkdir -p /myriadfs/home/ucju659/SUMSTATS/prsice_ready

TMP_MAP="/tmp/target_chrpos_to_snp.rsIDtarget.${USER}.$$"

# .bim columns:
# 1 = CHR
# 2 = SNP ID, here rsID
# 3 = genetic distance
# 4 = BP
# 5 = allele 1
# 6 = allele 2
#
# This maps CHR:BP from the GWAS file to the rsID used in the target .bim.
awk 'BEGIN{OFS="\t"} {print $1":"$4, $2}' "$TARGET_BIM" > "$TMP_MAP"

awk '
BEGIN {
  OFS="\t"
}
NR==FNR {
  snp[$1] = $2
  next
}
FNR==1 {
  print "MarkerName","CHR","POS","A1","A2","N","BETA","SE","MAF","P"
  next
}
{
  # Input columns:
  # CHR BP A2 A1 N BETA SE MAF P INFO
  key = $1":"$2

  if (key in snp) {
    print snp[key], $1, $2, $4, $3, $5, $6, $7, $8, $9
    kept++
  } else {
    missed++
  }
}
END {
  print "rsID target - Kept SNPs: " kept+0 > "/dev/stderr"
  print "rsID target - Missed SNPs: " missed+0 > "/dev/stderr"
}
' "$TMP_MAP" <(zcat "$IN") | gzip -c > "$OUT"

rm -f "$TMP_MAP"

zcat "$OUT" | head
