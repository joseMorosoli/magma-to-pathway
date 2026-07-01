
# ==============================================================================
# Custom LDpred2 Pipeline for UCL Myriad
# ==============================================================================
#
# Authors:
#   Andrea G. Allegrini, Jose J. Morosoli
#
# Date created:
#   17 June 2026
#
# Description:
#   Custom R pipeline for generating LDpred2 polygenic scores on the
#   UCL Myriad computing cluster.
#
# Source:
#   Developed by José J. Morosoli based on the original LDpred2 pipeline
#   created by A. Allegrini (https://github.com/AndreAllegrini/LDpred2).
#
# Computing environment:
#   UCL Myriad
#   https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/
#
# Notes:
#   This is a project-specific adaptation and should be validated when used
#   with new GWAS summary statistics, genotype datasets, LD references, genome
#   builds, or software versions.
#
# Required input:
#   - LDpred2-ready GWAS summary statistics produced by script 0.
#   - bigsnpr .rds genotype object for the target sample.
#   - HapMap3+ LD reference directory containing map_hm3_plus.rds and LDref/.
#
# Main output:
#   - *_pred_inf.txt
#   - *_pred_auto.txt
#   - *_beta_inf.txt
#   - *_final_beta_auto.txt
#   - *_multi_beta_auto.rds
#   - *_auto_chain_diagnostics.txt
#   - run log and diagnostic plots
#
# Key project-specific changes retained here:
#   - explicit input checks before LDpred2-auto;
#   - p_init restricted to 1e-4 to 0.2;
#   - shrink_corr = 0.95;
#   - use_MLE = FALSE;
#   - allow_jump_sign = FALSE;
#   - diagnostic saving before convergence filtering.
#
# Do not commit generated scores, sparse backing files, logs or large RDS files.
# ==============================================================================

#! /usr/bin/env Rscript

required_packages <- c(
  "optparse", "ggplot2", "bigsnpr", "bigreadr", "data.table",
  "dplyr", "tidyr", "tibble", "vctrs", "cowplot", "bigstatsr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "
Install these in the Myriad R library used by R_LIBS before running.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(bigsnpr)
  library(bigreadr)
  library(data.table)
  library(cowplot)
  library(bigstatsr)
})

# this scripts uses info from the following resources:
#recommended GWAS QC
#see  https://github.com/privefl/paper-misspec/tree/main/code/prepare-sumstats
#see https://privefl.github.io/bigsnpr/articles/LDpred2.html
# ew LD ref - hap-map + 
#see also https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R


option_list = list(
  make_option(c("-s", "--sumstats"), type="character", default=NULL, 
              help="Name of GWAS summary statistics.\n
                Note sumstats should have *at least* the following header:\n
                case/control traits: CHR BP A2 A1 NCAS NCON BETA SE \n
                continuous traits: CHR BP A2 A1 N BETA SE\n", metavar="character"),
  
  make_option(c("-g", "--geno"), type="character",
              help="path/to/bigsnp.rds\n", metavar="character"),
  
  make_option(c("-t", "--type"), type="logical", 
              default = TRUE,
              help="Whether GWAS trait is case/control\n 
              (TRUE = binary, FALSE = continuous).\n 
              [default = %default]", metavar="logical"),
  
  make_option(c("-o", "--out"), type="character", 
              help="path/to/output_dir/.\n", metavar="character"),
  
  make_option(c("-d", "--sdir"), type="character", 
              help="path/to/sumstats_dir/.\n", metavar="character"),
  
  make_option(c("-m", "--misc"), type="character", 
              help = "path/to/hapmap3plus/map_hm3_plus.rds", metavar="character"),
  
  make_option(c("--maf"), type = "double", default = 0.01,
              help = "MAF threshold for QC. [default = %default]",
              metavar = "numeric"),
  
  make_option(c("--info"), type = "double", default = 0.6,
              help = "INFO threshold for QC. [default = %default]",
              metavar = "numeric"),
  
  make_option(c("-l", "--lift"), type = "character",
              default = "hg19",
              help = "Genome build of the test (target) genotype data.
              Summary statistics will be automatically lifted to this build if necessary.
              Options: hg18, hg19, hg38, or FALSE to disable lift. [default = %default]",
              metavar = "character"),
  
  make_option(c("-c", "--cores"), type="integer", 
              default = 8, 
              help="Number of cores. \n 
              [default = %default]", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out_path = opt$out
misc_path = opt$misc
NCORES = opt$cores
geno = opt$geno

if (is.null(out_path) || is.null(misc_path) || is.null(geno) || is.null(opt$sdir) || is.null(opt$sumstats)) {
  stop("Missing required command-line arguments. Use --help for required inputs.", call. = FALSE)
}

if (!grepl("/$", out_path)) out_path <- paste0(out_path, "/")
if (!grepl("/$", misc_path)) misc_path <- paste0(misc_path, "/")
if (!grepl("/$", opt$sdir)) opt$sdir <- paste0(opt$sdir, "/")

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# checks

# Sumstats not found
sumstats_path <- file.path(opt$sdir, opt$sumstats)
if (!file.exists(sumstats_path)) {
  stop("Sumstats file not found at: ", sumstats_path, call. = FALSE)
}

lift_val <- opt$lift
if (lift_val == "FALSE") lift_val <- FALSE
if (!isFALSE(lift_val)) {
  lift_val <- tolower(lift_val)# only if not FALSE
}


if (!isFALSE(lift_val) && !lift_val %in% c("hg18","hg19","hg38")) {
  stop("--lift must be one of: 'hg18', 'hg19', 'hg38', or FALSE. You gave: ", opt$lift, call. = FALSE)
}

#type
opt$type <- isTRUE(opt$type)


start <- Sys.time() #start time

file_log <- file.path(out_path, paste0(opt$sumstats, ".log")) 

file.create(file_log) #create log file

cat(" -------------------------------------------------","\n",
    " Generate LDpred2-auto and infinitesimal scores. ","\n",
    " bugs and questions: a.allegrini@ucl.ac.uk ","\n",
    "-------------------------------------------------","\n",
    " ","\n", sep = " ", file=file_log, append=TRUE)

#print  options
sink(file = file_log, append = TRUE)
#print_help(opt_parser)

cat("#####","SELECTED OPTIONS:","#####\n"," ", "\n", file=file_log, append=TRUE)
opt
cat("##########\n"," ", "\n", file=file_log, append=TRUE)
sink()


cat(paste0("Analyses started at ", start),"\n"," ","\n", file=file_log,append=TRUE)

cat("R version: ", R.version.string, "\n",
    "bigsnpr version: ", as.character(packageVersion("bigsnpr")), "\n",
    "bigstatsr version: ", as.character(packageVersion("bigstatsr")), "\n",
    file = file_log, append = TRUE, sep = "")


cat("Reading: ", opt$sumstats, " sumstats.", "\n"," ","\n", sep='',file=file_log,append=TRUE) 

# load sumstats and convert to LDpred header format + calculate effective sample size 

sumstats <- bigreadr::fread2(input = sumstats_path)

#sumstats <- bigreadr::fread2(input = paste0(sumstatDir,opt$s))

cat("Loaded sumstats have: ", dim(sumstats)[1], " rows and ", dim(sumstats)[2], " columns.", "\n"," ","\n", sep="",file=file_log,append=TRUE) 


cat("Sumstats header is: ", paste(names(sumstats), collapse = " ", sep = " "), "\n"," ","\n", sep="",file=file_log,append=TRUE) 

# basic QC 

#to add discard sample N whe .6 < .9 quant 
#aslo add per variant effective sample size if more than .5 discarded
#https://github.com/privefl/bigsnpr/issues/281


cat("Starting sumstats QC...","\n"," ","\n", file=file_log,append=TRUE)  


# filter out duplicates
dups <- vctrs::vec_duplicate_detect(sumstats[, c("CHR", "BP")])

sumstats <- sumstats[!dups, ]

#dup_count <- nrow(sumstats[duplicated(sumstats[, c("chr","pos")]), ])


cat("There were: ", sum(dups)," duplicated physical positions in GWAS data.\n" ," ","\n", sep='',file=file_log,append=TRUE)


# filter on MAF if present 
if ("MAF" %in% colnames(sumstats)) {
  
  # Ensure MAF = 0.5
  sumstats$MAF <- ifelse(sumstats$MAF <= 0.5, sumstats$MAF, 1 - sumstats$MAF)
  
  n_maf_disc <- sum(sumstats$MAF < opt$maf, na.rm = TRUE)
  
  cat("N = ", n_maf_disc,
      " variants discarded because MAF < ", opt$maf, "\n",
      sep = '', file = file_log, append = TRUE)
  
  sumstats <- subset(sumstats, MAF >= opt$maf | is.na(MAF))
  
} else {
  cat("No MAF column provided.\n", file = file_log, append = TRUE)
}

# filter on INFO if present 
if ("INFO" %in% colnames(sumstats)) {
  
  n_info_disc <- sum(sumstats$INFO < opt$info, na.rm = TRUE)
  
  cat("N = ", n_info_disc,
      " variants discarded because INFO < ", opt$info, "\n",
      sep = '', file = file_log, append = TRUE)
  
  sumstats <- subset(sumstats, INFO >= opt$info | is.na(INFO))
  
} else {
  cat("No INFO column provided.\n"," ","\n", file = file_log, append = TRUE)
}


# effective sample size
if(opt$type == T){ 
  
  #check colummn names are correct
  if (!all(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for case/control traits is:\n        ",
                paste(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE"), collapse=", ")))}
  
  #select columns and rename
  sumstats <- sumstats[,c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE")]
  names(sumstats) <- c("chr", "pos", "a0", "a1", "Ncas","Ncon","beta","beta_se")
  
  Ncas <- sumstats$Ncas 
  Ncon <- sumstats$Ncon
  
  sumstats$n_eff <- 4 / (1 / Ncas + 1 / Ncon)
  
}else{ #if not
  
  #check colummn names are correct
  if (!all(c("CHR","BP","A2","A1","N","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for continuous traits is:\n        ",
                paste(c("CHR","BP","A2","A1","N","BETA","SE"), collapse=", ")))}
  
  #select columns and rename
  sumstats <- sumstats[,c("CHR","BP","A2","A1","N","BETA","SE")]
  names(sumstats) <- c("chr", "pos", "a0", "a1", "n_eff","beta","beta_se")
  
}


#sometimes this is an issue:
if(!is.integer(sumstats$chr)){
  sumstats$chr <- as.integer(sumstats$chr)
}


#also this
sumstats$a0 <- toupper(sumstats$a0)
sumstats$a1 <- toupper(sumstats$a1)

map_ldref <- readRDS(paste0(misc_path,"map_hm3_plus.rds")) #read reference map


# update build if needed using positions for hapmap3+ set provided in ref file

# check overlap with all builds

overlap_hg18 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos_hg18))
overlap_hg19 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos))
overlap_hg38 <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos_hg38))

overlaps <- c(hg18 = overlap_hg18, hg19 = overlap_hg19, hg38 = overlap_hg38)
max_build <- names(which.max(overlaps))

cat("Detected GWAS build (max overlap):", max_build, "\n"," ","\n", 
    "Overlaps: hg18 =", overlap_hg18, ", hg19 =", overlap_hg19, ", hg38 =", overlap_hg38, "\n"," ","\n",
    file = file_log, append = TRUE)



# target set build
if (isFALSE(lift_val)) {
  
  cat("lift option = FALSE: no lift performed", "\n"," ","\n",
      file = file_log, append = TRUE)
  
} else {
  map_ldref$pos_hg19 <- map_ldref$pos # create hg19 pos column for consistency - easy fix otherwise merge breaks
  target_build <- lift_val  # user specified build
  
  
  # Harmonize sumstats to target_build
  if (max_build != target_build) {
    # Merge sumstats with reference to lift positions
    
    cat("Sumstats build (",max_build ,") different from target build. Lifting to:", target_build, "\n"," ","\n", 
        file = file_log, append = TRUE)
    
    sumstats <- merge(
      sumstats,
      map_ldref[, c("chr", "pos", "pos_hg18","pos_hg19", "pos_hg38")],
      by.x = c("chr", "pos"),
      by.y = c("chr", paste0("pos_", max_build)),
      all = FALSE,
      sort = FALSE)
    
    sumstats$pos <- sumstats[[paste0("pos_", target_build)]]
    
    
    # Update map_ldref positions to target_build
    map_ldref$pos <- if (target_build == "hg19") map_ldref$pos else map_ldref[[paste0("pos_", target_build)]]
    
    post_overlap <- sum(paste0(sumstats$chr, ":", sumstats$pos) %in% paste0(map_ldref$chr, ":", map_ldref$pos))
    cat("Overlap after lift:", post_overlap, "\n", file=file_log, append=TRUE)
  } else {
    
    cat("GWAS already in target build; no lift needed\n", file=file_log, append=TRUE)
    
  }
  
}

info_snp <- tibble::as_tibble(snp_match(sumstats, map_ldref, return_flip_and_rev = TRUE))


cat("There were: ", sum( vctrs::vec_duplicate_detect(sumstats[, c("chr","pos")]))," duplicated physical positions in GWAS data.\n" ," ","\n", sep='',file=file_log,append=TRUE)

# match sumstats with reference data
cat(" N = ", dim(info_snp)[1], "SNPs have been matched with reference data (i.e. HapMap3 + )\n",
    "N =  ",sum(info_snp$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(info_snp$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

#drop NAs and make sure order is the same with SD file below
info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))

#chi-squared GWAS
chi2 <- with(info_snp, (beta / beta_se)^2)
cat("Mean Chi^2 = ",mean(chi2,na.rm=T),".\n",sep='',file=file_log,append=TRUE)  


#lambda gc 
lgc <- median(chi2,na.rm=T)/qchisq(0.5,1)
cat("Lambda GC = ",lgc,".\n",sep='',file=file_log,append=TRUE)  

#SD REFERENCE

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

if(opt$type == T){  #if GWAS trait is binary 
  
  cat("treating GWAS trait as binary.\n",sep='',file=file_log,append=TRUE)  
  
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2)) #sumstats sd
  
}else{ #if it is continuous 
  
  cat("treating GWAS trait as continuous.\n",sep='',file=file_log,append=TRUE)  
  
  #https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats-bbj/height.R
  #assumes sd(y) = 1 (beta std 1)
  #if not below sd(y) is reestimated 
  
  sd_ss = with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))
  
  sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
  
}

is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

cat(" ","\n", file=file_log,append=TRUE)

#make (temp) out dir
dir.create(paste0(out_path,opt$sumstats,'_',format(Sys.time(), '%d%B%Y')))

tmp <- paste0(out_path,opt$sumstats,'_',format(Sys.time(), '%d%B%Y'),'/')

#plot SD SS vs SD REF
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave(paste0(tmp,"sd_", opt$sumstats,".png"), width = 10, height = 7)


df_beta <- info_snp[!is_bad, ] #remove bad SNPs

cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)

cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)

# Warning if more than 50% variants have discordant SD
# see also https://github.com/privefl/bigsnpr/issues/281
if(sum(is_bad, na.rm = T) > (length(is_bad)*0.5)) {
  
  cat("WARNING: More than half the variants had a discordant SD. Imputing Neff for binary trait.\n
      Double check your input sumstats: reference population/per-variant sample size/effective N.\n\n",file=file_log,append=TRUE) # see: https://github.com/privefl/bigsnpr/issues/281.
  if(opt$type == T){
    #see eq 4 and 5 here  https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub#sec3.2
    
    info_snp$n_eff_imp <- (4 / sd_ldref^2 - info_snp$beta^2) / info_snp$beta_se^2
    
    #correcting difference in per allele effect size
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  #doi: 10.1093/bioinformatics/btw613
    
    info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
    
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    
  }else{
    
    cat("WARNING: More than half the variants had a discordant SD. Imputing Neff and estimating sd(y) for continuous trait.\n
      Double check your input sumstats: reference population/per-variant sample size/effective N.\n\n",file=file_log,append=TRUE) 
    
    info_snp$n_eff_imp <- (1 / (sd_ldref^2 - info_snp$beta^2)) / info_snp$beta_se^2
    
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  
    
    info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
    
    sd_y = with(info_snp, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01))) # https://github.com/privefl/bigsnpr/issues/349
    
    sd_ss = with(info_snp, sd_y / sqrt(n_eff * beta_se^2 + beta^2))
    
    sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
  }
  
  cat("N = ",sum(is.na(info_snp$n_eff_imp)), " SNPs have been filtered to correct for low sample size.\n"," ","\n",sep='',file=file_log,append=TRUE)
  
  cat("new iputed median N = ", median(info_snp$n_eff_imp, na.rm=T) ,"\n", file=file_log,append=TRUE) 
  
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  
  cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)
  
  
  #plot SD SS vs SD REF
  qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
    theme_bigstatsr() +
    coord_equal() +
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from allele frequencies of the LD reference",
         y = "Standard deviations derived from the summary statistics",
         color = "Removed?")
  
  ggsave(paste0(tmp,"sd_", opt$sumstats,".medianN.png"), width = 10, height = 7)
  
  df_beta <- info_snp[!is_bad, ] #remove bad SNPs
  
  cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)
  
}

#load test data
obj.bigsnp <- snp_attach(paste0(geno))

# shortcut for geno test set
G <- obj.bigsnp$genotypes

map_test <- dplyr::transmute(obj.bigsnp$map,
                             chr = as.integer(chromosome), pos = physical.pos,
                             a0 = allele2, a1 = allele1)

#remove duplicates (if any) in test data 
dups <- vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])
if (any(dups)) {
  cat("There are: ", sum( vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])), " duplicated physical positions in test data.\n" ," ","\n", sep='',file=file_log,append=TRUE)
  cat("removing duplicated...\n"," ","\n",file=file_log,append=TRUE) 
  map_test <-  map_test[!dups, ] }

# match with variants in test data and prepare files for matrix multiplication later
cat("Matching with test data...\n"," ","\n",file=file_log,append=TRUE)

map_pgs <- df_beta[1:4]; map_pgs$beta <- 1

map_pgs2 <- snp_match(map_pgs, map_test, return_flip_and_rev = T)

# match sumstats with reference data
cat(" N = ", dim(map_pgs2)[1], "SNPs have been matched with test data\n",
    "N =  ",sum(map_pgs2$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(map_pgs2$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

#in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_test[, c("chr", "pos")])
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_pgs2[, c("chr", "pos")])

#cat("in_test vector = TRUE is: ", sum(in_test), " long.\n"," ","\n",sep='',file=file_log,append=TRUE)

df_beta <- df_beta[in_test, ]

# The score-weight vectors below assume that df_beta and map_pgs2 are in
# exactly the same variant order. Stop rather than silently calculate
# misaligned scores if that assumption is violated.
if (nrow(df_beta) != nrow(map_pgs2)) {
  stop("Mismatch after target-data matching: df_beta has ", nrow(df_beta),
       " rows but map_pgs2 has ", nrow(map_pgs2), " rows.", call. = FALSE)
}

if (!all(df_beta$chr == map_pgs2$chr & df_beta$pos == map_pgs2$pos)) {
  stop("Variant-order mismatch between df_beta and map_pgs2 after target-data matching.",
       call. = FALSE)
}

cat("There are: ", dim(df_beta)[1], " SNPs in common with test data.\n"," ","\n",sep='',file=file_log,append=TRUE)

# LDSC
cat("Running LDSC...\n"," ","\n",file=file_log,append=TRUE)

(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))

h2_est <- ldsc[["h2"]]

cat(paste(names(ldsc),ldsc, sep ="=", collapse="; "),"\n"," ","\n", file=file_log, append=TRUE)

#ldsc ratio 
ldsc_int <- ldsc[["int"]]
mchi2 <- mean(chi2,na.rm=T)

ldscRATIO <- (ldsc_int - 1) /  (mchi2 - 1)

cat("Ratio = ",ldscRATIO,".\n"," ","\n",file=file_log,append=TRUE)  

cat("Running sparse matrix...\n"," ","\n",file=file_log,append=TRUE)

#sparse matrix
for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  
  corr_chr <- readRDS(paste0(misc_path,"LDref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, paste0(tmp,"corr_chr"), compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3  # file size in GB



cat("Running LDpred2-inf...\n"," ","\n",file=file_log,append=TRUE)

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)

beta_inf2 <- cbind(df_beta[,1:4], beta_inf) #bind with chr:pos:a0:a1

#save weights
write.table(beta_inf2, paste0(tmp,opt$sumstats,"_beta_inf.txt"), col.names=T,row.names=F,quote=F)

pred_inf <- big_prodVec(G, beta_inf * map_pgs2$beta,
                        ind.col = map_pgs2[["_NUM_ID_"]],
                        ncores = NCORES)


final_pred_inf <- cbind(obj.bigsnp$fam,pred_inf)		#bind with test map				

#save final scores
write.table(final_pred_inf, paste0(tmp,opt$sumstats,"_pred_inf.txt"), col.names=T, row.names=F, quote = F)


cat("Running LDpred2-auto...\n"," ","\n", file=file_log,append=TRUE)

# ------------------------------------------------------------------
# LDpred2-auto input checks
# ------------------------------------------------------------------

cat("Checking LDpred2-auto inputs...\n", file = file_log, append = TRUE)

if (nrow(df_beta) != nrow(corr) || nrow(df_beta) != ncol(corr)) {
  stop(
    "Dimension mismatch before LDpred2-auto: df_beta has ", nrow(df_beta),
    " rows, whereas corr is ", nrow(corr), " x ", ncol(corr), ".",
    call. = FALSE
  )
}

required_numeric <- c("beta", "beta_se", "n_eff", "ld")

for (v in required_numeric) {
  if (!v %in% names(df_beta)) {
    stop("Required column missing from df_beta: ", v, call. = FALSE)
  }
  
  n_bad <- sum(!is.finite(df_beta[[v]]))
  
  cat(v, ": non-finite values = ", n_bad, "\n",
      sep = "", file = file_log, append = TRUE)
  
  if (n_bad > 0) {
    stop("Non-finite values found in df_beta$", v, call. = FALSE)
  }
}

if (any(df_beta$beta_se <= 0)) {
  stop("df_beta$beta_se contains zero or negative values.", call. = FALSE)
}

if (any(df_beta$n_eff <= 0)) {
  stop("df_beta$n_eff contains zero or negative values.", call. = FALSE)
}

if (!is.finite(h2_est) || h2_est <= 0 || h2_est >= 1) {
  stop("Invalid LDpred2-auto starting heritability: ", h2_est,
       call. = FALSE)
}

if (!is.numeric(NCORES) || length(NCORES) != 1 || NCORES < 1) {
  stop("The number of cores must be a positive integer.", call. = FALSE)
}

cat("LDpred2-auto input checks passed.\n", file = file_log, append = TRUE)
cat(
  "LDpred2-auto settings: p_init = 1e-4 to 0.2; ",
  "shrink_corr = 0.95; use_MLE = FALSE; ",
  "allow_jump_sign = FALSE.\n",
  file = file_log, append = TRUE, sep = ""
)

# The official convergence guidance recommends a maximum initial p of
# 0.2 and trying use_MLE = FALSE when chains are unstable.
set.seed(1)

multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
  allow_jump_sign = FALSE,
  shrink_corr = 0.95,
  use_MLE = FALSE,
  ncores = NCORES
)

# Save all chains before convergence filtering so failed runs can be
# inspected without rerunning the expensive model.
saveRDS(
  multi_auto,
  file = paste0(tmp, opt$sumstats, "_multi_beta_auto.rds")
)

# ------------------------------------------------------------------
# Chain diagnostics
# ------------------------------------------------------------------

first_or_na <- function(x) {
  if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
}

chain_diagnostics <- do.call(
  rbind,
  lapply(seq_along(multi_auto), function(i) {
    auto <- multi_auto[[i]]
    
    finite_corr <- auto$corr_est[is.finite(auto$corr_est)]
    
    data.frame(
      chain = i,
      p_init = first_or_na(auto$p_init),
      p_est = first_or_na(auto$p_est),
      h2_est = first_or_na(auto$h2_est),
      path_p_finite = sum(is.finite(auto$path_p_est)),
      path_p_total = length(auto$path_p_est),
      path_h2_finite = sum(is.finite(auto$path_h2_est)),
      path_h2_total = length(auto$path_h2_est),
      beta_finite = sum(is.finite(auto$beta_est)),
      beta_total = length(auto$beta_est),
      corr_finite = length(finite_corr),
      corr_total = length(auto$corr_est),
      corr_range = if (length(finite_corr) >= 2) {
        diff(base::range(finite_corr))
      } else {
        NA_real_
      }
    )
  })
)

write.table(
  chain_diagnostics,
  file = paste0(tmp, opt$sumstats, "_auto_chain_diagnostics.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(chain_diagnostics)

cat(
  "LDpred2-auto chains with fully finite beta estimates: ",
  sum(chain_diagnostics$beta_finite == chain_diagnostics$beta_total),
  " of ", nrow(chain_diagnostics), "\n",
  "LDpred2-auto chains with finite final p and h2: ",
  sum(is.finite(chain_diagnostics$p_est) &
        is.finite(chain_diagnostics$h2_est)),
  " of ", nrow(chain_diagnostics), "\n",
  sep = "", file = file_log, append = TRUE
)

# Plot the chain with the largest number of finite trajectory values.
plot_chain <- which.max(
  chain_diagnostics$path_p_finite + chain_diagnostics$path_h2_finite
)
auto <- multi_auto[[plot_chain]]

p_trace <- data.frame(
  iteration = seq_along(auto$path_p_est),
  value = auto$path_p_est
)
p_trace <- p_trace[is.finite(p_trace$value) & p_trace$value > 0, , drop = FALSE]

h2_trace <- data.frame(
  iteration = seq_along(auto$path_h2_est),
  value = auto$path_h2_est
)
h2_trace <- h2_trace[is.finite(h2_trace$value), , drop = FALSE]

if (nrow(p_trace) > 0 && nrow(h2_trace) > 0) {
  p_plot <- ggplot(p_trace, aes(x = iteration, y = value)) +
    geom_line() +
    theme_bigstatsr() +
    scale_y_log10() +
    labs(y = "p", x = "Iteration")
  
  if (is.finite(auto$p_est)) {
    p_plot <- p_plot + geom_hline(yintercept = auto$p_est, color = "blue")
  }
  
  h2_plot <- ggplot(h2_trace, aes(x = iteration, y = value)) +
    geom_line() +
    theme_bigstatsr() +
    labs(y = "h2", x = "Iteration")
  
  if (is.finite(auto$h2_est)) {
    h2_plot <- h2_plot + geom_hline(yintercept = auto$h2_est, color = "blue")
  }
  
  chain_plot <- plot_grid(p_plot, h2_plot, ncol = 1, align = "hv")
  
  ggsave(
    paste0(tmp, "auto_chains_", opt$sumstats, ".png"),
    plot = chain_plot,
    width = 12,
    height = 10
  )
}

# A chain is eligible only when the final hyperparameters, all beta
# estimates and its correlation range are finite.
valid_chain <- with(
  chain_diagnostics,
  is.finite(p_est) &
    is.finite(h2_est) &
    beta_finite == beta_total &
    is.finite(corr_range)
)

if (!any(valid_chain)) {
  cat(
    "ERROR: All LDpred2-auto chains diverged. Diagnostics were saved to ",
    paste0(tmp, opt$sumstats, "_auto_chain_diagnostics.txt"), "\n",
    file = file_log, append = TRUE, sep = ""
  )
  
  # Remove the large temporary sparse matrix while retaining the saved
  # chain object and diagnostics for inspection.
  if (file.exists(paste0(tmp, "corr_chr.sbk"))) {
    file.remove(paste0(tmp, "corr_chr.sbk"))
  }
  file.copy(file_log, tmp, overwrite = TRUE)
  
  stop(
    "All LDpred2-auto chains diverged. Inspect the saved chain diagnostics, ",
    "summary-statistic QC and GWAS-LD reference matching.",
    call. = FALSE
  )
}

# Filter outlier chains and average the remaining finite beta estimates.
chain_range <- chain_diagnostics$corr_range

range_threshold <- 0.95 * quantile(
  chain_range[valid_chain],
  probs = 0.95,
  na.rm = TRUE
)

keep <- which(valid_chain & chain_range > range_threshold)

if (length(keep) == 0) {
  if (file.exists(paste0(tmp, "corr_chr.sbk"))) {
    file.remove(paste0(tmp, "corr_chr.sbk"))
  }
  file.copy(file_log, tmp, overwrite = TRUE)
  
  stop("No LDpred2-auto chains passed convergence filtering.", call. = FALSE)
}

cat(
  "Retaining LDpred2-auto chains: ", paste(keep, collapse = ", "), "\n",
  "Correlation-range threshold: ", range_threshold, "\n",
  sep = "", file = file_log, append = TRUE
)

beta_matrix <- vapply(
  multi_auto[keep],
  function(auto) auto$beta_est,
  numeric(nrow(df_beta))
)

final_beta_auto <- rowMeans(beta_matrix)

if (any(!is.finite(final_beta_auto))) {
  if (file.exists(paste0(tmp, "corr_chr.sbk"))) {
    file.remove(paste0(tmp, "corr_chr.sbk"))
  }
  file.copy(file_log, tmp, overwrite = TRUE)
  
  stop("The averaged LDpred2-auto effects contain non-finite values.",
       call. = FALSE)
}

#save final betas
final_beta_auto2 <- cbind(df_beta[,1:4], final_beta_auto) #bind with chr:pos:a0:a1
write.table(final_beta_auto2, paste0(tmp,opt$sumstats,"_final_beta_auto.txt"), col.names=T,row.names=F,quote=F)

#calculate score (matrix multiplication)
final_pred_auto <- big_prodVec(G, final_beta_auto * map_pgs2$beta, 
                               ind.col = map_pgs2[["_NUM_ID_"]],
                               ncores = NCORES)

final_pred_auto <- cbind(obj.bigsnp$fam,final_pred_auto)		#bind with test fam				

write.table(final_pred_auto, paste0(tmp,opt$sumstats,"_pred_auto.txt"), col.names=T, row.names=F, quote = F)

#remove sbk
file.remove(paste0(tmp, "corr_chr.sbk"))

end <- Sys.time()

cat(paste0("Analyses ended at ", end)," ","\n", file=file_log,append=TRUE)

cat(paste0("(Analyses took: ", round(as.numeric(difftime(end, start , units="mins")),digits=1) ," minutes)"),"\n"," ","\n", file=file_log,append=TRUE)

cat("###END###\n",file=file_log,append=TRUE)

#move logs and figures in out folder 
file.copy(file_log, tmp, overwrite = TRUE)
file.remove(file_log)