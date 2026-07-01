# Pathway PGS, MAGMA and LDpred2 pipeline

This repository contains a reproducible workflow for preparing GWAS summary statistics, running MAGMA gene and pathway analyses, and generating LDpred2 polygenic scores on the UCL Myriad HPC environment.

The current positive-control examples are height and BMI. The scripts are deliberately written as **one-trait-at-a-time** scripts to make debugging and re-running safer.

## Workflow summary

```text
0_format-sumstats-ldpred2.sh
    Raw GWAS summary statistics -> LDpred2-ready summary statistics

1_munge-sumstats-magma.R
    LDpred2-ready summary statistics + MAGMA reference .bim -> MAGMA-ready summary statistics

2_create-annotation-file.sh
    MAGMA reference .bim + NCBI gene locations -> reusable MAGMA .genes.annot file

3_MAGMA-run-gene-level.sh
    MAGMA-ready summary statistics + .genes.annot -> .genes.raw and .genes.out

4_run-MAGMA-pathways.sh
    .genes.raw + pathway .gmt files -> MAGMA .gsa.out pathway results

5_extract-sig-pathways.R
    .gsa.out files -> FDR/Bonferroni-adjusted tables and top pathways

A_ldpred2_auto_inf_qc_lift2_custom.R
    Custom LDpred2-inf/auto R pipeline with additional stability checks

B_subLDPred2.mcs.uclhg.sh
    Myriad job-submission script for A_ldpred2_auto_inf_qc_lift2_custom.R
```

## Intended Myriad paths

The MAGMA scripts currently assume:

```bash
/myriadfs/home/ucju659/SOFTWARE/MAGMA
```

with this structure:

```text
SOFTWARE/MAGMA/
├── v1.10/
│   └── magma
├── g1000_eur/
│   ├── g1000_eur.bed
│   ├── g1000_eur.bim
│   └── g1000_eur.fam
├── gene_locations/
│   ├── NCBI37.3.gene.loc
│   └── GRCh37_35kb_10kb.genes.annot
├── pathways/
│   ├── c5.go.v2026.1.Hs.entrez.gmt
│   └── c2.cp.reactome.v2026.1.Hs.entrez.gmt
├── munged/
├── results/
└── logs/
```

The LDpred2 scripts currently assume:

```bash
/myriadfs/home/ucju659/uclhg-mcs-pgs
/myriadfs/home/ucju659/SUMSTATS/ldpred2_ready
/myriadfs/home/ucju659/misc/hapmap3plus
/myriadfs/home/ucju659/REFERENCE/UCLhg/MCS
```

Edit paths at the top of each script before running on a different system.

## Step 0: Format raw GWAS summary statistics for LDpred2

Run interactively on Myriad:

```bash
qrsh -pe smp 1 -l mem=16G,h_rt=0:30:00 -now no
bash scripts/0_format-sumstats-ldpred2.sh
```

Edit only the user settings at the top:

```bash
IN=...
OUTDIR=...
OUT=...
N_GWAS=...
MIN_MAF=...
MIN_INFO=...
```

Output columns:

```text
CHR BP A2 A1 N BETA SE MAF P INFO
```

`A1` is the effect allele and `A2` is the other allele. Ambiguous A/T and C/G SNPs are removed.

## Step 1: Munge summary statistics for MAGMA

Run locally or on Myriad:

```bash
Rscript scripts/1_munge-sumstats-magma.R
```

This script maps variants by `CHR:BP` to the MAGMA reference `.bim` file and uses the `.bim` SNP ID as the MAGMA `SNP` column.

Required output columns for MAGMA:

```text
SNP P N
```

Extra audit columns are retained.

## Step 2: Create MAGMA gene annotation

Run once per reference panel, genome build, gene-location file and annotation window:

```bash
bash scripts/2_create-annotation-file.sh
```

Output:

```text
gene_locations/GRCh37_35kb_10kb.genes.annot
```

The current annotation window is 35 kb upstream and 10 kb downstream.

## Step 3: Run MAGMA gene-level analysis

Submit one trait at a time on Myriad:

```bash
qsub scripts/3_MAGMA-run-gene-level.sh
```

Edit only:

```bash
MUNGED_FILE=...
GENE_OUT=...
```

Outputs:

```text
<GENE_OUT>.genes.raw
<GENE_OUT>.genes.out
<GENE_OUT>.log
```

The `.genes.raw` file is the input for pathway analysis.

## Step 4: Run MAGMA pathway analysis

Run after `.genes.raw` exists:

```bash
qsub scripts/4_run-MAGMA-pathways.sh
```

This submits GO and Reactome as separate MAGMA gene-set analyses within one Myriad job. It does **not** rerun the gene-level step.

Outputs:

```text
<trait>_GO.gsa.out
<trait>_Reactome.gsa.out
```

## Step 5: Extract significant pathways

Run:

```bash
Rscript scripts/5_extract-sig-pathways.R
```

This adds Benjamini-Hochberg FDR and Bonferroni correction separately within each pathway collection.

Outputs for each collection:

```text
<prefix>_all_with_FDR.tsv
<prefix>_nominal_P05.tsv
<prefix>_FDR05.tsv
<prefix>_Bonferroni05.tsv
<prefix>_top50.tsv
```

## Step A/B: Run LDpred2 PGS

Prepare `sumstats_list.csv` in the LDpred2 working directory. By default, the submission script expects no header:

```text
GIANT_UKBB_BMI_2018_ALL_SITES.ldpred2.gz,FALSE
GIANT_HEIGHT_YENGO_2022_EUR.ldpred2.gz,FALSE
```

Submit:

```bash
qsub scripts/B_subLDPred2.mcs.uclhg.sh
```

If there are two rows in `sumstats_list.csv`, change the SGE array line in `B_subLDPred2.mcs.uclhg.sh` to:

```bash
#$ -t 1-2
```

## Reproducibility checklist

For each trait, record:

- GWAS source and citation
- download date
- ancestry
- genome build
- sample size handling
- effect allele definition
- filters applied
- summary-statistic columns used
- reference panel used
- MAGMA version
- LDpred2/bigsnpr version
- pathway database, collection and version
- gene identifier type
- SNP-to-gene window
- Myriad job ID
- output file names

For pathway analyses, record pathway provenance explicitly, for example:

```text
MSigDB C5 Gene Ontology v2026.1, Entrez IDs, GRCh37, 35 kb upstream / 10 kb downstream, 1000 Genomes EUR reference.
MSigDB C2 Reactome v2026.1, Entrez IDs, GRCh37, 35 kb upstream / 10 kb downstream, 1000 Genomes EUR reference.
```

## Interpretation cautions

MAGMA gene-set results are competitive pathway-level statistical tests based on SNP-to-gene mapping assumptions. They do not prove that a pathway, gene or biological mechanism is causal.

Pathway-based PGS/PRSet-style analyses produce individual-level genetic liability scores restricted to defined gene sets or pathways. Their interpretation depends on pathway definition, LD, ancestry, SNP-to-gene mapping, sample size, p-value thresholding and covariate structure.

Positive-control traits such as height and BMI are useful for checking that the pipeline behaves sensibly, but results still depend on GWAS power, ancestry matching, phenotype measurement and sample overlap.

## No-exit Bash behaviour

The Bash scripts in this repository are currently written in **no-exit mode**. Failed checks print `WARNING:` messages and skip unsafe analysis commands rather than terminating the shell. This is intentional for interactive Myriad work, especially if a script is accidentally sourced rather than executed.

Important trade-off: this is safer for an interactive shell, but less strict for submitted jobs. A Myriad job may finish without scheduler failure even if the analysis was skipped. Always inspect the scheduler log and confirm that the expected output files were created and are non-empty.

Recommended usage remains:

```bash
bash scripts/2_create-annotation-file.sh
qsub scripts/3_MAGMA-run-gene-level.sh
qsub scripts/4_run-MAGMA-pathways.sh
```

Avoid sourcing these scripts with `source script.sh` unless you are deliberately running them line by line.

## Interactive Myriad note

The Bash scripts deliberately avoid custom Bash helper functions such as `warn()`. Warnings are printed using plain `echo "WARNING: ..." >&2` statements. This makes the scripts easier to run line-by-line in an interactive Myriad session without needing to define helper functions first.
