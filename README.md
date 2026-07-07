# MAGMA pathway discovery and pathway-based polygenic scoring workflow

This repository contains a reproducible workflow for moving from GWAS summary statistics to:

1. MAGMA gene-level results.
2. MAGMA gene-set/pathway results.
3. Selected pathway sets for follow-up scoring.
4. PRSice/PRSet clumping-and-thresholding pathway scores.
5. LDpred2 posterior-beta pathway scores.
6. Genome-wide background scores excluding SNPs mapped to the selected pathway-gene set.

The workflow was developed for use on the UCL Myriad high-performance computing environment, but the scripts are written as editable templates and can be adapted to other Unix-like systems.

## Repository scope

This repository is for pathway discovery, pathway/background score construction, and optional post-GWAS tissue/cell-type annotation. It does **not** perform MCS phenotype construction, child-mother-father score merging, role-specific PC residualisation, or downstream trio model estimation.

For the downstream MCS pathway-score merge and trio-analysis workflow, see:

```text
https://github.com/joseMorosoli/pathway-based-trios.git
```

The main scientific use case is pathway-based polygenic scoring. Tissue/cell-type scripts are included as optional exploratory utilities and are not required for the core pathway PRSet/LDpred2 pipeline.

## Current repository layout

Scripts currently live in the repository root, not in a separate `scripts/` directory. Commands in this README therefore call scripts as `./script_name.R` or `bash script_name.sh` from the repository root.

```text
magma-to-pathway/
├── archive/                         # Old scripts, exploratory attempts, or inactive materials
├── data/                            # Local metadata/manifests or small non-sensitive helper inputs
├── logs/                            # Local job logs; do not commit substantive logs
├── munged/                          # Locally generated formatted summary statistics; do not commit
├── results/                         # Locally generated MAGMA/PRSice/LDpred2 outputs; do not commit
├── tissue-cell-analyses/            # Optional tissue/cell-type analyses and bridge scripts
├── LICENSE
├── .gitignore
├── README.md
├── magma-to-pathway.Rproj
├── 0a_format-sumstats-ldpred2.sh
├── 0b_format-sumstats-prsice-rsid.sh
├── 0c_update-r-packages-magma.R
├── 0d_record-software-versions.R
├── 1_munge-sumstats-magma.R
├── 2_create-annotation-file.sh
├── 3_MAGMA-run-gene-level.sh
├── 4a_make-short-id-gmts.R
├── 4b_run-MAGMA-pathways.sh
├── 5a_extract-sig-pathways.R
├── 5b_check-selected-pathways.R
├── 6a_make-selected-pathway-gmt.R
├── 6b_create-nonpathway-snp-list.sh
├── 6c_score-nonpathway-PRSice.sh
├── 6d_score-nonpathway-LDpred2.sh
├── 7a_score-selected-pathways-PRSice.sh
├── 7b_score-selected-pathways-LDpred2.sh
├── A_ldpred2_auto_inf_qc_lift2_custom.R
└── B_subLDPred2.mcs.uclhg.sh
```

The following local R session files may appear during development but should not normally be committed:

```text
.RData
.Rhistory
.RDataTmp
```

## Data access and privacy

This repository should contain code, documentation, and lightweight non-sensitive metadata only. It should not contain raw GWAS summary statistics, MCS phenotype files, genotype files, target PLINK/bigsnpr files, generated polygenic scores, MAGMA outputs, or individual-level pathway scores.

Generated score files are for private/internal use only and must be stored on approved secure servers. Users are responsible for ensuring that their use of MCS data, genotype data, GWAS summary statistics, and derived scores complies with the relevant data-access agreements.

Do not commit any individual-level data, genotype files, summary statistics, score files, logs containing identifiable or sensitive paths, or large intermediate outputs to GitHub.

Recommended use of local folders:

```text
data/      small non-sensitive manifests, toy examples, or documentation inputs only
munged/    formatted GWAS files created by 0a/0b/1; normally gitignored
results/   MAGMA, PRSice, PRSet, LDpred2 outputs; normally gitignored
logs/      Myriad job logs; normally gitignored
archive/   old or superseded scripts kept for audit, not part of the active pipeline
```

## Main software

The workflow was developed with:

- MAGMA v1.10.
- PRSice-2 / PRSet v2.3.3.
- R with packages including `data.table`, `bigsnpr`, `bigstatsr`, `GenomicRanges`, `GenomeInfoDb`, `rtracklayer`, `MungeSumstats`, and supporting tidyverse/plotting packages.
- MSigDB Human v2026.1.Hs pathway files.
- LDpred2/bigsnpr posterior-beta scoring.

Exact R package versions should be recorded at run time with:

```bash
Rscript --vanilla 0d_record-software-versions.R
```

If a command-line software version script is added later, edit all local paths before running it. Do not assume paths from the development environment are portable.

## External resources and provenance

### Pathway databases

The pathway inputs use MSigDB Human v2026.1.Hs files:

```text
c2.cp.reactome.v2026.1.Hs.entrez.gmt
c2.cp.reactome.v2026.1.Hs.symbols.gmt
c5.go.bp.v2026.1.Hs.symbols.gmt
c5.go.cc.v2026.1.Hs.symbols.gmt
c5.go.mf.v2026.1.Hs.symbols.gmt
c5.go.v2026.1.Hs.entrez.gmt
c5.go.v2026.1.Hs.symbols.gmt
```

MAGMA gene-set analysis uses Entrez-ID GMTs. PRSice/PRSet and LDpred2 pathway scoring use gene-symbol GMTs because the SNP-to-gene mapping step uses a GRCh37 GTF with gene names.

### Genome build and SNP-to-gene mapping

The current project convention is:

```text
Genome build: GRCh37
SNP-to-gene window: 35 kb upstream, 10 kb downstream
```

The MAGMA annotation step creates a reusable `.genes.annot` file. The downstream pathway and background scoring scripts use the same 35 kb / 10 kb convention when mapping SNPs to genes from the GTF.

The project has used a GRCh37 GTF whose filename indicates Ensembl GRCh37 release 87 (`Homo_sapiens.GRCh37.87.gtf`). Optional PRSet scoring of tissue/cell-type-derived GMTs can also use GENCODE Release 19 GRCh37 if the gene identifiers are compatible. Verify the exact GTF source, release, genome build, and gene identifier field before reporting final analyses.

### MAGMA reference data

MAGMA software and reference files should be obtained from the CNCR MAGMA website:

```text
https://cncr.nl/research/magma/
```

The MAGMA website provides 1000 Genomes Phase 3 reference data files with SNP locations in human genome Build 37.

### LDpred2 LD reference

The LDpred2 HapMap3+ LD reference used by this project is available from Figshare:

```text
https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061
```

### Optional tissue/cell-type resources

The optional tissue/cell-type workflow may use:

- FUMA GTEx v8 general gene expression file.
- FUMA GTEx v8 specific gene expression file.
- MAGMA.Celltyping / `ewceData` cell-type datasets, for example `ctd_allKI`, `ctd_DRONC_human`, `ctd_AIBS`, and `ctd_Tasic`.
- A GRCh37-compatible GTF for optional PRSet scoring of top-decile tissue/cell-type gene sets.

These resources must be recorded with: database/source, version, genome build where applicable, download date, gene identifier type, filtering criteria, and any threshold used to define a binary gene set.

### GWAS examples used in this project

The scripts are written one trait at a time. Worked traits include:

- BMI: GIANT + UK Biobank 2018 BMI summary statistics from Zenodo: `https://zenodo.org/records/1251813`
- Height: GIANT/Yengo 2022 height GWAS summary statistics and PGS weights, including `Yengo.2022.height.GWAS.EUR.gz`.
- F4 psychiatric/internalising factor: Grotzinger et al., *Nature*, DOI `10.1038/s41586-025-09820-3`.

Users should check the README or metadata accompanying each GWAS before running the scripts, especially for sample size, ancestry, genome build, allele definitions, and whether per-SNP sample size is available.

## Important terminology

Throughout this repository, “non-pathway” means **outside the selected pathway-gene set**, not outside all known biological pathways.

Operationally, the non-pathway/background score is a genome-wide background score excluding SNPs mapped to genes in the selected pathway set. It is the complement of the selected pathway SNP union in the target genotype data, using the same genome build and SNP-to-gene window.

Do not interpret the non-pathway score as a generic “all non-biological-pathway SNPs” score.

For optional tissue/cell-type work, distinguish:

```text
MAGMA tissue/cell-type gene-property analysis:
  A post-GWAS enrichment/annotation analysis using continuous tissue or cell-type expression/specificity values.

Tissue/cell-type-derived PRSet scores:
  Optional exploratory PRSet scores created only after converting continuous tissue/cell-type values into binary gene sets, for example top 10% expressed/specific genes.
```

A top-decile tissue/cell-type PRSet score is not a continuous tissue-weighted score. It is a PRS restricted to SNPs mapped to genes passing a chosen tissue/cell-type threshold.

## Pipeline overview

Run scripts one trait at a time by editing the user-settings block at the top of each script. The order below assumes all external software, reference files, GWAS summary statistics, target genotype files, and pathway files have already been obtained and stored securely.

### 0. Setup and summary-statistic formatting

#### `0c_update-r-packages-magma.R`

Optional helper to install or update R packages needed for the local MAGMA/pathway workflow. Use carefully on Myriad. For stable analyses, record package versions and avoid changing the environment mid-project.

#### `0d_record-software-versions.R`

Record R version and R package versions for reproducibility.

Run from the repository root:

```bash
Rscript --vanilla 0d_record-software-versions.R
```

#### `0a_format-sumstats-ldpred2.sh`

Convert raw GWAS summary statistics into a simple LDpred2-ready format.

Output columns:

```text
CHR BP A2 A1 N BETA SE MAF P INFO
```

Conventions:

- `A1` is the effect allele.
- `A2` is the other/non-effect allele.
- Ambiguous A/T and C/G SNPs are removed.
- Duplicate physical positions are removed, keeping the first occurrence.
- If the raw file has no per-SNP `N`, the user must provide a defensible fallback `N_GWAS` value.

#### `0b_format-sumstats-prsice-rsid.sh`

Convert LDpred2-ready summary statistics to PRSice-ready format by mapping `CHR:BP` to the rsID used in the target `.bim` file.

Output columns:

```text
MarkerName CHR POS A1 A2 N BETA SE MAF P
```

`MarkerName` should match SNP IDs in the target PLINK files.

### 1. MAGMA gene-level analysis

#### `1_munge-sumstats-magma.R`

Munge summary statistics into MAGMA-ready format.

For standard MAGMA gene-level analyses with per-SNP sample size, the expected minimal columns are:

```text
SNP P N
```

For GWAS where SNP-level sample size is not available, do not fabricate per-SNP `N`. Use a documented fixed effective sample size only when this is defensible and report it in analysis notes.

#### `2_create-annotation-file.sh`

Create the reusable MAGMA SNP-to-gene annotation file.

Current convention:

```text
Genome build: GRCh37
Window: 35 kb upstream, 10 kb downstream
```

Rerun this step only if the reference `.bim`, gene-location file, genome build, or window changes.

#### `3_MAGMA-run-gene-level.sh`

Run MAGMA gene-level analysis for one trait at a time.

Expected outputs:

```text
<GENE_OUT>.genes.raw
<GENE_OUT>.genes.out
<GENE_OUT>.log
```

### 2. MAGMA pathway analysis and pathway selection

#### `4a_make-short-id-gmts.R`

Create short-ID GMT files to avoid long pathway names being truncated in MAGMA output, while preserving a mapping back to full pathway names.

#### `4b_run-MAGMA-pathways.sh`

Run MAGMA gene-set/pathway analyses using GO and Reactome as separate gene-set collections.

Expected outputs:

```text
<TRAIT>_GO.gsa.out
<TRAIT>_Reactome.gsa.out
```

Multiple-testing correction and selected-pathway extraction are handled by later scripts.

#### `5a_extract-sig-pathways.R`

Extract significant pathway results using the selected correction threshold. The project convention has used options such as:

```text
Bonferroni05
FDR05
nominal_P05
top50
```

#### `5b_check-selected-pathways.R`

Check selected pathways before score construction. This step is intended to catch missing pathway IDs, naming problems, or mismatches between MAGMA output and source GMT files.

#### `6a_make-selected-pathway-gmt.R`

Create selected pathway GMT files for downstream PRSice/PRSet and LDpred2 scoring.

Typical selected GMT output:

```text
<TRAIT>_selected_<SELECTION>_GO.symbols.gmt
<TRAIT>_selected_<SELECTION>_Reactome.symbols.gmt
<TRAIT>_selected_<SELECTION>_GO_Reactome.symbols.gmt
```

The combined GO+Reactome selected GMT is used for main selected-pathway scoring and background-score construction.

### 3. Background/non-selected-pathway scores

#### `6b_create-nonpathway-snp-list.sh`

Create the selected pathway SNP union and its complement in the target genotype data.

Expected outputs:

```text
pathway_union_snps.txt
nonpathway_snps.extract.txt
nonpathway_snp_definition_metadata.tsv
selected_pathway_genes_from_gmt.txt
selected_pathway_genes_found_in_gtf.txt
```

The metadata file records the selected GMT, target genotype file, GTF, window size, number of selected pathways, number of genes, number of pathway-union SNPs, and number of non-pathway SNPs.

#### `6c_score-nonpathway-PRSice.sh`

Create a PRSice C+T background score excluding selected-pathway SNPs.

The script uses:

```text
--extract nonpathway_snps.extract.txt
--bar-levels 1
--fastscore
--no-regress
--all-score
--print-snp
```

The script also creates or copies a p1-only `.all_score` file for downstream use:

```text
<TRAIT>_selected_<SELECTION>_nonpathway_PRSice.p1_only.all_score
```

#### `6d_score-nonpathway-LDpred2.sh`

Create an LDpred2 background score by restricting LDpred2 posterior betas to the non-selected-pathway SNP complement.

Expected outputs:

```text
<TRAIT>_selected_<SELECTION>_nonpathway_LDpred2.tsv.gz
<TRAIT>_selected_<SELECTION>_nonpathway_LDpred2.metadata.tsv
<TRAIT>_selected_<SELECTION>_nonpathway_LDpred2.snps_used.tsv.gz
<TRAIT>_selected_<SELECTION>_nonpathway_LDpred2.rds
```

### 4. Selected pathway scores

#### `7a_score-selected-pathways-PRSice.sh`

Create PRSice/PRSet clumping-and-thresholding scores for each selected pathway.

The script uses:

```text
--gtf <GRCh37 GTF>
--msigdb <selected pathway GMT>
--wind-5 35kb
--wind-3 10kb
--clump-kb 250
--clump-r2 0.1
--clump-p 1
--bar-levels 1
--no-regress
--all-score
--print-snp
```

The intended downstream file is the p1-only slim file:

```text
<TRAIT>_selected_<SELECTION>_PRSet_CT_GRCh37_35kb_10kb.p1_only.all_score
```

#### `7b_score-selected-pathways-LDpred2.sh`

Create one LDpred2 posterior-beta score per selected pathway.

Expected outputs:

```text
<TRAIT>_selected_<SELECTION>_LDpred2_pathway_scores.tsv.gz
<TRAIT>_selected_<SELECTION>_LDpred2_pathway_scores.metadata.tsv
<TRAIT>_selected_<SELECTION>_LDpred2_pathway_scores.snps_used.tsv.gz
<TRAIT>_selected_<SELECTION>_LDpred2_pathway_scores.rds
```

The metadata file records pathway-level counts of genes, genes found in the GTF, target SNPs mapped to the pathway, and LDpred2 SNPs used.

## Auxiliary LDpred2 scripts and credit

`A_ldpred2_auto_inf_qc_lift2_custom.R` and `B_subLDPred2.mcs.uclhg.sh` are auxiliary LDpred2 scripts adapted from Andrea G. Allegrini’s LDpred2 pipeline:

```text
https://github.com/AndreAllegrini/LDpred2
```

They have been modified for this project to improve compatibility with recent LDpred2-auto behaviour and to reduce errors in automated posterior-beta selection. Users should validate these scripts when changing GWAS inputs, target genotype data, LD references, genome builds, or software versions.

## Minimal command order

From the repository root, a typical pathway workflow is:

```bash
# 0. Optional: record software versions
Rscript --vanilla 0d_record-software-versions.R

# 1. Format summary statistics
bash 0a_format-sumstats-ldpred2.sh
bash 0b_format-sumstats-prsice-rsid.sh

# 2. MAGMA gene-level and pathway analyses
Rscript --vanilla 1_munge-sumstats-magma.R
bash 2_create-annotation-file.sh
bash 3_MAGMA-run-gene-level.sh
Rscript --vanilla 4a_make-short-id-gmts.R
bash 4b_run-MAGMA-pathways.sh

# 3. Select pathways and make selected GMTs
Rscript --vanilla 5a_extract-sig-pathways.R
Rscript --vanilla 5b_check-selected-pathways.R
Rscript --vanilla 6a_make-selected-pathway-gmt.R

# 4. Create background and selected-pathway scores
bash 6b_create-nonpathway-snp-list.sh
bash 6c_score-nonpathway-PRSice.sh
bash 6d_score-nonpathway-LDpred2.sh
bash 7a_score-selected-pathways-PRSice.sh
bash 7b_score-selected-pathways-LDpred2.sh
```

## Quality-control checklist

Before using any generated score file, check:

- The selected pathway counts match expectations.
- Selected pathway names match the source GMT files.
- The GTF gene identifier type matches the selected GMT gene identifier type.
- `nonpathway_snp_definition_metadata.tsv` contains non-zero selected genes, pathway-union SNPs, and non-pathway SNPs.
- PRSice `.all_score` files have the expected number of columns.
- PRSice selected-pathway files use the p1-only slim files for downstream merging.
- LDpred2 score files have non-zero variance and plausible SNP counts in metadata.
- Height scores predict measured height in the target cohort as a positive-control check.
- BMI scores predict measured BMI in the target cohort as a positive-control check.

For PRSice p1-only selected-pathway files, expected columns are:

```text
2 ID columns + number of selected pathways
```

For PRSice p1-only background files, expected columns are:

```text
FID IID <one score column>
```
## Suggested Myriad setup

Example R 4.5.1 setup:

```bash
module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0
unset R_LIBS
export R_LIBS_USER="/myriadfs/home/<USER>/MyRLibs/R-4.5.1"
```

Example threading controls:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```
## References and resources

- MAGMA software and documentation: https://cncr.nl/research/magma/
- de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set analysis of GWAS data. *PLoS Computational Biology*. 2015;11(4):e1004219. DOI: 10.1371/journal.pcbi.1004219.
- PRSice-2 documentation: https://choishingwan.github.io/PRSice/
- Choi SW, O'Reilly PF. PRSice-2: Polygenic Risk Score software for biobank-scale data. *GigaScience*. 2019. DOI: 10.1093/gigascience/giz082.
- Choi SW, García-González J, et al. Pathway-based polygenic risk score analyses and software. *PLOS Genetics*. 2023. DOI: 10.1371/journal.pgen.1010624.
- LDpred2/bigsnpr documentation: https://privefl.github.io/bigsnpr/articles/LDpred2.html
- Privé F, Arbel J, Vilhjálmsson BJ. LDpred2: better, faster, stronger. *Bioinformatics*. 2021. DOI: 10.1093/bioinformatics/btaa1029.
- LDpred2 HapMap3+ LD reference: https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061
- MSigDB collections: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- MAGMA.Celltyping documentation: https://neurogenomics.github.io/MAGMA_Celltyping/
- FUMA Downloads and tutorials: https://fuma.ctglab.nl/
- GIANT BMI 2018 Zenodo record: https://zenodo.org/records/1251813
- GIANT/Yengo 2022 height data: https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files
- Yengo L, et al. A saturated map of common genetic variants associated with human height. *Nature*. 2022. DOI: 10.1038/s41586-022-05275-y.
- Grotzinger AD, et al. Mapping the genetic landscape across 14 psychiatric disorders. *Nature*. DOI: 10.1038/s41586-025-09820-3.
- 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature*. 2015. DOI: 10.1038/nature15393.

## Licence

This repository is released under the MIT Licence. See `LICENSE`.
