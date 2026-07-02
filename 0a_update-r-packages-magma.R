qrsh -l h_rt=4:00:00,mem=8G,tmpfs=10G -pe smp 2 -now no

module -f unload compilers mpi gcc-libs
module load r/4.5.1-openblas/gnu-10.2.0
unset R_LIBS
export R_LIBS_USER="/myriadfs/home/ucju659/MyRLibs/R-4.5.1"

mkdir -p "$R_LIBS_USER"

R -- vanilla

.libPaths(c("/myriadfs/home/ucju659/MyRLibs/R-4.5.1", .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = 2
)

# For basic MAGMA

install.packages(
  c(
    "BiocManager",
    "data.table",
    "optparse",
    "R.utils",
    "bigsnpr",
    "bigstatsr",
    "msigdbr",
    "dplyr",
    "readr",
    "ggplot2"
  ),
  dependencies = c("Depends", "Imports", "LinkingTo")
)

BiocManager::install(
  version = "3.22",
  ask = FALSE,
  update = FALSE
)

BiocManager::install(
  c(
    "GenomicRanges",
    "GenomeInfoDb",
    "rtracklayer",
    "MungeSumstats",
    "SNPlocs.Hsapiens.dbSNP144.GRCh37",
    "BSgenome.Hsapiens.1000genomes.hs37d5"
  ),
  ask = FALSE,
  update = FALSE
)

# For MAGMA celltyping
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", dependencies = c("Depends", "Imports", "LinkingTo"))
}

remotes::install_github(
  "neurogenomics/MAGMA_Celltyping",
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never",
  build_vignettes = FALSE
)

# Check installation
pkgs <- c(
  "data.table",
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "bigsnpr",
  "bigstatsr",
  "MungeSumstats",
  "SNPlocs.Hsapiens.dbSNP144.GRCh37",
  "BSgenome.Hsapiens.1000genomes.hs37d5"
)

print(sapply(pkgs, requireNamespace, quietly = TRUE))
print(.libPaths())
sessionInfo()

# Try to install missing packages
BiocManager::install(
  c(
    "SNPlocs.Hsapiens.dbSNP144.GRCh37",
    "BSgenome.Hsapiens.1000genomes.hs37d5"
  ),
  ask = FALSE,
  update = FALSE
)