# Install every R package the analysis uses into the project's renv library, then record them in
# renv.lock. Run once, after tools/bootstrap.sh; everything compiles from source (several hours):
#
#   sbatch scripts/slurm/setup/install-r-packages.sh
#
# Afterwards the template's usual flow applies (renv::install(), use the package, renv::snapshot()).
# Packages that fail are listed at the end: fix the cause (often a C library missing from
# environment.yml), then run this again; packages already installed are skipped.

ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(Ncpus = ncpu)
Sys.setenv(MAKEFLAGS = paste0("-j", ncpu))
# Packages that still ask for C++11 (presto) get C++14, which current RcppArmadillo needs
Sys.setenv(R_MAKEVARS_USER = normalizePath("scripts/R/Makevars"))
# R 4.5 goes with Bioconductor 3.22
renv::settings$bioconductor.version("3.22")

bioc <- c(
  "AnnotationDbi", "anndataR", "BiocGenerics", "BiocParallel", "biomaRt",
  "BSgenome.Mmusculus.UCSC.mm10", "celldex", "ComplexHeatmap", "decoupleR", "destiny",
  "dittoSeq", "DropletUtils", "edgeR", "EnsDb.Mmusculus.v79", "ensembldb", "fgsea",
  "GenomeInfoDb", "GenomicRanges", "glmGamPoi", "IRanges", "JASPAR2024", "limma", "miloR",
  "motifmatchr", "multtest", "Nebulosa", "org.Hs.eg.db", "org.Mm.eg.db", "orthogene",
  "plyranges", "S4Vectors", "scDblFinder", "SingleCellExperiment", "SingleR", "slingshot",
  "SummarizedExperiment", "TFBSTools", "UCell", "zellkonverter",
  # WGCNA / hdWGCNA dependencies
  "impute", "preprocessCore", "GeneOverlap"
)

cran <- c(
  "ade4", "ape", "base64enc", "BimodalIndex", "circlize", "colorspace", "correlation", "cowplot",
  "data.table", "dendsort", "devtools", "diptest", "dplyr", "DT", "e1071", "energy",
  "fastcluster", "future", "GGally", "ggalluvial", "ggnewscale", "ggplot2", "ggpointdensity",
  "ggraph", "ggrepel", "ggridges", "ggwordcloud", "glmpca", "glue", "graphlayouts", "harmony",
  "htmltools", "httr", "igraph", "infotheo", "knitr", "magrittr", "MASS", "Matrix",
  "matrixStats", "mclust", "metap", "msigdbr", "paletteer", "patchwork", "pbapply", "pheatmap",
  "phateR", "plotly", "Polychrome", "presto", "psych", "purrr", "qs2", "RANN", "RColorBrewer",
  "RcppML", "readr", "remotes", "reshape2", "reticulate", "Rfast2", "rlang", "RSQLite", "scales",
  "scatterpie", "sccore", "scCustomize", "seriation", "Seurat", "SeuratObject", "Signac", "spdep",
  "stringr", "SuperCell", "tibble", "tidygraph", "tidyr", "tidyverse", "viridis", "WGCNA", "xfun",
  "yaml", "hrbrthemes",
  # SuperCell and SCENT dependencies that aren't resolved on their own
  "WeightedCluster", "corpcor", "Hmisc", "qlcMatrix"
)

# Not on CRAN or Bioconductor 3.22 (any more): name = the package, value = where it comes from
other <- c(
  smoother     = "smoother@1.1",                        # CRAN archive (destiny needs it)
  symphony     = "immunogenomics/symphony",
  Rmagic       = "KrishnaswamyLab/MAGIC/Rmagic",
  deMULTIplex2 = "Gartner-Lab/deMULTIplex2",
  DElegate     = "cancerbits/DElegate",
  SCENT        = "aet21/SCENT",
  CytoTRACE2   = "digitalcytometry/cytotrace2/cytotrace2_r",
  hdWGCNA      = "smorabit/hdWGCNA@dev",
  Azimuth      = "satijalab/azimuth",
  CellChat     = "jinworks/CellChat",
  nichenetr    = "saeyslab/nichenetr",
  SCENIC       = "aertslab/SCENIC",
  miloDE       = "MarioniLab/miloDE",
  ggforce      = "thomasp85/ggforce"                    # dev version: the CRAN release had a broken
                                                        # check_labeller reference against ggplot2 (2026-04)
)

installed <- function(pkg) nzchar(system.file(package = pkg))
failed <- character()
install <- function(spec, pkg = spec) {
  if (installed(pkg)) return(invisible())
  message("\n>>> ", spec)
  tryCatch(renv::install(spec, prompt = FALSE),
           error = function(e) message("FAILED: ", spec, ": ", conditionMessage(e)))
  if (!installed(pkg)) failed <<- c(failed, spec)
}
for (p in bioc) install(paste0("bioc::", p), p)
for (p in cran) install(p)
for (p in names(other)) install(other[[p]], p)

# Record what the code uses (and _dependencies.R lists) in renv.lock
renv::snapshot(prompt = FALSE, force = length(failed) > 0)

if (length(failed)) {
  message("\nNot installed (", length(failed), "): ", paste(failed, collapse = ", "))
  quit(status = 1)
}
message("\nAll packages installed; renv.lock updated.")
