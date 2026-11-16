# Dev and infra packages renv should track (never sourced)
library(IRkernel)
library(languageserver)
library(vscDebugger)
library(jsonlite)
library(rlang)
library(R6)
library(ggplot2)
library(yaml)
library(rmarkdown)
library(knitr)

# Pipeline and reports (used by _targets.R / tar_quarto(), listed so renv records them)
library(targets)
library(tarchetypes)
library(quarto)
library(visNetwork)

# Used behind the scenes by packages the analysis calls (no code loads them by name)
library(qs2)        # targets' format = "qs"
library(presto)     # fast Wilcoxon tests for Seurat::FindAllMarkers
library(glmGamPoi)  # Seurat::SCTransform (v2)
library(Rfast2)     # Seurat::RunMoransI
library(metap)      # Seurat::FindConservedMarkers
library(multtest)   # metap
