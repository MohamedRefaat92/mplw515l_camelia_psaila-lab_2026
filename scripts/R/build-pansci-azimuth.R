#!/usr/bin/env Rscript
# RETIRED: the pansci_stromal / pansci_immune targets this reads are commented out in _targets.R
# (and so are their projections); kept for when the PanSci references come back.
#
# Build Azimuth references from PanSci stromal and immune WT datasets.
#
# Usage, from the project root:
#   env/bin/Rscript scripts/R/build-pansci-azimuth.R
#

library(targets)

# Source utilities
source("src/R/scrna-seq_utilities.R")

# ── Build PanSci Stromal WT reference ──────────────────────────────────────

message(strrep("=", 70))
message("Loading pansci_stromal reference from targets...")
pansci_stromal_ref <- tar_read(pansci_stromal)

message("Building Azimuth reference (stromal)...")
message("  Reference contains ", ncol(pansci_stromal_ref), " cells")
message("  Reference contains ", length(unique(pansci_stromal_ref$cell_type)), " cell types")

# Build Azimuth reference
ref_path_stromal <- build_azimuth_ref(
  srt         = pansci_stromal_ref,
  ref_name    = "pansci_stromal",
  refDR       = "refDR",
  refUMAP     = "refUMAP",
  refAssay    = "refAssay",
  metadata    = "cell_type",
  azimuth_dir = "external_data/azimuth"
)

message("Azimuth reference (stromal) saved to: ", ref_path_stromal)
message("Reference directory contents:")
print(list.files(ref_path_stromal, full.names = TRUE))

# Clean up
rm(pansci_stromal_ref)
gc()

# ── Build PanSci Immune WT reference ───────────────────────────────────────

message(strrep("=", 70))
message("Loading pansci_immune reference from targets...")
pansci_immune_ref <- tar_read(pansci_immune)

message("Building Azimuth reference (immune)...")
message("  Reference contains ", ncol(pansci_immune_ref), " cells")
message("  Reference contains ", length(unique(pansci_immune_ref$cell_type)), " cell types")

# Build Azimuth reference
ref_path_immune <- build_azimuth_ref(
  srt         = pansci_immune_ref,
  ref_name    = "pansci_immune",
  refDR       = "refDR",
  refUMAP     = "refUMAP",
  refAssay    = "refAssay",
  metadata    = "cell_type",
  azimuth_dir = "external_data/azimuth"
)

message("Azimuth reference (immune) saved to: ", ref_path_immune)
message("Reference directory contents:")
print(list.files(ref_path_immune, full.names = TRUE))

message(strrep("=", 70))
message("Done! Built 2 Azimuth references:")
message("  - ", ref_path_stromal)
message("  - ", ref_path_immune)
