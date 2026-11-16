#!/usr/bin/env Rscript
# RETIRED: the pansci_stromal / pansci_immune targets this reads are commented out in _targets.R
# (and so are their projections); kept for when the PanSci references come back.
#
# Test loading PanSci immune RDS file
#

cat("Testing PanSci immune loading...\n")
cat("R version:", R.version.string, "\n")

# Check memory
cat("\nMemory info:\n")
system("free -h")

cat("\nFile info:\n")
file_path <- "external_data/pansci_immune_wt.rds"
cat("  Path:", file_path, "\n")
cat("  Exists:", file.exists(file_path), "\n")
cat("  Size:", round(file.size(file_path) / 1024^3, 2), "GB\n")

cat("\nAttempting to load RDS...\n")
tryCatch({
  srt <- readRDS(file_path)
  cat("SUCCESS: Loaded RDS\n")
  cat("  Class:", class(srt), "\n")
  cat("  Dims:", dim(srt), "\n")
  if ("Sub_cell_type_organ" %in% colnames(srt@meta.data)) {
    cat("  Has Sub_cell_type_organ column: YES\n")
    cat("  Unique cell types:", length(unique(srt$Sub_cell_type_organ)), "\n")
  } else {
    cat("  Has Sub_cell_type_organ column: NO\n")
    cat("  Available columns:", paste(head(colnames(srt@meta.data), 20), collapse = ", "), "\n")
  }
}, error = function(e) {
  cat("ERROR loading RDS:\n")
  cat("  ", conditionMessage(e), "\n")
  traceback()
})

cat("\nDone.\n")
