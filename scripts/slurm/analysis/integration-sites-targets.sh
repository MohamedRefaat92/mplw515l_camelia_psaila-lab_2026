#!/usr/bin/env bash
#SBATCH --job-name=is_targets
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/is-targets_%j.out
#SBATCH --error=output/logs/is-targets_%j.err
#
# The four derived targets the scarab reports read, and the two pre-QC objects they need.
#
#   sbatch scripts/slurm/analysis/integration-sites-targets.sh
#
# WHY THESE ARE TARGETS AND NOT A SCRIPT. The reports used to read the processed Seurat objects
# directly - 2 GB each, ten minutes per read - which was tolerable for one report and absurd for
# four. What a report actually needs from them is coordinates, a cluster label, a QC flag, the
# counts matrix and the cell-cell graph, so those are extracted once here and read as small targets.
# Nothing is exported to a flat file: a target is already a cached, provenance-tracked object, and a
# script writing CSVs beside it would be a second source of truth for the same numbers.
#
#   is_cell_embedding  all four objects' UMAP + PCA coordinates, metadata and QC flag, ~MB
#   is_hsc_counts      the RNA counts of the cells that pass QC, kept sparse
#   is_hsc_graph       the SNN graph over those same cells, in the same order
#   is_hsc_loadings    genes x 30 PCA loadings
#
# WHY IT NEEDS 12 HOURS AND 192 GB. mpn_raw_processed_hsc and mpn_raw_processed_mnc are stale: their
# input mpn_raw_hsc was rebuilt on 2026-11-20 at 09:54, after they were built on 2026-11-10, so
# targets re-runs SCTransform on both. That is the expensive part. It is done properly rather than
# read around with shortcut = TRUE, because is_cell_embedding holds pre-QC and post-QC rows in ONE
# table: reading the stale objects would put rows built from a superseded input beside rows built
# from the current one, and three tables in scarab-qc compare exactly those two groups.
#
# shortcut = TRUE is still passed, but only to stop targets walking FURTHER up: mpn_raw_hsc and
# mpn_processed_hsc_fullcc are current, and this job has no business rebuilding the whole project.
set -euo pipefail
source scripts/slurm/env.sh

echo "=== outdated before ==="
env/bin/Rscript -e 'print(targets::tar_outdated(names = c("is_cell_embedding", "is_hsc_counts",
  "is_hsc_graph", "is_hsc_loadings"), reporter = "silent"))'

env/bin/Rscript -e 'targets::tar_make(
  names = c("mpn_raw_processed_hsc", "mpn_raw_processed_mnc",
            "is_cell_embedding", "is_hsc_counts", "is_hsc_graph", "is_hsc_loadings"),
  shortcut = TRUE, reporter = "summary")'

echo "=== what was built ==="
env/bin/Rscript -e '
m <- targets::tar_meta(fields = c("time", "bytes", "seconds", "error"))
n <- c("mpn_raw_processed_hsc", "mpn_raw_processed_mnc", "is_cell_embedding",
       "is_hsc_counts", "is_hsc_graph", "is_hsc_loadings")
print(as.data.frame(m[m$name %in% n, ]), row.names = FALSE)
e <- targets::tar_read(is_cell_embedding)
cat("\nis_cell_embedding:", nrow(e), "rows\n")
print(table(e$panel, e$stage))
stopifnot(sum(grepl("^PC_", names(e))) == 30L)
cat("PC_1..PC_30 present\n")'

echo "=== done ==="
