#!/usr/bin/env Rscript
# Record every table already in an output directory in its manifest.tsv.
#
#   env/bin/Rscript scripts/R/record-output-manifest.R output/objects/integration-sites
#
# The pipeline writes manifest.tsv as it goes (write_csv_checked() -> record_in_manifest()), so this
# is not part of a normal run. It exists for two cases:
#
#   - outputs produced before the manifest existed, which is why it was written: the tables in
#     output/objects/integration-sites/ are from job 35742986 and predate record_in_manifest();
#   - a manifest lost or truncated while the tables themselves are intact.
#
# What it attests is weaker than what the pipeline attests, and the difference matters. The pipeline
# records a file at the moment it writes it, having just counted its rows; this records whatever is
# on disk now. It can therefore bless a file that was already wrong. Run it only on outputs whose
# provenance you know by other means - here, the clip-table gate passed on these files and every one
# of them was written through write_csv_checked() - and re-run the pipeline rather than this script
# if there is any doubt.
suppressMessages({
  library(readr)
})
source("src/R/integration-sites.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: record-output-manifest.R <output directory>", call. = FALSE)
}
dir <- args[[1]]
if (!dir.exists(dir)) stop("no such directory: ", dir, call. = FALSE)

files <- list.files(dir, pattern = "\\.csv$|\\.csv\\.gz$", full.names = TRUE)
if (!length(files)) stop("no CSVs in ", dir, call. = FALSE)

# The row count is read back from the file rather than passed in, because there is no frame in
# memory here - that is exactly the weakness noted above, and it is recorded as the source so the
# manifest says how each entry was obtained.
for (f in files) {
  rows <- length(readLines(f, warn = FALSE)) - 1L
  record_in_manifest(f, rows)
  message(sprintf("  %-40s %7d rows", basename(f), rows))
}

man <- file.path(dir, "manifest.tsv")
message("\nwrote ", man, " (", length(files), " files)")
# colClasses, not the default: a commit hash is often a valid number. The short form of 8e55518 is
# read as 8 x 10^55518, i.e. Inf, which is how this printout first came out.
shown <- utils::read.delim(man, colClasses = "character", na.strings = character())
shown$commit <- substr(shown$commit, 1, 8)
print(shown[, c("file", "rows", "bytes", "commit")], row.names = FALSE)
