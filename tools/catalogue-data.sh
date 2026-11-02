#!/usr/bin/env bash
# Catalogue every input file under data/ (following linked folders) into
# metadata/auto_catalogue_data.tsv: machine-generated and git-ignored.
#
# It starts from your reviewed copy, metadata/catalogue_data.tsv (committed): your descriptions
# are kept, file facts (path, format, size, owner, md5, ...) refreshed, new files added, vanished
# ones marked "missing". Review the generated file, then accept it:
#   cp metadata/auto_catalogue_data.tsv metadata/catalogue_data.tsv   (describe new rows, commit)
#
#   tools/catalogue-data.sh            full scan; md5 only for new or changed files
#   tools/catalogue-data.sh --no-md5   quick scan without checksums
set -euo pipefail
cd "$(dirname "$0")/.."   # run from the project root, wherever this was called from

md5=TRUE
case "${1:-}" in
  --no-md5) md5=FALSE ;;
  "")       ;;
  *)        echo "usage: tools/catalogue-data.sh [--no-md5]" >&2; exit 2 ;;
esac

# shellcheck disable=SC2016  # single-quoted R code; $ is R's list operator
env/bin/Rscript -e '
  source("tools/R/catalogue.R")
  r <- update_catalogue(md5 = as.logical(commandArgs(TRUE)[1]))
  t <- r$table
  cat(sprintf("metadata/auto_catalogue_data.tsv: %d file(s), %d present, %d missing\n",
              nrow(t), sum(t$status == "present"), sum(t$status == "missing")))
  show <- function(label, x) if (length(x)) cat(label, paste(x, collapse = "\n             "), "\n")
  if (!r$reviewed) {
    cat("  no reviewed copy yet (metadata/catalogue_data.tsv)\n")
  } else {
    cat("  compared with metadata/catalogue_data.tsv:\n")
    show("    new:     ", r$added)
    show("    changed: ", r$changed)
    show("    missing: ", r$missing)
    if (!length(c(r$added, r$changed, r$missing))) cat("    same files and checksums\n")
  }
  todo <- sum(t$status == "present" & t$kind == "")
  if (todo) cat(sprintf("  %d file(s) without a description (kind, sample, ...)\n", todo))
  cat("review it, then: cp metadata/auto_catalogue_data.tsv metadata/catalogue_data.tsv\n")
' "$md5"
