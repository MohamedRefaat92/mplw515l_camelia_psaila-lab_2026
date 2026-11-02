#!/usr/bin/env bash
# Create small mock inputs to try the example pipeline:
#   counts.tsv (data)          written OUTSIDE the project and linked into data/
#   metadata/samples.csv       the sample sheet, written only if it doesn't exist yet
#
#   tools/mock-data.sh [dir]     default dir: ../<project>-mock-data (next to the project folder)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
dir=${1:-"$(dirname "$root")/$(basename "$root")-mock-data"}

cd "$root"
had_sheet=false; [[ -f metadata/samples.csv ]] && had_sheet=true
env/bin/Rscript -e 'source("tools/R/mock-data.R"); write_mock_data(commandArgs(TRUE)[1])' "$dir" >/dev/null
echo "mock counts written to $dir"
if $had_sheet; then
  echo "metadata/samples.csv already exists, left as is"
else
  echo "mock sample sheet written to metadata/samples.csv"
fi
if [[ -e data/counts.tsv || -L data/counts.tsv ]]; then
  echo "data/counts.tsv already exists, left as is"
else
  tools/link-input.sh "$dir/counts.tsv"
fi
tools/catalogue-data.sh
