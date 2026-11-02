#!/usr/bin/env bash
# Run the targets pipeline, then snapshot every report it (re-)rendered into output/reports/archive/.
#
#   tools/run-pipeline.sh
set -euo pipefail
cd "$(dirname "$0")/.."   # run from the project root, wherever this was called from

marker=$(mktemp)          # reports rendered after this moment were produced by this run
trap 'rm -f "$marker"' EXIT

env/bin/Rscript -e 'targets::tar_make()'

rendered=$(find output/reports/latest -name '*.html' -newer "$marker" 2>/dev/null || true)
if [[ -z $rendered ]]; then
  echo "No report was re-rendered (everything up to date); nothing to snapshot."
  exit 0
fi
while IFS= read -r html; do
  tools/snapshot-report.sh "$(dirname "$html")"
done <<<"$rendered"
