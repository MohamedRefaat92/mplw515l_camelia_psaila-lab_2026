#!/usr/bin/env bash
# Copy a rendered report folder to the dated archive, where it is never modified again:
#   output/reports/archive/<YYYY-MM-DD>/<report>_<YYYY-MM-DD_HH-MM>/
# Share or zip the archived folder as a whole (the HTML links to the files next to it).
#
#   tools/snapshot-report.sh <report name or folder>   e.g. example-report
set -euo pipefail
cd "$(dirname "$0")/.."   # run from the project root, wherever this was called from

if [[ $# -ne 1 ]]; then
  echo "usage: tools/snapshot-report.sh <report name or folder under output/reports/latest/>" >&2
  exit 2
fi

src=$1
if [[ ! -d $src ]]; then   # a bare name: find its folder under output/reports/latest/
  src=$(find output/reports/latest -type d -name "$1" -not -name '*_files' 2>/dev/null | head -1)
fi
if [[ -z $src || ! -d $src ]]; then
  echo "no rendered report '$1' under output/reports/latest/ (run tools/run-pipeline.sh first)" >&2
  exit 1
fi

report=$(basename "$src")
stamp=$(date +%Y-%m-%d_%H-%M)
dest="output/reports/archive/${stamp%_*}/${report}_${stamp}"
if [[ -e $dest ]]; then
  echo "already snapshotted this minute: $dest" >&2
  exit 1
fi

mkdir -p "$(dirname "$dest")"
cp -R "$src" "$dest"
echo "snapshot: $dest"
