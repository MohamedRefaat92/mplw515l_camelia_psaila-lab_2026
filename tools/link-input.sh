#!/usr/bin/env bash
# Make input files or folders available to the project as symlinks in data/.
# Each pipeline run records the inputs (resolved path, format, owner, size, time) in its reports.
#
#   tools/link-input.sh <path> [<path> ...]      each linked under its own name (e.g. /data/run1/*.fastq.gz)
#   tools/link-input.sh <path> --name <name>     one path, linked under another name
# A path that can't be linked (missing, inside the project, name taken) is reported and skipped.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
usage="usage: tools/link-input.sh <path> [<path> ...]   or   tools/link-input.sh <path> --name <name>"

paths=() name=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --name) [[ $# -ge 2 ]] || { echo "$usage" >&2; exit 2; }; name=$2; shift 2 ;;
    *)      paths+=("$1"); shift ;;
  esac
done
if [[ ${#paths[@]} -eq 0 || ( -n $name && ${#paths[@]} -ne 1 ) ]]; then
  echo "$usage" >&2
  exit 2
fi

linked=0 skipped=0
for path in "${paths[@]}"; do
  if [[ ! -e $path ]]; then
    echo "skipped: no such file or folder: $path" >&2; skipped=$((skipped + 1)); continue
  fi
  target="$(cd "$(dirname "$path")" && pwd -P)/$(basename "$path")"   # absolute, so the link works from anywhere
  link_name=${name:-$(basename "$path")}
  link="$root/data/$link_name"
  if [[ $target == "$root"/* ]]; then
    echo "skipped: $target is inside the project; inputs should live outside it" >&2; skipped=$((skipped + 1)); continue
  fi
  if [[ -e $link || -L $link ]]; then
    echo "skipped: data/$link_name already exists (-> $(readlink "$link" || echo '?')); remove it first, or use --name" >&2
    skipped=$((skipped + 1)); continue
  fi
  ln -s "$target" "$link"
  echo "data/$link_name -> $target"
  linked=$((linked + 1))
done

echo "$linked linked, $skipped skipped"
if [[ $linked -gt 0 ]]; then
  echo "next: tools/catalogue-data.sh, review metadata/auto_catalogue_data.tsv, copy it to metadata/catalogue_data.tsv"
fi
[[ $skipped -eq 0 ]]   # exit status 1 if anything was skipped
