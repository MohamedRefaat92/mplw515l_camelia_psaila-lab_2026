#!/usr/bin/env bash
# Link each tube's Cell Ranger results (scripts/slurm/cellranger/cellranger-multi.sh) into data/ as
# data/cellranger_<tube>: its per-sample folder (matrix, BAM, web summary, metrics).
# Run once the jobs have finished and their web summaries look right.
#
#   scripts/sh/link-cellranger.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/slurm/cellranger/env.sh

unfinished=0
for tube in $(tail -n +2 "$CR_LIBRARIES" | cut -d, -f1); do
  result=$CR_OUT/$tube/outs/per_sample_outs/$tube
  if [[ ! -d $result/sample_filtered_feature_bc_matrix ]]; then
    echo "not finished: $tube ($result)" >&2
    unfinished=$((unfinished + 1))
  elif [[ -L data/cellranger_$tube ]]; then
    echo "already linked: data/cellranger_$tube"
  else
    tools/link-input.sh "$result" --name "cellranger_$tube"
  fi
done

if [[ $unfinished -gt 0 ]]; then
  echo "$unfinished tube(s) not finished yet: link them once they are" >&2
  exit 1
fi
echo "the first tools/catalogue-data.sh checksums the BAMs too (~430 GB, about an hour)"
