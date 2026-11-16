#!/usr/bin/env bash
# Download the mouse cisTarget resources SCENIC needs (mm10 refseq_r80, motif collection
# v10nr_clust) into a folder outside the project if they're missing, then link the folder as
# data/scenic (scripts/slurm/analysis/scenic.sh reads it). A copy already in the data hub is linked by
# scripts/sh/link-references.sh instead.
#
#   scripts/sh/get-scenic-resources.sh /well/ludwig/users/nrs939/data_hub/scenic
#
# A download that stops halfway is resumed on the next run.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd -P)
cd "$root"
dest=${1:?usage: scripts/sh/get-scenic-resources.sh <folder outside the project>}
mkdir -p "$dest"
dest=$(cd "$dest" && pwd -P)
if [[ $dest == "$root" || $dest == "$root"/* ]]; then
  echo "error: $dest is inside the project; inputs should live outside it" >&2
  exit 2
fi

base=https://resources.aertslab.org/cistarget
db=$base/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based
urls=(
  "$db/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "$db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
  "$base/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"   # must match the motif collection above
  "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/mm_mgi_tfs.txt"
)

for url in "${urls[@]}"; do
  f=$(basename "$url")
  if [[ -s $dest/$f ]]; then
    echo "present:     $dest/$f"
  else
    # into .part first, so a file under its real name is always complete; -C - resumes a partial one
    echo "downloading: $f"
    curl -fL -C - --retry 5 --retry-delay 10 -o "$dest/$f.part" "$url"
    mv "$dest/$f.part" "$dest/$f"
  fi
done

if [[ -L data/scenic ]]; then
  echo "linked:      data/scenic -> $(readlink data/scenic)"
else
  tools/link-input.sh "$dest" --name scenic
fi
