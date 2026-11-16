#!/usr/bin/env bash
#SBATCH --job-name=bwa_index
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=output/logs/bwa-index_%j.out
#SBATCH --error=output/logs/bwa-index_%j.err
#
# bwa index of the Cell Ranger reference's genome, for the integration-site analysis: the soft-clipped
# tails of reads running off the 3' end of the reporter contig are realigned against it to find where
# the vector integrated.
#
#   sbatch scripts/slurm/setup/bwa-index.sh
#
# Indexes $CR_REF/fasta/genome.fa, which is the mouse genome *plus* the reporter contig, so a clip
# landing on the transgene is recognisable as vector rather than quietly mistaken for host sequence.
# Written to $CR_REF/bwa/genome.* so it travels with the reference it indexes and is rebuilt with it.
# Skipped if the index exists. Takes ~1 h, mostly single-threaded; the index is ~5.5 GB.
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

fasta=$CR_REF/fasta/genome.fa
prefix=$CR_REF/bwa/genome

[[ -s $fasta ]] || { echo "ERROR: no genome at $fasta (run mkref.sh first)" >&2; exit 1; }
if [[ -s $prefix.bwt ]]; then
  echo "[SKIP] $prefix.* exists (delete it to rebuild)"
  exit 0
fi

mkdir -p "$(dirname "$prefix")"
echo ">>> indexing $fasta -> $prefix.*"
"$PROJ"/env/bin/bwa index -p "$prefix" "$fasta"   # explicit path: jobs do not activate ./env

for ext in amb ann bwt pac sa; do
  [[ -s $prefix.$ext ]] || { echo "ERROR: $prefix.$ext missing after the build" >&2; exit 1; }
done
echo "=== done: $(du -sh "$(dirname "$prefix")" | cut -f1) in $(dirname "$prefix") ==="
