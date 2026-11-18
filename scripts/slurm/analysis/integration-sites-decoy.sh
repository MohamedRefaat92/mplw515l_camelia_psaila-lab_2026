#!/usr/bin/env bash
#SBATCH --job-name=is_decoy
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=output/logs/integration-sites-decoy_%j.out
#SBATCH --error=output/logs/integration-sites-decoy_%j.err
#
# Error bound for the genome-mapped channel, genome-wide.
#
#   sbatch scripts/slurm/analysis/integration-sites-decoy.sh            # whole genome
#   sbatch scripts/slurm/analysis/integration-sites-decoy.sh chr19      # one chromosome
#
# The 100 decoy flanks bound the ASSIGNMENT tier. They say nothing about the channel that calls four
# of the fifteen sites on its own and carries IS10's 99 cells, whose claim is different: that a read
# the aligner placed at MAPQ 255, whose soft-clipped end matches the last min_genomic bases of the
# provirus terminus within the absorption window, is a junction read. This measures how often that is
# true of a terminus that cannot be there -- the same 13 bases permuted to the same base composition.
#
# Why it is a separate job rather than a stage of the main one: it costs what the genomic scan costs,
# and a control must not sit on the critical path of a run someone is waiting for.
#
# Cost, from the main job's measurement (35742986: 8 BAMs forked, 1 h 49 m, 13.7 GB peak): one
# whole-BAM pass is ~105 min a tube, so forked one core per tube this is ~1 h 50 m. It is NOT double
# that, because the real counts are read from the genomic stage's cache rather than rescanned -- so
# run the main job first, or this falls back to rescanning and takes twice as long.
#
# Stays on short and NOT the _el9 partitions, for the renv library-path reason in the main job file.
#
# Measured so far, chr19 in Rx_4: 27.6M alignments containing a real 29-read site, zero decoy hits.
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

rscript=$PROJ/env/bin/Rscript
out=$PROJ/output/objects/integration-sites
region=${1:-}

[[ -s $out/work/tubes.txt ]] || {
  echo "ERROR: $out/work/tubes.txt missing - run the extract stage first" >&2; exit 1; }
if [[ ! -s $out/work/genomic_Rx_4.csv ]]; then
  echo "WARNING: no genomic cache; the real pass will be rescanned and this will take ~2x longer" >&2
fi

echo "=== scrambled-terminus decoy over the genome-mapped channel ${region:+for $region} ==="
$rscript scripts/R/call-integration-sites.R genomic-decoy \
  "$out" "$PROJ/data" "$region" "$CR_REF/fasta/genome.fa"

echo "=== done ==="
cat "$out/genomic_decoy_rate.csv"
