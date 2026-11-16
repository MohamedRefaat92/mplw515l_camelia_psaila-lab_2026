#!/usr/bin/env bash
#SBATCH --job-name=integration_sites
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=output/logs/integration-sites_%j.out
#SBATCH --error=output/logs/integration-sites_%j.err
#
# Lentiviral integration sites from read-through transcripts.
#
#   sbatch scripts/slurm/analysis/integration-sites.sh
#
# Transcription of the cassette does not stop cleanly at the polyA site in the 3' LTR's R region: a
# minority of transcripts read through to the end of the provirus and on into flanking host DNA. Such
# a read aligns to the reporter contig as far as its last base and soft-clips the rest, and that
# clipped tail is host sequence naming where the vector integrated.
#
# Three steps, because only the middle one needs bwa:
#   1. extract  src/R/integration-sites.R reads each tube's BAM and writes the clipped tails
#   2. align    bwa places those tails on the genome
#   3. call     the same R code clusters the placements into sites and writes the tables
#
# The reasoning behind each filter and threshold lives in src/R/integration-sites.R, next to the code
# that applies it. Every threshold is in IS_DEFAULTS there; override one by editing it, so the job
# file and the analysis cannot disagree.
#
# Stays on short -- NOT the _el9 partitions. renv keys its library path on the OS release
# (renv/library/linux-rocky-8.10/...), so on a Rocky 9 node it finds no library, bootstraps an empty
# one and the first library() call fails. Anything here that uses ./env must stay on EL8 nodes.
# short and long are the same 115 nodes and differ only in their limit (30 h vs 10 days), so short
# is the whole pool for a job this size.
#
# Needs the bwa index: sbatch scripts/slurm/setup/bwa-index.sh (once, ~50 min).
# Reads data/cellranger_<tube>/sample_alignments.bam, so it is covered by the project's provenance.
# Writes output/objects/integration-sites/, which notebooks/integration-sites/ reads.
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

bwa=$PROJ/env/bin/bwa
rscript=$PROJ/env/bin/Rscript
index=$CR_REF/bwa/genome
gtf=$CR_REF/genes/genes.gtf.gz
out=$PROJ/output/objects/integration-sites
work=$out/work
threads=${SLURM_CPUS_PER_TASK:-8}

for f in "$index.bwt" "$gtf"; do
  [[ -s $f ]] || { echo "ERROR: $f missing (run scripts/slurm/setup/bwa-index.sh first)" >&2; exit 1; }
done

tubes=$(tail -n +2 "$CR_LIBRARIES" | cut -d, -f1 | tr '\n' ' ')

echo "=== 1/3 extract the clipped tails ==="
# shellcheck disable=SC2086  # tubes is a deliberate word-split list
$rscript scripts/R/call-integration-sites.R extract "$out" "$PROJ/data" $tubes

echo "=== 2/3 place them on the genome ==="
# bwa mem with a shortened seed: a 25-48 base query will not seed at the default 19.
$bwa mem -k 16 -t "$threads" "$index" "$work/discovery.fa" > "$work/discovery.sam"

# bwa aln for the short tier: it was written for reads this length and reports calibrated MAPQ
# there, where mem becomes erratic. -n 1 allows one mismatch, -l 15 seeds the whole read.
if [[ -s $work/support.fa ]]; then
  $bwa aln -n 1 -o 0 -l 15 -t "$threads" "$index" "$work/support.fa" > "$work/support.sai"
  $bwa samse "$index" "$work/support.sai" "$work/support.fa" > "$work/support.sam"
fi

echo "=== 3/3 call the sites ==="
$rscript scripts/R/call-integration-sites.R call "$out" "$gtf"

echo "=== done ==="
ls -l "$out"/*.csv
