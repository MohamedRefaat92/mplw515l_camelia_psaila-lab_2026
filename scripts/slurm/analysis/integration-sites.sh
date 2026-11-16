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

echo "=== 1/5 extract the clipped tails ==="
# shellcheck disable=SC2086  # tubes is a deliberate word-split list
$rscript scripts/R/call-integration-sites.R extract "$out" "$PROJ/data" $tubes

echo "=== 2/5 place them on the genome ==="
# bwa mem with a shortened seed: a 25-48 base query will not seed at the default 19.
$bwa mem -k 16 -t "$threads" "$index" "$work/discovery.fa" > "$work/discovery.sam"

# The 5' channel, same treatment as discovery: these clips are long enough to place on their own.
if [[ -s $work/five_prime.fa ]]; then
  $bwa mem -k 16 -t "$threads" "$index" "$work/five_prime.fa" > "$work/five_prime.sam"
fi

echo "=== 3/5 call the sites and build the mini-reference ==="
$rscript scripts/R/call-integration-sites.R sites "$out" "$gtf" "$CR_REF/fasta/genome.fa"

echo "=== 4/5 assign clips against the mini-reference ==="
# A separate, tiny index: the called flanks only. -k 8 -T 12 because the search space is a few kb
# rather than 2.7 Gb, so a far shorter clip is decisive here -- and because the real discriminator
# is position, not score, read_mini_hits() requires the clip to land on the junction base exactly.
$bwa index "$work/site_flanks.fa" 2>/dev/null
$bwa mem -k 8 -T 12 -a -t "$threads" "$work/site_flanks.fa" "$work/assign.fa" > "$work/assign.sam"

echo "=== 5/5 cell assignments and the TSD check ==="
$rscript scripts/R/call-integration-sites.R assign "$out"

echo "=== done ==="
ls -l "$out"/*.csv
