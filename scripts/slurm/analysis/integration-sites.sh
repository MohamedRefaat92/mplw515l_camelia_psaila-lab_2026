#!/usr/bin/env bash
#SBATCH --job-name=integration_sites
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
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
# EDIT THIS FILE BEFORE YOU SUBMIT IT, NOT AFTER. Slurm copies the batch script into its spool at
# submit time, so a job already queued runs the version it was submitted with, while the R files it
# calls are read at run time and so are the current ones. On 2026-11-17 this file gained the genome
# FASTA argument on the assign line six minutes after the job was submitted: the job ran the old
# line, tsd_check() got no reference to compare against, and the sequence half of the strongest check
# in the analysis wrote NA with no error and no warning. Resubmit rather than edit in place.
#
# Measured: 1 h 49 m end to end for eight tubes (job 35742986), 13.7 GB peak.
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
#
# -T 20 is not cosmetic. The default is 30, and with -A 1 a perfect L-base clip scores L, so every
# clip shorter than 30 fell below threshold and was dropped with no error and no record: placement
# rate was exactly 0.00 at 25, 26, 27, 28 and 29 bp -- 244 clips, a quarter of discovery.fa --
# while IS_DEFAULTS$min_discovery advertised 25. Set it explicitly and let MAPQ carry specificity,
# which it does: 42 of the surviving hits already fall below MAPQ 30.
$bwa mem -k 16 -T 20 -t "$threads" "$index" "$work/discovery.fa" > "$work/discovery.sam"

# The 5' channel, same treatment as discovery: these clips are long enough to place on their own.
if [[ -s $work/five_prime.fa ]]; then
  $bwa mem -k 16 -T 20 -t "$threads" "$index" "$work/five_prime.fa" > "$work/five_prime.sam"
fi

echo "=== 3/5 the genome-mapped channel ==="
# The only expensive stage: a streaming pass over every alignment in each BAM, ~104 min a tube, run
# one core per tube. It needs no aligner afterwards -- these reads already carry their host
# coordinate at MAPQ 255 -- and it is kept separate so it can be rerun without redoing the rest.
# It is also the stage that finds the most sites: 13 of the 15 on its own, four of them invisible to
# the 3' channel because their clips align at MAPQ 0.
# shellcheck disable=SC2086  # tubes is a deliberate word-split list
$rscript scripts/R/call-integration-sites.R genomic "$out" "$PROJ/data" $tubes

# The genome-mapped channel's error bound. Commented out because it repeats the most expensive
# stage in the pipeline -- the decoy pass costs what the real pass costs, so enabling it roughly
# doubles the job -- and because a control must never be on the critical path of a run that someone
# is waiting for. Uncomment, or run the stage on its own, when the bound needs tightening:
#
#   sbatch --time=12:00:00 --wrap "scripts/R/call-integration-sites.R genomic-decoy $out $PROJ/data"
#
# Measured so far: chr19 in Rx_4, 27.6M alignments containing a real 29-read site, zero decoy hits.
# $rscript scripts/R/call-integration-sites.R genomic-decoy "$out" "$PROJ/data" chr19 "$CR_REF/fasta/genome.fa"

echo "=== 4/5 call the sites, and build the flanks plus the decoys ==="
$rscript scripts/R/call-integration-sites.R sites "$out" "$gtf" "$CR_REF/fasta/genome.fa"

echo "=== 5/5 assign clips, cell table, TSD, decoy error rate ==="
# bwa no longer runs against the flanks. It contributed nothing here -- MAPQ against a handful of
# short sequences is near-meaningless and was already ignored, -a was needed to stop its ranking
# discarding true hits, and -T 12 silently voided the advertised 8 bp floor (placement rate 0.00 at
# 8, 9, 10 and 11 bp, 169 reads written and dropped). assign_clips_exact() does the positional test
# directly: exhaustive, deterministic, seconds, and the floor is now whatever min_assign says.
# min_assign is 10, set by the decoy flanks: zero decoy hits at 10 bp and above, two below it.
$rscript scripts/R/call-integration-sites.R assign "$out" "$CR_REF/fasta/genome.fa"

echo "=== done ==="
ls -l "$out"/*.csv
