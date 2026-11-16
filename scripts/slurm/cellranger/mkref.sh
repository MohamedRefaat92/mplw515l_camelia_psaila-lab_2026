#!/usr/bin/env bash
#SBATCH --job-name=cellranger_mkref
#SBATCH --partition=short,long,short_el9,long_el9
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/cellranger-mkref_%j.out
#SBATCH --error=output/logs/cellranger-mkref_%j.err
#
# Build the reference cellranger multi runs against ($CR_REF): the 10x prebuilt mouse reference
# ($CR_REF_SRC) with the reporter cassette ($CR_REPORTER.fa and .gtf) appended as one extra contig.
# Cell Ranger 10.1.0, from output/containers/cellranger.sif (build it first, see that recipe).
#
#   sbatch scripts/slurm/cellranger/mkref.sh
#
# Appended, not passed as a second --genome: a second genome would prefix every gene name with it
# (GRCm39_Actb, reporter_mCherry). Takes ~1-2 h, almost all of it the STAR index.
# Skipped if $CR_REF exists, so delete it to rebuild. Paths: scripts/slurm/cellranger/env.sh
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

if [[ -d $CR_REF ]]; then
  echo "[SKIP] $CR_REF exists (delete it to rebuild)"
  exit 0
fi
for f in "$CR_REF_SRC/fasta/genome.fa" "$CR_REF_SRC/genes/genes.gtf.gz" "$CR_REPORTER.fa" "$CR_REPORTER.gtf"; do
  [[ -s $f ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# The cassette's contig name, which its GTF's first column must match
contig=$(sed -n '1s/^>\([^[:space:]]*\).*/\1/p' "$CR_REPORTER.fa")
awk -F'\t' -v c="$contig" '!/^#/ && $1 != c { print "ERROR: " FILENAME " line " FNR " is on " $1 ", not " c > "/dev/stderr"; exit 1 }' "$CR_REPORTER.gtf"

# mkref writes its logs into the working directory, so build in a scratch folder next to the result
work=$(dirname "$CR_REF")/.mkref.$SLURM_JOB_ID
mkdir -p "$work"
trap 'rm -rf "$work"' EXIT

echo ">>> joining the 10x reference and $contig"
cat "$CR_REF_SRC/fasta/genome.fa" "$CR_REPORTER.fa" >"$work/genome.fa"
zcat "$CR_REF_SRC/genes/genes.gtf.gz" >"$work/genes.gtf"
grep -v '^#' "$CR_REPORTER.gtf" >>"$work/genes.gtf"

export EXTRA_BINDS="$(dirname "$(dirname "$CR_REF_SRC")"),$(dirname "$CR_REPORTER")"
cd "$work"
in_container cellranger cellranger mkref \
  --genome "$(basename "$CR_REF")" \
  --fasta "$work/genome.fa" \
  --genes "$work/genes.gtf" \
  --output-dir "$CR_REF" \
  --nthreads "$SLURM_CPUS_PER_TASK" \
  --memgb $((SLURM_MEM_PER_NODE / 1024 - 8))

# The cassette has to have survived into the built reference, or the run would silently lose it
grep -qw "$contig" "$CR_REF/star/chrName.txt" || { echo "ERROR: $contig is not in $CR_REF" >&2; exit 1; }
echo "=== built $CR_REF (contigs: $(wc -l <"$CR_REF/star/chrName.txt")) ==="
