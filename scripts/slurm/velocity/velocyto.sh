#!/usr/bin/env bash
#SBATCH --job-name=velocyto
#SBATCH --partition=short,long
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/velocyto_%A_%a.out
#SBATCH --error=output/logs/velocyto_%A_%a.err
#
# Spliced / unspliced counts for one tube: array index i is tube Rx_i.
# Reads data/cellranger_Rx_i/ (BAM and cell barcodes; scripts/sh/link-cellranger.sh) and writes
# output/objects/velocyto/Rx_i/Rx_i.loom, which scripts/slurm/velocity/scvelo.sh reads.
#
#   sbatch scripts/slurm/velocity/velocyto.sh             all 8 tubes
#   sbatch --array=3 scripts/slurm/velocity/velocyto.sh   one tube
#
# The BAM is first sorted by cell barcode into the output folder; velocyto (patched, see
# scripts/containers/velocyto-run.py) finds it there instead of writing a 53 GB copy next to the
# Cell Ranger output. The sorted copy is deleted once the loom is written.
# Needs output/containers/velocyto.sif: sbatch scripts/slurm/setup/build-container.sh scripts/containers/velocyto.def
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

tube=Rx_${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch scripts/slurm/velocity/velocyto.sh}
cr=data/cellranger_$tube
work=output/objects/velocyto/$tube
if [[ -s $work/$tube.loom ]]; then
  echo "[SKIP] $work/$tube.loom exists"
  exit 0
fi
if [[ ! -f $cr/sample_alignments.bam ]]; then
  echo "ERROR: $cr/sample_alignments.bam missing: link the Cell Ranger results (scripts/sh/link-cellranger.sh)" >&2
  exit 1
fi
mkdir -p "$work"

# velocyto needs an uncompressed GTF: the annotation Cell Ranger aligned against
gtf=$work/genes.gtf
[[ -s $gtf ]] || zcat "$CR_REF/genes/genes.gtf.gz" >"$gtf"

ln -sfn "$(readlink -f "$cr/sample_alignments.bam")" "$work/sample_alignments.bam"
sorted=$work/cellsorted_sample_alignments.bam
if [[ ! -s $sorted ]]; then
  in_container velocyto samtools sort -t CB -O bam -m 4000M -@ "$SLURM_CPUS_PER_TASK" \
    -T "$work/tmp_sort" -o "$sorted.part" "$work/sample_alignments.bam"
  mv "$sorted.part" "$sorted"
fi

in_container velocyto velocyto run --sampleid "$tube" --outputfolder "$work" \
  --bcfile "$cr/sample_filtered_feature_bc_matrix/barcodes.tsv.gz" \
  --samtools-memory 4000 -@ "$SLURM_CPUS_PER_TASK" \
  "$work/sample_alignments.bam" "$gtf"

rm -f "$sorted" "$gtf"
echo "=== $tube done: $work/$tube.loom ==="
