#!/usr/bin/env bash
#SBATCH --job-name=cellranger_multi
#SBATCH --partition=short,long,short_el9,long_el9
#SBATCH --array=1-8%4
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=30:00:00
#SBATCH --output=output/logs/cellranger-multi_%A_%a.out
#SBATCH --error=output/logs/cellranger-multi_%A_%a.err
#
# Cell Ranger multi for one tube: array index i is tube Rx_i of metadata/cellranger/libraries.csv.
# Cell Ranger 10.1.0, from output/containers/cellranger.sif (build it first, see that recipe), with
# the setup of the first processing of these data (rio548, January 2026, Cell Ranger 9.0.0): the
# GRCm39 + mCherry + DRE reference, the hashtags read as Antibody Capture with the feature references
# in metadata/cellranger/, BAMs kept (RNA velocity needs them).
#
#   sbatch scripts/slurm/cellranger/cellranger-multi.sh             all 8 tubes, at most 4 at a time
#   sbatch --array=3 scripts/slurm/cellranger/cellranger-multi.sh   one tube
#
# The FASTQs are read where the sequencing provider delivered them ($FASTQ_ROOT); a tube's run writes
# its multi config and then cellranger multi into $CR_OUT/Rx_i/ (resubmitting an interrupted tube
# resumes where it stopped). Cell Ranger 9.0.0 took ~20 h per tube on 32 cores, at most ~120 GB of
# memory, and wrote ~60 GB of results, 53 GB of it the BAM.
# Then link the results into data/: scripts/sh/link-cellranger.sh
# Paths and the Cell Ranger version: scripts/slurm/cellranger/env.sh
set -euo pipefail
source scripts/slurm/env.sh
source scripts/slurm/cellranger/env.sh

tube=Rx_${SLURM_ARRAY_TASK_ID:?submit as an array job: sbatch scripts/slurm/cellranger/cellranger-multi.sh}
row=$(awk -F, -v t="$tube" '$1 == t' "$CR_LIBRARIES")
[[ -n $row ]] || { echo "ERROR: $tube is not in $CR_LIBRARIES" >&2; exit 1; }
IFS=, read -r _ feature_ref gex_id gex_dir cite_id cite_dir <<<"$row"

out=$CR_OUT/$tube
if [[ -d $out/outs/per_sample_outs/$tube/sample_filtered_feature_bc_matrix ]]; then
  echo "[SKIP] $tube: already processed ($out)"
  exit 0
fi

gex_fastqs=$FASTQ_ROOT/$gex_dir
cite_fastqs=$FASTQ_ROOT/$cite_dir
for dir in "$gex_fastqs" "$cite_fastqs"; do
  [[ -d $dir ]] || { echo "ERROR: $tube: no FASTQ folder at $dir" >&2; exit 1; }
done

# The config Cell Ranger is given; it saves the one it used as outs/config.csv. Kept beside the
# results and not inside $out: cellranger multi refuses an --output-dir that exists but is not one
# of its own run folders ("is not a pipestance directory"), so $out must not be created here.
mkdir -p "$CR_OUT/configs"
config=$CR_OUT/configs/$tube.csv
cat >"$config" <<EOF
[gene-expression]
reference,$CR_REF
create-bam,true

[feature]
reference,$PROJ/metadata/cellranger/$feature_ref

[libraries]
fastq_id,fastqs,feature_types
$gex_id,$gex_fastqs,Gene Expression
$cite_id,$cite_fastqs,Antibody Capture
EOF

# Inside the container Cell Ranger sees the project, the results folder, the FASTQs and the
# reference; not the folders behind data/, which include its own (read-only) earlier results
export CONTAINER_BINDS=$PROJ EXTRA_BINDS="$CR_OUT,$FASTQ_ROOT:$FASTQ_ROOT:ro,$CR_REF:$CR_REF:ro"
version=$(in_container cellranger cellranger --version)   # "cellranger 10.1.0"
read -r _ found <<<"$version"
if [[ ${found#cellranger-} != "$CR_VERSION" ]]; then      # 9.x printed "cellranger cellranger-9.0.0"
  echo "ERROR: output/containers/cellranger.sif has $version, not Cell Ranger $CR_VERSION" >&2
  exit 1
fi
rm -f "$out/_lock"   # left behind when an earlier attempt was killed; never run two jobs for one tube
in_container cellranger cellranger multi --id "$tube" --csv "$config" --output-dir "$out" \
  --localcores "$SLURM_CPUS_PER_TASK" --localmem $((SLURM_MEM_PER_NODE / 1024 - 10)) --disable-ui

echo "=== $tube done: $out/outs/per_sample_outs/$tube ==="
