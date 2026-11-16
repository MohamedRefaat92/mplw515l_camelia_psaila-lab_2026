#!/usr/bin/env bash
#SBATCH --job-name=scenic
#SBATCH --partition=short,long
#SBATCH --array=1-4
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=30:00:00
#SBATCH --output=output/logs/scenic_%A_%a.out
#SBATCH --error=output/logs/scenic_%A_%a.err
#
# pySCENIC on the HSC compartment (src/python/analysis/scenic_run.py), one array task per subset:
#   1 mch_neg   2 mch_pos   3 wt_ctrl   4 full (all cells)
# Reads the pipeline's output/objects/hsc.h5ad and the mm10 cisTarget resources in data/scenic/
# (scripts/sh/get-scenic-resources.sh); writes output/objects/scenic/, which the pipeline's scenic_*
# targets and the scenic-hsc reports read.
#
#   sbatch scripts/slurm/analysis/scenic.sh             all four subsets
#   sbatch --array=4 scripts/slurm/analysis/scenic.sh   one subset
#
# Step 1 (h5ad → loom) runs in the project's Python environment, steps 2-7 in the pySCENIC
# container: sbatch scripts/slurm/setup/build-container.sh scripts/containers/pyscenic.def (once).
# A step whose output exists is skipped, so resubmitting resumes.
set -euo pipefail
source scripts/slurm/env.sh

subsets=(mch_neg mch_pos wt_ctrl full)
subset=${subsets[SLURM_ARRAY_TASK_ID - 1]:?array index must be 1-4}
export N_JOBS=$SLURM_CPUS_PER_TASK

uv run --no-sync python -m analysis.scenic_run prepare "$subset"
in_container pyscenic env N_JOBS="$N_JOBS" python src/python/analysis/scenic_run.py infer "$subset"
