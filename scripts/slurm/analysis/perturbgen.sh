#!/usr/bin/env bash
#SBATCH --job-name=perturbgen
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=output/logs/perturbgen_%j.out
#SBATCH --error=output/logs/perturbgen_%j.err
#
# PerturbGen in-silico TF perturbations on the HSC compartment (src/python/analysis/perturbgen_run.py),
# into output/objects/perturbgen/. Reads the pipeline's hsc.h5ad and mm_hs_ortho_long.csv, the SCENIC
# AUC matrix (scripts/slurm/analysis/scenic.sh) and the PerturbGen folder linked as data/perturbgen.
#
#   sbatch scripts/slurm/analysis/perturbgen.sh
#   TOP_N_TFS=5 EPOCHS_DEC=3 sbatch scripts/slurm/analysis/perturbgen.sh
#
# Its Python packages are a separate dependency group: uv sync --group perturbgen (once, and again
# after any plain `uv sync`, which removes them).
set -euo pipefail
source scripts/slurm/env.sh

export WANDB_MODE=disabled WANDB_SILENT=true
export N_WORKERS=${N_WORKERS:-$SLURM_CPUS_PER_TASK}
uv run --no-sync python -m analysis.perturbgen_run
