#!/usr/bin/env bash
#SBATCH --job-name=cellrank
#SBATCH --partition=short,long
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/cellrank_%j.out
#SBATCH --error=output/logs/cellrank_%j.err
#
# CellRank on the HSC compartment (src/python/analysis/cellrank_run.py), into
# output/objects/cellrank/<run tag>/. The input decides what drives the transitions:
#
#   INPUT_TYPE=cytotrace sbatch scripts/slurm/analysis/cellrank.sh                  CytoTRACE2 pseudotime (pipeline's hsc.h5ad)
#   INPUT_TYPE=cytotrace CLUSTER_RES=1.5 sbatch scripts/slurm/analysis/cellrank.sh
#   INPUT_TYPE=dpt CLUSTER_RES=1.5 sbatch scripts/slurm/analysis/cellrank.sh        DPT (notebooks/exploration/dpt_hsc.ipynb first)
#   INPUT_TYPE=velocity VELOCITY_MODE=dynamical sbatch scripts/slurm/analysis/cellrank.sh   (scripts/slurm/velocity/scvelo.sh first)
#
# The cellrank-hsc reports read the runs cytotrace, cytotrace_res1_5 and dpt_res1_5.
# Other settings (N_MACROSTATES, H5AD_PATH, OUT_DIR): the header of the Python file.
set -euo pipefail
source scripts/slurm/env.sh

export INPUT_TYPE=${INPUT_TYPE:-velocity}
export N_JOBS=$SLURM_CPUS_PER_TASK
uv run --no-sync python -m analysis.cellrank_run
