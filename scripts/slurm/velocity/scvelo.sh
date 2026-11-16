#!/usr/bin/env bash
#SBATCH --job-name=scvelo
#SBATCH --partition=short,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=output/logs/scvelo_%j.out
#SBATCH --error=output/logs/scvelo_%j.err
#
# scVelo on the HSC compartment (src/python/analysis/velocity_scvelo.py): the pipeline's
# output/objects/hsc.h5ad plus the velocyto looms (scripts/slurm/velocity/velocyto.sh), into
# output/objects/velocity/. Its hsc_velocity_<mode>.h5ad feeds scripts/slurm/analysis/cellrank.sh.
#
#   sbatch scripts/slurm/velocity/scvelo.sh                          dynamical model
#   VELOCITY_MODE=stochastic sbatch scripts/slurm/velocity/scvelo.sh
# Other settings (ROOT_CLUSTER, NEGATE_VELOCITY, ...): the header of the Python file.
set -euo pipefail
source scripts/slurm/env.sh

export VELOCITY_MODE=${VELOCITY_MODE:-dynamical}
export N_JOBS=$SLURM_CPUS_PER_TASK
uv run --no-sync python -m analysis.velocity_scvelo
