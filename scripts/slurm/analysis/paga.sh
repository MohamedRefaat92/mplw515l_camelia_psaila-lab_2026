#!/usr/bin/env bash
#SBATCH --job-name=paga
#SBATCH --partition=short,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=output/logs/paga_%j.out
#SBATCH --error=output/logs/paga_%j.err
#
# PAGA connectivity of the HSC clusters, overall and per condition (src/python/analysis/run_paga.py),
# into output/objects/paga/ and output/objects/paga_conditions/.
#
#   RESOLUTION=0.9 sbatch scripts/slurm/analysis/paga.sh
#
# The paga-traj-analysis-* reports read resolutions 0.9, 1 and 1.5 (written 1, not 1.0: it is
# part of the cluster column's name, sNN_pca.SCT_30_res.1).
set -euo pipefail
source scripts/slurm/env.sh

export RESOLUTION=${RESOLUTION:-0.9}
uv run --no-sync python -m analysis.run_paga
