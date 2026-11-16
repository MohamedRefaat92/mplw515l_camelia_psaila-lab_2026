#!/usr/bin/env bash
#SBATCH --job-name=install_r_packages
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=output/logs/install-r-packages_%j.out
#SBATCH --error=output/logs/install-r-packages_%j.err
#
# Install the analysis's R packages into the renv library and record them in renv.lock
# (scripts/R/install-r-packages.R). Once, after tools/bootstrap.sh; compiling takes several hours.
#
#   sbatch scripts/slurm/setup/install-r-packages.sh
#
# 13 packages come from GitHub: a token avoids its rate limit (export GITHUB_PAT=... before sbatch).
set -euo pipefail
source scripts/slurm/env.sh

env/bin/Rscript scripts/R/install-r-packages.R
