#!/usr/bin/env bash
#SBATCH --job-name=pipeline
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=5-00:00:00
#SBATCH --output=output/logs/pipeline_%j.out
#SBATCH --error=output/logs/pipeline_%j.err
#
# The targets pipeline and its reports (tools/run-pipeline.sh) as a cluster job: the Seurat
# processing, the reference projections and the reports need more memory and time than an
# interactive session. targets resumes where it stopped, so resubmitting continues the run.
#
#   sbatch scripts/slurm/pipeline/run-pipeline.sh
#
# A target that fails doesn't stop the others (error = "continue" in _targets.R). Afterwards:
#   env/bin/Rscript -e 'targets::tar_errored()'           what failed
#   env/bin/Rscript -e 'targets::tar_meta(fields = error)' why
set -euo pipefail
source scripts/slurm/env.sh

tools/run-pipeline.sh
