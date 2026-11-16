#!/usr/bin/env bash
#SBATCH --job-name=build_container
#SBATCH --partition=short,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=output/logs/build-container_%j.out
#SBATCH --error=output/logs/build-container_%j.err
#
# Build an Apptainer image from a recipe into output/containers/<recipe name>.sif (skipped if it
# exists). Arguments after the recipe go to apptainer build (e.g. the recipe's --build-arg).
#
#   sbatch scripts/slurm/setup/build-container.sh scripts/containers/pyscenic.def
#   sbatch scripts/slurm/setup/build-container.sh scripts/containers/velocyto.def
#   sbatch scripts/slurm/setup/build-container.sh scripts/containers/cellranger.def \
#     --build-arg CELLRANGER_TAR=/path/to/cellranger-10.1.0.tar.gz
#
# Built into a temporary file first, so a failed build never leaves a half-written image behind;
# the temporary file is deleted if the build fails. The recipe's %test runs at the end of the build.
set -euo pipefail
source scripts/slurm/env.sh

def=${1:?usage: sbatch scripts/slurm/setup/build-container.sh <recipe.def> [apptainer build options]}
shift
sif=output/containers/$(basename "$def" .def).sif
if [[ -s $sif ]]; then
  echo "[SKIP] $sif exists (delete it to rebuild)"
  exit 0
fi

mkdir -p output/containers
tmp=$sif.building.$SLURM_JOB_ID
trap 'rm -f "$tmp"' EXIT   # removes a failed build's leftover; after the mv there is nothing to remove
apptainer build --force "$@" "$tmp" "$def"
mv "$tmp" "$sif"

ls -lh "$sif"
apptainer inspect "$sif" | grep -E "base.digest|deffile.from"
