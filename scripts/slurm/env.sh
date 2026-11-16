# shellcheck shell=bash
# Shared settings for the job files in scripts/slurm/: each one sources this first.
# Jobs are submitted from the project root; their #SBATCH log paths are relative to it.
if [[ ! -f _targets.R || ! -d scripts/slurm ]]; then
  echo "ERROR: submit from the project root (cd there, then sbatch scripts/slurm/...)" >&2
  exit 1
fi
PROJ=$(pwd -P)
export PROJ
mkdir -p output/logs   # SLURM can't create the log folder itself

export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$HOME/.apptainer/cache}
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-$HOME/.apptainer/tmp}
mkdir -p "$APPTAINER_TMPDIR"

# Folders a container may see: the project (read-write) and the folder behind each link in data/
# (read-only), derived from the links so relinking an input needs no edit here.
CONTAINER_BINDS=$PROJ
while IFS= read -r target; do
  CONTAINER_BINDS+=",$target:$target:ro"
done < <(for link in data/*; do
           [[ -L $link ]] || continue
           target=$(readlink -f "$link")
           [[ -d $target ]] || target=$(dirname "$target")
           echo "$target"
         done | sort -u)
export CONTAINER_BINDS

# in_container <image> <command ...>: run a command inside output/containers/<image>.sif.
# EXTRA_BINDS (comma-separated) adds folders outside the project and data/, e.g. a reference.
in_container() {
  local sif=$PROJ/output/containers/$1.sif
  shift
  if [[ ! -s $sif ]]; then
    echo "ERROR: $sif missing: sbatch scripts/slurm/setup/build-container.sh scripts/containers/$(basename "$sif" .sif).def" >&2
    return 1
  fi
  apptainer exec --cleanenv --bind "$CONTAINER_BINDS${EXTRA_BINDS:+,$EXTRA_BINDS}" "$sif" "$@"
}
