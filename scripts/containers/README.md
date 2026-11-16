# scripts/containers/

Apptainer / Singularity definition files (`.def`) for tools conda doesn't provide (e.g. Cell Ranger).

The `.def` recipe is committed; the built image (`.sif`) is large and machine-specific, so it is
git-ignored. Here the images are built into `output/containers/<name>.sif`: like everything else the
project generates, they live in `output/`, and can be rebuilt from their recipe.

These are for tools whose dependencies don't fit the project's environment. Build one as a job,
giving its recipe: `sbatch scripts/slurm/setup/build-container.sh scripts/containers/<name>.def`.
The job files run them through `in_container` (`scripts/slurm/env.sh`).

| Recipe | Image | Used by |
| --- | --- | --- |
| `pyscenic.def` | pySCENIC 0.12.1 (aertslab's image, pinned by digest) | `scripts/slurm/analysis/scenic.sh` |
| `velocyto.def` | velocyto 0.17.17 + samtools on Python 3.10 (base pinned by digest); installs `velocyto-run.py`, velocyto's run command with one fix (see its header) | `scripts/slurm/velocity/velocyto.sh` |
| `cellranger.def` | Cell Ranger 10.1.0 on Debian (base pinned by digest), from a release you download | `scripts/slurm/cellranger/cellranger-multi.sh` |

Cell Ranger: 10x distributes it only from its download page, behind its licence, so download the
release (`.tar.gz`, from https://www.10xgenomics.com/support/software/cell-ranger/downloads), keep it
outside the project, then
`sbatch scripts/slurm/setup/build-container.sh scripts/containers/cellranger.def --build-arg CELLRANGER_TAR=/path/to/cellranger-10.1.0.tar.gz`.
The licence doesn't allow redistributing it, so `cellranger.sif` stays private: never push it to a registry.
