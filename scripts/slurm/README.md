# scripts/slurm/

Cluster job files (`#SBATCH` headers) for work that runs on the cluster: Cell Ranger, alignment,
or the pipeline itself (`tools/run-pipeline.sh` inside a job).

Keeping job files with the project records how the inputs in `data/` were produced.
Rules: write logs to `output/logs/`, never into `data/`; resource requests (time, memory, cores)
go in the `#SBATCH` header, not in the command line, so the file alone documents the run.

## In this project

Submit from the project root (`sbatch scripts/slurm/<folder>/<job>.sh`); each job first sources
`env.sh` (project-root check, `output/logs/`, and `in_container`, which runs a command in one of the
images in `scripts/containers/` with the project and the folders behind `data/` visible). How to run
each job, and with which settings: its header.

| Folder | Jobs | Reads → writes |
| --- | --- | --- |
| `cellranger/` | `cellranger-multi.sh` (one array task per tube; paths in `cellranger/env.sh`) | the FASTQ tar archives → `cellranger_multi/Rx_<n>/` in the data hub |
| `setup/` | `install-r-packages.sh`, `build-container.sh <recipe.def>` | → `renv/library` + `renv.lock`; → `output/containers/<name>.sif` |
| `pipeline/` | `run-pipeline.sh` | the targets pipeline and its reports |
| `velocity/` | `velocyto.sh` (per tube), then `scvelo.sh` | `data/cellranger_Rx_<n>/` BAMs, `hsc.h5ad` → `output/objects/velocyto/`, `velocity/` |
| `analysis/` | `cellrank.sh`, `paga.sh`, `scenic.sh` (per subset), `perturbgen.sh` | `output/objects/hsc.h5ad` → `output/objects/<tool>/` |

The `analysis/` and `velocity/` jobs start from the pipeline's `output/objects/hsc.h5ad`; the reports
that read their results render on the next pipeline run.
