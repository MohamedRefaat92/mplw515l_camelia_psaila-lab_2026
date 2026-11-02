# scripts/slurm/

Cluster job files (`#SBATCH` headers) for work that runs on the cluster: Cell Ranger, alignment,
or the pipeline itself (`tools/run-pipeline.sh` inside a job).

Keeping job files with the project records how the inputs in `data/` were produced.
Rules: write logs to `output/logs/`, never into `data/`; resource requests (time, memory, cores)
go in the `#SBATCH` header, not in the command line, so the file alone documents the run.
