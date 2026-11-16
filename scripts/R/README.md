# scripts/R/

The project's R scripts that aren't analysis code (the analysis functions are in `src/R/`). Run them
from the project root with `env/bin/Rscript scripts/R/<script>.R`.

| Script | Does |
| --- | --- |
| `install-r-packages.R` (+ `Makevars`) | installs the analysis's R packages into renv and writes `renv.lock`; run as `sbatch scripts/slurm/setup/install-r-packages.sh` |
| `rebuild-refannot-mnc-stroma.R` | rebuilds just the MNC and Stroma reference-annotation targets |
| `build-pansci-azimuth.R`, `test-pansci-load.R` | retired with the PanSci references (their targets are commented out in `_targets.R`) |
