# scripts/sh/

The project's own shell scripts: downloads, format conversions, runs of command-line tools. Analysis
code goes in `src/` (R functions in `src/R/`, Python in `src/python/analysis/`).

Start each script with `set -euo pipefail` and `cd "$(dirname "$0")/../.."` (the project root), so it
runs the same from anywhere. The template's own scripts are in `tools/`: don't add yours there.

## In this project

| Script | Does |
| --- | --- |
| `link-references.sh [hub]` | links the reference datasets from the data hub into `data/`; copies the sample sheet into `metadata/samples.csv` |
| `link-cellranger.sh` | links each tube's Cell Ranger results as `data/cellranger_Rx_<n>` |
| `get-scenic-resources.sh <folder>` | downloads the mm10 cisTarget databases if missing, links them as `data/scenic` |
| `get-10x-reference.sh <name>` | downloads and unpacks a 10x prebuilt reference into `data_hub/10x_resources/reference_genomes/raw/` |
