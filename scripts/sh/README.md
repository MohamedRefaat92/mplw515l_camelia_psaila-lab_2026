# scripts/sh/

The project's own shell scripts: downloads, format conversions, runs of command-line tools. Analysis
code goes in `src/` (R functions in `src/R/`, Python in `src/python/analysis/`).

Start each script with `set -euo pipefail` and `cd "$(dirname "$0")/../.."` (the project root), so it
runs the same from anywhere. The template's own scripts are in `tools/`: don't add yours there.
