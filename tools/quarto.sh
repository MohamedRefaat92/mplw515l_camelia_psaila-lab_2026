#!/usr/bin/env bash
# Run ./env's Quarto from a terminal. conda's Quarto needs the variables its activation script sets
# (QUARTO_SHARE_PATH, QUARTO_DENO, ...), and this project never activates ./env; R gets them from .Rprofile.
# Reports normally render through the pipeline (tools/run-pipeline.sh); use this for ad-hoc renders.
#
#   tools/quarto.sh render notebooks/<report>/<report>.qmd
#   tools/quarto.sh --version
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)

activate="$root/env/etc/conda/activate.d/quarto.sh"
if [[ ! -f $activate ]]; then
  echo "Quarto is not installed in ./env (run tools/bootstrap.sh)" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$activate"
export PATH="$root/env/bin:$PATH"
exec "$root/env/bin/quarto" "$@"
