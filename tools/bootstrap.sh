#!/usr/bin/env bash
# Build (or repair) all project environments. Safe to re-run.
#
#   tools/bootstrap.sh            exact rebuild from conda-<platform>.lock when one exists
#   tools/bootstrap.sh --resolve  re-solve the conda env from environment.yml instead
#
# Needs: Miniforge (mamba + conda) and uv on PATH. macOS/Linux; on Windows use WSL.
set -euo pipefail
cd "$(dirname "$0")/.."   # run from the project root, wherever this was called from
# Kernel names allow only letters, digits, '.', '_' and '-'
name=$(basename "$PWD" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9._-]/-/g')

# Keep mamba's package cache in ~/miniforge3, not in the project
unset MAMBA_ROOT_PREFIX MAMBA_EXE

resolve=false
[[ "${1:-}" == "--resolve" ]] && resolve=true

case "$(uname -s)-$(uname -m)" in
  Darwin-arm64)  platform=osx-arm64 ;;
  Darwin-x86_64) platform=osx-64 ;;
  Linux-x86_64)  platform=linux-64 ;;
  Linux-aarch64) platform=linux-aarch64 ;;
  *) echo "Unsupported platform: $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac
lock="conda-${platform}.lock"

echo "==> [1/5] conda env (./env): runtime, C libraries, CLI tools, JupyterLab"
if [[ -d env/conda-meta ]]; then
  echo "    ./env exists, skipping (delete it to rebuild)"
elif [[ -f "$lock" && "$resolve" == false ]]; then
  mamba create -y -p ./env --file "$lock"
else
  mamba env create -y -p ./env -f environment.yml
  echo "    Solved fresh. Record it with: conda list -p ./env --explicit --md5 > $lock"
fi

echo "==> [2/5] JupyterLab config: project kernels only, R language server from ./env"
# jupyterlab-lsp would otherwise run whichever Rscript comes first on PATH
PROJECT_NAME="$name" env/bin/python - <<'EOF'
import json, os, pathlib
name = os.environ["PROJECT_NAME"]
root = pathlib.Path.cwd()
config = {
    "KernelSpecManager": {"allowed_kernelspecs": [f"{name}-r", f"{name}-py"]},
    "LanguageServerManager": {"language_servers": {"r-languageserver": {
        "argv": [str(root / "env/bin/Rscript"), "--slave", "-e", "languageserver::run()"],
        "display_name": f"languageserver ({name})",
        "languages": ["r"],
        "mime_types": ["text/x-rsrc"],
        "version": 2,
        "env": {"RENV_PROJECT": str(root), "R_PROFILE_USER": str(root / ".Rprofile")},
    }}},
}
path = root / "env/etc/jupyter/jupyter_server_config.json"
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(config, indent=2) + "\n")
EOF

echo "==> [3/5] Python packages (uv)"
uv sync --locked

echo "==> [4/5] R packages (renv)"
env/bin/Rscript -e 'renv::restore(prompt = FALSE)'

echo "==> [5/5] Jupyter kernels"
rm -rf env/share/jupyter/kernels/*   # packages (e.g. ipykernel) ship their own kernelspecs; keep only ours
env/bin/Rscript tools/R/register-r-kernel.R "$name"
# ipykernel warns that ./env "may not be found": that's the .venv's view; env/bin/jupyter finds it
uv run --locked python -m ipykernel install --prefix ./env \
  --name "$name-py" --display-name "Python ($name · uv)"

echo "Done. Start JupyterLab with: env/bin/jupyter lab"
