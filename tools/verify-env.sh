#!/usr/bin/env bash
# verify-env.sh: read-only check that the environment matches its declared setup
# (environment.yml, lockfiles, kernels, editor settings, project structure) and that the layers are
# wired together correctly. It changes nothing; tools/bootstrap.sh is what repairs.
#
#   tools/verify-env.sh          fast checks
#   tools/verify-env.sh --deep   also start each kernel and run code
set -uo pipefail                       # no -e: one failed check must not end the report

deep=false
for arg in "$@"; do
  case $arg in
    --deep)    deep=true ;;
    -h|--help) sed -n '2,7p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         echo "unknown option: $arg (try --help)" >&2; exit 2 ;;
  esac
done

cd "$(dirname "$0")/.." || exit 1      # run from the project root, wherever this was called from
root=$PWD

n_pass=0 n_warn=0 n_fail=0
section() { printf '\n== %s\n' "$1"; }
info()    { printf '  · %s\n' "$1"; }
pass()    { printf '  ✔ %s\n' "$1"; n_pass=$((n_pass + 1)); }
warn()    { printf '  ! %s\n    → %s\n' "$1" "$2"; n_warn=$((n_warn + 1)); }
fail()    { printf '  ✖ %s\n    → %s\n' "$1" "$2"; n_fail=$((n_fail + 1)); }
short() { head -5 | sed 's#.*/##; s#\.conda.*##; s#\.tar\.bz2.*##' | paste -sd, - | sed 's/,/, /g'; }

section "0. Prerequisites"

for tool in mamba conda uv; do
  if command -v "$tool" >/dev/null; then
    pass "$tool found: $(command -v "$tool")"
  else
    fail "$tool not on PATH" "install Miniforge (mamba, conda) / uv, or fix PATH"
  fi
done

have_lock=false
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64)  platform=osx-arm64 ;;
  Darwin-x86_64) platform=osx-64 ;;
  Linux-x86_64)  platform=linux-64 ;;
  Linux-aarch64) platform=linux-aarch64 ;;
  *) platform="" ;;
esac
if [[ -z $platform ]]; then
  fail "unsupported platform: $(uname -s)-$(uname -m)" \
       "supported: osx-arm64, osx-64, linux-64, linux-aarch64 (on Windows use WSL)"
else
  lock="conda-${platform}.lock"
  if [[ ! -f $lock ]]; then
    warn "no lockfile for $platform" "bootstrap will solve from environment.yml"
  elif grep -q '^@EXPLICIT' "$lock"; then
    pass "lockfile: $lock"
    have_lock=true
  else
    fail "$lock has no @EXPLICIT header, so it can't recreate the env" \
         "regenerate: conda list -p ./env --explicit --md5 > $lock"
  fi
fi

have_env=false
if [[ -d env/conda-meta ]]; then
  have_env=true
  pass "./env exists"
else
  fail "no ./env (conda env)" "run tools/bootstrap.sh"
fi

section "1. Runtime (./env)"
if ! $have_env || ! $have_lock; then
  info "skipped: needs ./env and a usable lockfile"
else
  locked=$(grep '^https' "$lock" | sort)
  installed=$(conda list -p ./env --explicit --md5 | grep '^https' | sort)
  missing=$(comm -23 <(echo "$locked") <(echo "$installed"))
  extra=$(comm -13 <(echo "$locked") <(echo "$installed"))
  if [[ -z $missing ]] && [[ -z $extra ]]; then
    pass "./env matches $lock ($(grep -c '^https' "$lock") packages)"
  else
    msg=""
    [[ -n $missing ]] && msg+="missing: $(short <<<"$missing")  "
    [[ -n $extra ]]   && msg+="extra: $(short <<<"$extra")"
    fail "./env differs from $lock: $msg" \
         "if you changed the env on purpose, relock; else rm -rf env && tools/bootstrap.sh"
  fi

  names=$(sed -n '/^dependencies:/,$p' environment.yml |
    grep -E '^[[:space:]]*-[[:space:]]' |
    sed -E 's/^[[:space:]]*-[[:space:]]*//; s/[[:space:]]*#.*//; s/[=<>!~ ].*//')
  unlocked=""
  for pkg in $names; do
    if ! grep -qE "/${pkg}-[0-9]" "$lock"; then
      unlocked="${unlocked}${unlocked:+, }$pkg"
    fi
  done

  if [[ -z $unlocked ]]; then
    pass "all $(wc -w <<<"$names" | tr -d ' ') packages in environment.yml are locked"
  else
    fail "in environment.yml but not in $lock: $unlocked" \
         "mamba env update -p ./env -f environment.yml, then: conda list -p ./env --explicit --md5 > $lock"
  fi
fi

section "2. Runtime rules"

if $have_env; then
  meta=$(ls env/conda-meta)
  conda_r=$(grep '^r-' <<<"$meta" | grep -v '^r-base-[0-9]')
  if [[ -z $conda_r ]]; then
    rbase_ver=$(sed -nE 's/^r-base-([0-9][^-]*)-.*/\1/p' <<<"$meta")
    pass "only r-base $rbase_ver from conda; R packages come from renv"
  else
    names=$(sed -E 's/-[0-9].*$//' <<<"$conda_r" | paste -sd' ' -)
    fail "R packages installed by conda: $names" \
         "mamba remove -p ./env $names; drop them from environment.yml; renv::install() them instead"
  fi
else
  info "skipped: needs ./env"
fi

if [[ -z ${MAMBA_ROOT_PREFIX:-} ]]; then
  pass "MAMBA_ROOT_PREFIX unset: mamba uses its default cache (~/miniforge3/pkgs)"
elif [[ $MAMBA_ROOT_PREFIX == "$root"* ]]; then
  fail "MAMBA_ROOT_PREFIX points into the project: $MAMBA_ROOT_PREFIX" \
       "quit VS Code fully (⌘Q), reopen, check: echo \$MAMBA_ROOT_PREFIX"
else
  info "MAMBA_ROOT_PREFIX set outside the project: $MAMBA_ROOT_PREFIX"
fi

if [[ -d .micromamba ]]; then
  warn ".micromamba/ exists in the project ($(du -sh .micromamba | cut -f1))" \
       "rm -rf .micromamba once MAMBA_ROOT_PREFIX no longer points here"
fi

section "3. R (./env)"

have_r=false
r_version="" r_home="" r_lib="" r_renv_project="" r_cc="" r_synced="" renv_active=false
if $have_env; then
  # shellcheck disable=SC2016  # single-quoted R code; $ is R's list operator
  r_facts=$(PATH=/usr/bin:/bin env/bin/Rscript -e '
    cat("version=", format(getRversion()), "\n", sep = "")
    cat("home=", R.home(), "\n", sep = "")
    cat("lib=", .libPaths()[1], "\n", sep = "")
    cat("renv_project=", Sys.getenv("RENV_PROJECT"), "\n", sep = "")
    cc_name <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"), stdout = TRUE)
    cat("cc=", Sys.which(strsplit(cc_name, " ")[[1]][1]), "\n", sep = "")
    s <- tryCatch({ invisible(capture.output(x <- renv::status())); x }, error = function(e) NULL)
    cat("synced=", if (is.null(s)) "unknown" else s$synchronized, "\n", sep = "")
  ' 2>/dev/null)

  if [[ -z $r_facts ]]; then
    fail "env/bin/Rscript produced no output" "run it by hand to see the error: env/bin/Rscript -e 1"
  else
    while IFS='=' read -r key value; do
      [[ $key =~ ^[a-z_]+$ ]] || continue   # ignore stray output (e.g. renv startup warnings)
      printf -v "r_$key" '%s' "$value"     # creates r_version, r_home, ... for step 5
      info "$key: $value"
    done <<<"$r_facts"
    have_r=true
  fi
else
  info "skipped: needs ./env"
fi

if $have_r; then
  lock_r=$(env/bin/python -c 'import json; print(json.load(open("renv.lock"))["R"]["Version"])')
  if [[ $r_version == "$lock_r" ]];then
    pass "R $r_version matches renv.lock"
  else
    fail "R is $r_version but renv.lock was made with R $lock_r" \
         "if R was upgraded on purpose: reinstall packages, then renv::snapshot(); else rebuild ./env from the lockfile"
  fi
  if [[ $r_home == "$root/env/"* ]]; then
    pass "running R is ./env's R"
  else
    fail "the R that ran is not ./env's R: $r_home" \
         "env/bin/Rscript should launch ./env's R; rebuild the env (rm -rf env && tools/bootstrap.sh)"
  fi
  if [[ $r_lib == "$root/renv/library/"* ]];then
    renv_active=true
    pass "renv library active: ${r_lib#"$root/"}"
  else
    fail "R uses $r_lib, not the project's renv library" \
         "check .Rprofile: missing, guarded out, or activation failed."
  fi
  if [[ $r_renv_project == "$root" ]]; then
    pass "renv activated this project"
  elif [[ -z $r_renv_project ]]; then
    fail "renv not activated (RENV_PROJECT is empty)" \
         "check .Rprofile: missing, guarded out, or activation failed"
  else
    fail "renv activated a different project: $r_renv_project" \
         "a stale RENV_PROJECT? check this shell (echo \$RENV_PROJECT) or the kernel's env, and unset it"
  fi
  if [[ $r_synced == TRUE ]]; then
    pass "renv is synchronized"
  elif [[ $r_synced == unknown ]]; then
    info "renv sync: skipped, renv could not be loaded"
  else
    fail "renv library is out of sync with renv.lock" \
         "run env/bin/Rscript -e 'renv::status()' for details: renv::restore() if packages are missing, renv::snapshot() if the lockfile is out of date"
  fi
  if [[ $r_cc == "$root/env/"* ]];then
    pass "C compiler: ${r_cc#"$root/"}"
  elif [[ -z $r_cc ]];then
    fail "C compiler not found on R's PATH" "check the PATH block in .Rprofile; rerun tools/bootstrap.sh"
  else
    warn "C compiler outside ./env: $r_cc" "packages would build with a compiler conda didn't provide; check PATH order"
  fi
  if $renv_active; then
    libs=$(find renv/library -mindepth 2 -maxdepth 5 -name renv | sed 's#/renv$##')
    active=${r_lib#"$root/"}
    foreign=$(grep -vxF -- "$active" <<<"$libs")
    if [[ -z $foreign ]]; then
      pass "one renv library: $active"
    else
      fail "renv libraries from another R: $(paste -sd' ' - <<<"$foreign")" \
           "delete them (rm -rf <folder>) and find out which R opened this project"
    fi
  else
    info "foreign-library check: skipped, renv library not active"
  fi
fi
if [[ ! -f .Rprofile ]]; then
  fail ".Rprofile missing, so R starts without the guard, PATH setup or renv" \
       "restore it: git checkout .Rprofile"
elif [[ -z $(tail -c1 .Rprofile) ]]; then
  pass ".Rprofile ends with a newline"
else
  fail ".Rprofile does not end with a newline, so R silently skips its last line" "echo >> .Rprofile"
fi

info "R on this shell's PATH: $(command -v R || echo none)"

section "4. Python (uv)"
if ! command -v uv >/dev/null; then
  info "skipped: uv not on PATH"
else
  # 3a: uv.lock ↔ pyproject.toml
  if uv lock --check --offline >/dev/null 2>&1; then
    pass "uv.lock matches pyproject.toml"
  else
    fail "uv.lock is out of date with pyproject.toml" "uv lock, then commit uv.lock"
  fi
  # .venv exists
  have_venv=false
  if [[ -x .venv/bin/python ]]; then
    have_venv=true
    pass ".venv exists"
  else
    fail ".venv missing" "uv sync"
  fi
  if $have_venv; then
    if uv sync --check --offline >/dev/null 2>&1; then
      pass ".venv matches uv.lock"
    else
      fail ".venv is out of sync with uv.lock" "uv sync"
    fi
    want=$(<.python-version)
    have=$(.venv/bin/python -c 'import platform; print(platform.python_version())')
    if [[ $have == "$want" || $have == "$want".* ]]; then
      pass "Python $have matches .python-version ($want)"
    else
      fail "Python $have in .venv, but .python-version pins $want" "rm -rf .venv && uv sync"
    fi
  else
    info "skipped: .venv checks need .venv"
  fi
fi

section "5. Jupyter (./env)"
kernels=""
if ! $have_env; then
  info "skipped: needs ./env"
else
  # Facts from the kernelspecs and JupyterLab's config, one per line, fields separated by |
  j_facts=$(env/bin/python - 2>/dev/null <<'EOF'
import json, pathlib
for f in sorted(pathlib.Path("env/share/jupyter/kernels").glob("*/kernel.json")):
    k = json.loads(f.read_text())
    print("|".join(["kernel", f.parent.name, k.get("language", ""), k["argv"][0],
                    k.get("env", {}).get("RENV_PROJECT", "")]))
try:
    c = json.loads(pathlib.Path("env/etc/jupyter/jupyter_server_config.json").read_text())
except FileNotFoundError:
    c = {}
print("allowed|" + " ".join(c.get("KernelSpecManager", {}).get("allowed_kernelspecs", [])))
lsp = c.get("LanguageServerManager", {}).get("language_servers", {}).get("r-languageserver", {})
print("lsp|" + (lsp.get("argv") or [""])[0])
EOF
)
  if [[ -z $j_facts ]]; then
    fail "could not read the Jupyter config with env/bin/python" "run the Python block in this section by hand to see the error"
  else
    allowed="" lsp=""
    while IFS='|' read -r kind name lang prog renv; do
      case $kind in
        kernel)
          kernels+="$name "
          # 4a/4b: the program exists and belongs to the right layer (uv owns Python, conda owns R)
          if [[ ! -x $prog ]]; then
            fail "kernel $name: program missing: $prog" "tools/bootstrap.sh re-registers the kernels"
          elif [[ $lang == python && $prog == "$root/.venv/"* ]] || [[ $lang == R && $prog == "$root/env/"* ]]; then
            pass "kernel $name → ${prog#"$root/"}"
          else
            fail "kernel $name ($lang) runs $prog, not this project's $lang layer" \
                 "Python kernels must use .venv, R kernels ./env: rm -rf env/share/jupyter/kernels/$name, or tools/bootstrap.sh"
          fi
          # 4c: R kernels point renv at this project
          if [[ $lang == R && $renv != "$root" ]]; then
            fail "kernel $name points renv at ${renv:-nothing}, not this project" "folder moved or renamed? tools/bootstrap.sh"
          fi
          ;;
        allowed) allowed=$name ;;
        lsp)     lsp=$name ;;
      esac
    done <<<"$j_facts"

    if [[ -z $kernels ]]; then
      fail "no kernels registered in ./env" "tools/bootstrap.sh"
    fi

    # 4d/4e: the allowlist and the registered kernels agree
    if [[ -z $allowed ]]; then
      warn "no kernel allowlist: JupyterLab offers every kernel, including the env's own python3" "tools/bootstrap.sh writes it"
    else
      allow_ok=true
      for k in $allowed; do
        if [[ " $kernels " != *" $k "* ]]; then
          fail "allowlisted kernel $k does not exist" "tools/bootstrap.sh"
          allow_ok=false
        fi
      done
      for k in $kernels; do
        if [[ " $allowed " != *" $k "* ]]; then
          warn "kernel $k is hidden in JupyterLab (not allowlisted)" "remove it (rm -rf env/share/jupyter/kernels/$k) or allowlist it"
          allow_ok=false
        fi
      done
      if $allow_ok; then
        pass "JupyterLab allowlist matches the kernels: $allowed"
      fi
    fi

    # 4f: JupyterLab's R language server is ./env's
    if [[ -z $lsp ]]; then
      warn "no R language server configured for JupyterLab" "jupyterlab-lsp would use the first Rscript on PATH; tools/bootstrap.sh"
    elif [[ $lsp == "$root/env/"* ]]; then
      pass "JupyterLab R language server: ${lsp#"$root/"}"
    else
      fail "JupyterLab R language server runs $lsp" "tools/bootstrap.sh writes the LSP config"
    fi
  fi
fi

section "6. Editors (.vscode/settings.json)"
if [[ ! -f .vscode/settings.json ]]; then
  info "skipped: no .vscode/settings.json"
elif ! $have_env; then
  info "skipped: needs ./env (its Python reads the file)"
else
  case $platform in osx-*) os_key=mac ;; *) os_key=linux ;; esac
  # settings.json is JSONC (comments, trailing commas), so strip those before parsing.
  e_facts=$(OS_KEY=$os_key env/bin/python - 2>&1 <<'EOF'
import json, os, re
def strip_jsonc(t):
    out, i, n, in_str = [], 0, len(t), False
    while i < n:
        c = t[i]
        if in_str:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(t[i + 1]); i += 2; continue
            if c == '"':
                in_str = False
            i += 1
        elif c == '"':
            in_str = True; out.append(c); i += 1
        elif t.startswith("//", i):          # comment to end of line (not inside strings)
            j = t.find("\n", i); i = n if j < 0 else j
        elif t.startswith("/*", i):
            j = t.find("*/", i + 2); i = n if j < 0 else j + 2
        else:
            out.append(c); i += 1
    return re.sub(r",(\s*[}\]])", r"\1", "".join(out))   # trailing commas
try:
    s = json.loads(strip_jsonc(open(".vscode/settings.json").read()))
except Exception as e:
    print("parse|error|" + str(e).replace("\n", " ")); raise SystemExit
print("parse|ok|")
root, key = os.getcwd(), os.environ["OS_KEY"]
term_key = {"mac": "osx"}.get(key, key)
def get(*names):
    for name in names:
        if name in s:
            return str(s[name]).replace("${workspaceFolder}", root)
    return ""
print("rpath|" + get("r.executablePath", "r.rpath." + key) + "|")
print("rterm|" + get("r.consolePath", "r.rterm." + key) + "|")
print("termpath|" + str(s.get("terminal.integrated.env." + term_key, {}).get("PATH", "")).replace("${workspaceFolder}", root) + "|")
print("condadiscovery|" + str(s.get("positron.r.interpreters.condaDiscovery", "")).lower() + "|")
EOF
)
  while IFS='|' read -r kind value detail; do
    case $kind in
      parse)
        if [[ $value == ok ]]; then
          pass "settings.json parses (comments allowed)"
        else
          fail "settings.json does not parse, so editors silently ignore it: $detail" \
               "look for a stray character near that position (e.g. a leftover '...')"
        fi ;;
      rpath|rterm)
        what=$([[ $kind == rpath ]] && echo "R extension's R" || echo "R terminal")
        if [[ -z $value ]]; then
          warn "VS Code: $what not set" "the R extension would use the first R on PATH; set r.rpath / r.rterm"
        elif [[ $value == "$root/env/"* && -x $value ]]; then
          pass "VS Code: $what → ${value#"$root/"}"
        else
          fail "VS Code: $what is $value, not ./env's R" "point it at \${workspaceFolder}/env/bin/R"
        fi ;;
      termpath)
        if [[ $value == "$root/env/bin:"* ]]; then
          pass "integrated terminals put env/bin first on PATH"
        else
          warn "integrated terminals don't put env/bin first on PATH" "set terminal.integrated.env.$os_key PATH to \${workspaceFolder}/env/bin:\${env:PATH}"
        fi ;;
      condadiscovery)
        if [[ $value == true ]]; then
          pass "Positron: conda R discovery on, so ./env's R is listed"
        else
          info "Positron: condaDiscovery off; ./env's R won't be listed (VS Code users can ignore this)"
        fi ;;
      *)
        [[ -n $kind ]] && fail "could not read settings.json: $kind" "run the Python block in this section by hand" ;;
    esac
  done <<<"$e_facts"
fi

# Extensions that interfere with this setup (VS Code's extension folder; read-only listing).
# Uninstalled extensions keep their folder until VS Code cleans up; it lists them in .obsolete.
ext_dir=$HOME/.vscode/extensions
for ext in corker.vscode-micromamba donjayamanne.python-environment-manager; do
  installed=false
  for d in "$ext_dir/$ext"-*; do
    [[ -d $d ]] || continue
    grep -qF "\"$(basename "$d")\":true" "$ext_dir/.obsolete" 2>/dev/null || installed=true
  done
  if $installed; then
    case $ext in
      corker.*)       why="it points MAMBA_ROOT_PREFIX into the project" ;;
      donjayamanne.*) why="deprecated duplicate of ms-python.vscode-python-envs; confuses interpreter pickers" ;;
    esac
    warn "VS Code extension $ext is installed: $why" "uninstall it"
  fi
done

section "7. Project structure"
# Folders and files the template's workflow relies on
missing=""
for p in data metadata src/R src/python/analysis notebooks tools tools/R _targets.R _targets.yaml _quarto.yml; do
  [[ -e $p ]] || missing+="$p "
done
if [[ -z $missing ]]; then
  pass "folders and pipeline files present (data, metadata, src/R, src/python/analysis, notebooks, tools, tools/R, _targets.*, _quarto.yml)"
else
  fail "missing: $missing" "restore them from the template (git checkout -- <path>) or recreate them"
fi

# Everything generated stays under output/, including the targets store
if [[ -f _targets.yaml ]] && grep -qE '^[[:space:]]*store:[[:space:]]*output/_targets[[:space:]]*$' _targets.yaml; then
  pass "targets store: output/_targets"
else
  fail "_targets.yaml doesn't put the store in output/_targets" \
       "env/bin/Rscript -e 'targets::tar_config_set(store = \"output/_targets\")'"
fi

# data/ holds only symlinks to inputs that exist
broken="" copies="" n_links=0
for f in data/* data/.[!.]*; do
  [[ -e $f || -L $f ]] || continue
  [[ $f == data/README.md ]] && continue
  if [[ -L $f ]]; then
    n_links=$((n_links + 1))
    [[ -e $f ]] || broken+="${f#data/} "
  else
    copies+="${f#data/} "
  fi
done
if [[ -n $broken ]]; then
  fail "broken input links in data/: $broken" "the data moved or this is another machine: tools/link-input.sh <path> --name <name>"
elif [[ -n $copies ]]; then
  warn "real files in data/ (not symlinks): $copies" "keep data outside the project; link it with tools/link-input.sh"
else
  pass "data/: $n_links input link(s), all resolve"
fi

# The reviewed data catalogue lists every file in data/, and every file it lists is still there
reviewed=metadata/catalogue_data.tsv
if [[ ! -f $reviewed ]]; then
  warn "no reviewed data catalogue ($reviewed)" \
       "tools/catalogue-data.sh, review metadata/auto_catalogue_data.tsv, copy it to $reviewed, commit"
else
  listed=$(tail -n +2 "$reviewed" | cut -f1 | sort)
  on_disk=$(cd data && find -L . -type f ! -name README.md 2>/dev/null | sed 's#^\./##' | sort)
  uncatalogued=$(comm -13 <(echo "$listed") <(echo "$on_disk") | grep -v '^$' || true)
  gone=$(comm -23 <(echo "$listed") <(echo "$on_disk") | grep -v '^$' || true)
  if [[ -n $uncatalogued ]]; then
    warn "$(wc -l <<<"$uncatalogued" | tr -d ' ') file(s) in data/ not in the reviewed catalogue: $(head -3 <<<"$uncatalogued" | paste -sd' ' -)" \
         "tools/catalogue-data.sh, review, copy to $reviewed"
  elif [[ -n $gone ]]; then
    warn "$(wc -l <<<"$gone" | tr -d ' ') catalogued file(s) not in data/: $(head -3 <<<"$gone" | paste -sd' ' -)" \
         "relink them (tools/link-input.sh), or remove their rows if the project no longer uses them"
  else
    pass "reviewed data catalogue lists all $(grep -c . <<<"$listed") file(s) in data/"
  fi
fi

# Generated files and machine-specific links never belong in git
tracked=$( { git ls-files -- output data metadata/auto_* 2>/dev/null || true; } | grep -vxE '(data|output)/README\.md|(output/reports/(latest|archive)|output/objects)/README\.md' || true)
if [[ -z $tracked ]]; then
  pass "nothing under output/, data/ or metadata/auto_* is tracked by git"
else
  fail "tracked by git but generated or machine-specific: $(paste -sd' ' - <<<"$tracked")" \
       "git rm --cached <path> (keeps the file, stops tracking it)"
fi

# Stray files in the project root (e.g. from a mistyped redirect like `... > bash`)
expected=" .git .gitignore .python-version .Rprofile .Rhistory .renvignore .venv .vscode .quarto .DS_Store
  .ipynb_checkpoints __marimo__ __pycache__ README.md LICENSE _dependencies.R _targets.R _targets.yaml
  _quarto.yml environment.yml pyproject.toml uv.lock renv renv.lock env data metadata src notebooks output scripts tools "
strays=""
for f in * .[!.]*; do
  [[ -e $f || -L $f ]] || continue
  [[ $f == conda-*.lock ]] && continue
  [[ $expected == *[[:space:]]"$f"[[:space:]]* ]] || strays+="$f "
done
if [[ -z $strays ]]; then
  pass "no unexpected files in the project root"
else
  warn "unexpected in the project root: $strays" "move it into the structure, delete it, or add it to 'expected' in this script"
fi

if $deep; then
  section "8. Kernels, live (--deep)"
  if [[ -z $kernels ]]; then
    info "skipped: no kernels found in section 5"
  else
    sub=src; [[ -d $sub ]] || sub=.
    # Start each kernel from a subfolder and ask it which R library / Python it really uses.
    d_facts=$(cd "$sub" && "$root/env/bin/python" - "$kernels" 2>/dev/null <<'EOF'
import sys
import jupyter_client
from jupyter_client.kernelspec import KernelSpecManager
code = {"R": "cat(.libPaths()[1])", "python": "import sys; print(sys.executable, end='')"}
ksm = KernelSpecManager()
for name in sys.argv[1].split():
    try:
        lang = ksm.get_kernel_spec(name).language
    except Exception:
        print(f"{name}|error||no kernelspec found"); continue
    if lang not in code:
        print(f"{name}|error||no test for language {lang}"); continue
    try:
        km, kc = jupyter_client.manager.start_new_kernel(kernel_name=name, startup_timeout=60)
    except Exception as e:
        print(f"{name}|error||did not start: {e}".replace("\n", " ")); continue
    out = []
    def hook(msg):
        if msg["msg_type"] == "stream":
            out.append(msg["content"]["text"])
    try:
        status = kc.execute_interactive(code[lang], timeout=60, output_hook=hook)["content"]["status"]
    except Exception as e:
        status = f"no reply ({type(e).__name__})"
    finally:
        kc.stop_channels()
        km.shutdown_kernel(now=True)
    print(f"{name}|{lang}|{status}|{''.join(out).strip()}")
EOF
)
    if [[ -z $d_facts ]]; then
      fail "could not run the live kernel test" "run the Python block in this section by hand to see the error"
    fi
    while IFS='|' read -r name lang status value; do
      [[ -n $name ]] || continue
      if [[ $lang == error ]]; then
        fail "kernel $name: $value" "tools/bootstrap.sh; open it in JupyterLab to see the full error"
      elif [[ $status != ok ]]; then
        fail "kernel $name started, but the test code failed ($status)" "open it in JupyterLab and run the code by hand"
      elif [[ $lang == R && $value == "$root/renv/library/"* ]] || [[ $lang == python && $value == "$root/.venv/"* ]]; then
        pass "kernel $name, started from $sub/: ${value#"$root/"}"
      else
        fail "kernel $name, started from $sub/, uses ${value:-nothing}" \
             "R kernels must load renv, Python kernels .venv: check kernel.json and .Rprofile"
      fi
    done <<<"$d_facts"
  fi
fi

printf '\n%d passed, %d warnings, %d failed\n' "$n_pass" "$n_warn" "$n_fail"
[[ $n_fail -eq 0 ]]                    # exit code 1 if anything failed

