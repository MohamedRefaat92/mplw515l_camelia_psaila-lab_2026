"""
cellrank_run.py — CellRank trajectory analysis.

Input type and velocity mode determine which kernel is used and how
outputs are named — allowing multiple runs to coexist.

Env vars
--------
INPUT_TYPE       cytotrace | velocity | dpt    (default: velocity)
VELOCITY_MODE    stochastic | dynamical        (default: dynamical)
                   only used when INPUT_TYPE=velocity
H5AD_PATH        Path to input h5ad
                   cytotrace  → hsc.h5ad (from targets pipeline)
                   velocity   → hsc_velocity_{mode}.h5ad (from velocity_scvelo.py)
                   dpt        → output/objects/dpt/hsc_dpt.h5ad (from dpt_hsc.ipynb)
OUT_DIR          Root output dir (default: output/objects/cellrank)
N_MACROSTATES    Number of CFLARE eigenvectors  (default: 6)
N_JOBS           Parallel jobs                  (default: 8)

Run tag (used in all output names)
-----------------------------------
  cytotrace                → run_tag = "cytotrace"
  dpt                      → run_tag = "dpt"
  velocity + stochastic    → run_tag = "velocity_stochastic"
  velocity + dynamical     → run_tag = "velocity_dynamical"

Output layout
-------------
<OUT_DIR>/<run_tag>/
  checkpoints/
    adata_prepped.h5ad         — after neighbours + kernel computation
  figures/
    overview_umap.png
    kernel_projection.png
    cflare_spectrum.png
    terminal_states.png
    initial_states.png
    fate_probabilities.png
    fate_by_condition.png
    lineage_drivers/
      <state>.png
    gene_trends.png
  tables/
    transition_matrix.npz      — full cell×cell sparse transition matrix
    macrostate_memberships.csv — cell × macrostate soft assignment
    terminal_states.csv        — per-cell terminal state hard assignment
    initial_states.csv         — per-cell initial state hard assignment
    fate_probabilities.csv     — cell × terminal_state absorption probability
    fate_prob_per_cluster.csv  — cluster-mean fate probabilities
    lineage_drivers_<state>.csv — per-terminal-state ranked driver genes
  hsc_cellrank_{run_tag}.h5ad — final object with fate probs in obs
"""
import os, warnings
warnings.filterwarnings('ignore')

import anndata as ad
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import sparse

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_TYPE    = os.environ.get('INPUT_TYPE',    'velocity').lower()
VELOCITY_MODE = os.environ.get('VELOCITY_MODE', 'dynamical').lower()
CLUSTER_RES   = os.environ.get('CLUSTER_RES',   None)   # e.g. "1.5"; None → auto-pick
N_MACROSTATES = int(os.environ.get('N_MACROSTATES', 6))
N_JOBS        = int(os.environ.get('N_JOBS', 8))
CONDITION_COL = 'sample.name'
CYTOTRACE_COL = 'cytotrace2_score'

assert INPUT_TYPE in ('cytotrace', 'velocity', 'dpt'), \
    f"INPUT_TYPE must be 'cytotrace', 'velocity', or 'dpt', got {INPUT_TYPE!r}"
assert VELOCITY_MODE in ('stochastic', 'dynamical'), \
    f"VELOCITY_MODE must be 'stochastic' or 'dynamical', got {VELOCITY_MODE!r}"

# Run tag drives all output naming; includes resolution when explicitly set
_res_slug = CLUSTER_RES.replace('.', '_') if CLUSTER_RES else None
RUN_TAG = (
    ('cytotrace'    if INPUT_TYPE == 'cytotrace'
     else 'dpt'     if INPUT_TYPE == 'dpt'
     else f'velocity_{VELOCITY_MODE}')
    + (f'_res{_res_slug}' if _res_slug else '')
)

# Default h5ad path inferred from run type if not provided
_default_h5ad = (
    'output/objects/hsc.h5ad'                     if INPUT_TYPE == 'cytotrace'
    else 'output/objects/dpt/hsc_dpt.h5ad'              if INPUT_TYPE == 'dpt'
    else f'output/objects/velocity/hsc_velocity_{VELOCITY_MODE}.h5ad'
)
H5AD_PATH = Path(os.environ.get('H5AD_PATH', _default_h5ad))

_base_out   = Path(os.environ.get('OUT_DIR', 'output/objects/cellrank'))
OUT_DIR     = _base_out / RUN_TAG

CHKPT_DIR   = OUT_DIR / 'checkpoints'
FIG_DIR     = OUT_DIR / 'figures'
TABLE_DIR   = OUT_DIR / 'tables'
DRIVERS_DIR = FIG_DIR / 'lineage_drivers'
for d in (CHKPT_DIR, FIG_DIR, TABLE_DIR, DRIVERS_DIR):
    d.mkdir(parents=True, exist_ok=True)

GENES_OF_INTEREST = ['Gata1', 'Pf4', 'Vwf', 'Csf1r', 'Ly6c2',
                     'Cd34', 'Kit', 'Slamf1', 'Mpl']

print(f"cellrank    {cr.__version__}")
print(f"Input type: {INPUT_TYPE}")
print(f"Velocity mode: {VELOCITY_MODE if INPUT_TYPE == 'velocity' else 'N/A'}")
print(f"Run tag:    {RUN_TAG}")
print(f"H5AD:       {H5AD_PATH}")
print(f"Output:     {OUT_DIR}")

def save_fig(name, dpi=150):
    path = FIG_DIR / name
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close('all')
    print(f"  → {path}")

def _fix_nullable_strings(df):
    import pandas as pd
    df.index = pd.Index(df.index.astype(str).values.tolist(), dtype=object)
    for col in df.columns:
        dtype_name = str(df[col].dtype)
        if (dtype_name in ('string', 'string[pyarrow]', 'StringDtype')
                or 'StringDtype' in type(df[col].dtype).__name__
                or dtype_name.lower().startswith('string')):
            df[col] = df[col].astype(str).values.astype(object)
    return df

# ── 1. Load h5ad ──────────────────────────────────────────────────────────────
PREPPED_CHKPT = CHKPT_DIR / 'adata_prepped.h5ad'

if PREPPED_CHKPT.exists():
    print(f"\n=== 1. Load prepped checkpoint ===")
    adata = sc.read_h5ad(PREPPED_CHKPT)
    print(f"  {adata.shape}")
else:
    print(f"\n=== 1. Load h5ad ===")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"  {adata.shape}")

    # Cluster column — use CLUSTER_RES if set, else prefer res.0.9, else any res.
    if CLUSTER_RES:
        _res_pattern = f'res.{CLUSTER_RES}'
        CLUSTER_COL = next(
            (c for c in adata.obs.columns if _res_pattern in c),
            next((c for c in adata.obs.columns if 'res.' in c), 'seurat_clusters')
        )
    else:
        CLUSTER_COL = next(
            (c for c in adata.obs.columns if 'res.0.9' in c),
            next((c for c in adata.obs.columns if 'res.' in c), 'seurat_clusters')
        )
    adata.uns['_cellrank_cluster_col'] = CLUSTER_COL
    adata.obs[CLUSTER_COL] = adata.obs[CLUSTER_COL].astype(str).astype('category')
    print(f"  Cluster column: {CLUSTER_COL!r}")

    # Pseudotime column
    if INPUT_TYPE == 'cytotrace':
        assert CYTOTRACE_COL in adata.obs.columns, \
            f"'{CYTOTRACE_COL}' missing — rebuild hsc targets first."
        adata.obs['cytotrace2_pseudotime'] = 1.0 - adata.obs[CYTOTRACE_COL]
        print(f"  CytoTRACE2 pseudotime range: "
              f"{adata.obs['cytotrace2_pseudotime'].min():.3f} – "
              f"{adata.obs['cytotrace2_pseudotime'].max():.3f}")
    elif INPUT_TYPE == 'dpt':
        assert 'dpt_pseudotime' in adata.obs.columns, \
            "dpt_pseudotime missing — run dpt_hsc.ipynb first."
        print(f"  DPT pseudotime range: "
              f"{adata.obs['dpt_pseudotime'].min():.3f} – "
              f"{adata.obs['dpt_pseudotime'].max():.3f}")

    # Neighbours
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=30, n_pcs=30)
        print("  Computed neighbours from X_pca")
    else:
        print("  Neighbours already present")

    # Convert nullable / arrow string dtypes → numpy object so anndata can write
    adata.obs = _fix_nullable_strings(adata.obs)
    adata.var = _fix_nullable_strings(adata.var)
    adata.write_h5ad(PREPPED_CHKPT)
    print(f"  Checkpoint: {PREPPED_CHKPT}")

if CLUSTER_RES:
    _res_pattern = f'res.{CLUSTER_RES}'
    CLUSTER_COL = adata.uns.get('_cellrank_cluster_col', next(
        (c for c in adata.obs.columns if _res_pattern in c),
        next((c for c in adata.obs.columns if 'res.' in c), 'seurat_clusters')
    ))
else:
    CLUSTER_COL = adata.uns.get('_cellrank_cluster_col', next(
        (c for c in adata.obs.columns if 'res.0.9' in c),
        next((c for c in adata.obs.columns if 'res.' in c), 'seurat_clusters')
    ))
print(f"Cluster column: {CLUSTER_COL!r}")

# ── 2. Overview UMAPs ─────────────────────────────────────────────────────────
print(f"\n=== 2. Overview UMAPs ===")
panels = [CLUSTER_COL]
if INPUT_TYPE == 'cytotrace':
    panels.append('cytotrace2_pseudotime')
elif INPUT_TYPE == 'dpt':
    panels.append('dpt_pseudotime')
elif 'latent_time' in adata.obs.columns:
    panels.append('latent_time')
if CONDITION_COL in adata.obs.columns:
    panels.append(CONDITION_COL)

fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 5))
axes = [axes] if len(panels) == 1 else list(axes)
cmaps = {'cytotrace2_pseudotime': 'RdYlBu_r', 'dpt_pseudotime': 'viridis',
         'latent_time': 'gnuplot'}
for ax, col in zip(axes, panels):
    sc.pl.umap(adata, color=col, ax=ax, show=False,
               color_map=cmaps.get(col, None), title=col)
plt.tight_layout()
save_fig('overview_umap.png')

# ── 3. Build kernel ───────────────────────────────────────────────────────────
print(f"\n=== 3. Build kernel ({RUN_TAG}) ===")
from cellrank.kernels import VelocityKernel, PseudotimeKernel, ConnectivityKernel

ck = ConnectivityKernel(adata).compute_transition_matrix()

if INPUT_TYPE == 'velocity':
    vk = VelocityKernel(adata).compute_transition_matrix()
    combined = 0.8 * vk + 0.2 * ck
    print("  0.8 × VelocityKernel + 0.2 × ConnectivityKernel")
else:
    _pt_key = 'dpt_pseudotime' if INPUT_TYPE == 'dpt' else 'cytotrace2_pseudotime'
    _pt_label = 'DPT' if INPUT_TYPE == 'dpt' else 'CytoTRACE2'
    ptk = PseudotimeKernel(adata, time_key=_pt_key)
    ptk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)
    combined = 0.8 * ptk + 0.2 * ck
    print(f"  0.8 × PseudotimeKernel ({_pt_label}) + 0.2 × ConnectivityKernel")

print(combined)

# ── 4. Save transition matrix ─────────────────────────────────────────────────
print(f"\n=== 4. Save transition matrix ===")
T = combined.transition_matrix          # scipy sparse (csr or csr_array)
if not sparse.issparse(T):
    T = sparse.csr_matrix(T)
tm_path = TABLE_DIR / 'transition_matrix.npz'
sparse.save_npz(tm_path, T)
print(f"  Saved: {tm_path}  shape={T.shape}  nnz={T.nnz}")

# Save cell index so matrix rows/cols can be mapped back to obs_names
pd.Series(adata.obs_names, name='cell_id').to_csv(
    TABLE_DIR / 'transition_matrix_cells.csv', index=False)

# ── 4b. Cluster-level transition matrix ───────────────────────────────────────
# Aggregate T (n_cells × n_cells) to (n_clusters × n_clusters) without
# densifying.  Uses left-multiply by row-normalised cluster indicator L and
# right-multiply by cluster indicator R:
#   cluster_T = L @ T @ R
# where L[i, c] = 1/n_i  (mean over cells in cluster i)
#       R[c, j] = 1       (sum into cluster j)
# Result: cluster_T[i, j] = mean probability mass flowing from cluster i to j.
_cl_labels  = adata.obs[CLUSTER_COL].astype(str).values
_cl_unique  = sorted(set(_cl_labels))
_n_cl       = len(_cl_unique)
_n_cells    = T.shape[0]
_cl_to_idx  = {c: i for i, c in enumerate(_cl_unique)}

_rows_L, _cols_L, _data_L = [], [], []
_rows_R, _cols_R, _data_R = [], [], []
for _c, _ci in _cl_to_idx.items():
    _cell_idx = np.where(_cl_labels == _c)[0]
    _n_i      = len(_cell_idx)
    for _cell in _cell_idx:
        _rows_L.append(_ci);  _cols_L.append(_cell); _data_L.append(1.0 / _n_i)
        _rows_R.append(_cell); _cols_R.append(_ci);  _data_R.append(1.0)

_L = sparse.csr_matrix((_data_L, (_rows_L, _cols_L)), shape=(_n_cl, _n_cells))
_R = sparse.csr_matrix((_data_R, (_rows_R, _cols_R)), shape=(_n_cells, _n_cl))
_cl_T = (_L @ T @ _R).toarray()

cl_T_df = pd.DataFrame(_cl_T, index=_cl_unique, columns=_cl_unique)
cl_T_df.index.name = 'cluster'
cl_T_df.to_csv(TABLE_DIR / 'cluster_transition_matrix.csv')
print(f"  cluster_transition_matrix.csv  {cl_T_df.shape}")
print(cl_T_df.round(4))

# Kernel projection
ax = combined.plot_projection(basis='umap', color=CLUSTER_COL, show=False)
save_fig('kernel_projection.png')

# ── 5. CFLARE ─────────────────────────────────────────────────────────────────
print(f"\n=== 5. CFLARE (k={N_MACROSTATES}) ===")
from cellrank.estimators import CFLARE

c = CFLARE(combined)
c.fit(k=N_MACROSTATES)

c.plot_spectrum()
save_fig('cflare_spectrum.png')

c.predict(use=N_MACROSTATES, cluster_key=CLUSTER_COL)

# Terminal states
c.plot_macrostates(which='terminal', basis='umap', legend_loc='right margin')
save_fig('terminal_states.png')

# Initial states
try:
    c.plot_macrostates(which='initial', basis='umap', legend_loc='right margin')
    save_fig('initial_states.png')
except Exception as e:
    print(f"  [SKIP] initial states plot: {e}")

# ── 6. Save macrostate tables ─────────────────────────────────────────────────
print(f"\n=== 6. Save macrostate tables ===")

# Soft macrostate memberships (GPCCA only; CFLARE does not have this)
_mem = getattr(c, 'macrostates_memberships', None)
if _mem is not None:
    mem_df = pd.DataFrame(
        _mem.X if hasattr(_mem, 'X') else np.array(_mem),
        columns=_mem.names, index=adata.obs_names
    )
    mem_df.to_csv(TABLE_DIR / 'macrostate_memberships.csv')
    print(f"  macrostate_memberships.csv  {mem_df.shape}")
else:
    print("  [SKIP] macrostate_memberships: not available for CFLARE")

# Terminal state hard assignment per cell
if c.terminal_states is not None:
    ts_df = pd.DataFrame({
        'terminal_state': c.terminal_states,
        CLUSTER_COL:      adata.obs[CLUSTER_COL].values,
    }, index=adata.obs_names)
    if CONDITION_COL in adata.obs.columns:
        ts_df[CONDITION_COL] = adata.obs[CONDITION_COL].values
    ts_df.to_csv(TABLE_DIR / 'terminal_states.csv')
    print(f"  terminal_states.csv  {ts_df.shape}")

# Initial state assignment
try:
    if c.initial_states is not None:
        is_df = pd.DataFrame({'initial_state': c.initial_states},
                              index=adata.obs_names)
        is_df.to_csv(TABLE_DIR / 'initial_states.csv')
        print(f"  initial_states.csv  {is_df.shape}")
except Exception as e:
    print(f"  [SKIP] initial_states.csv: {e}")

# ── 7. Fate probabilities ─────────────────────────────────────────────────────
print(f"\n=== 7. Fate probabilities ===")
c.compute_fate_probabilities(solver='gmres', tol=1e-7, preconditioner='ilu')
c.plot_fate_probabilities(basis='umap', same_plot=False, ncols=3)
save_fig('fate_probabilities.png')

fp_arr = c.fate_probabilities
fp_df  = pd.DataFrame(
    fp_arr.X if hasattr(fp_arr, 'X') else np.array(fp_arr),
    columns=fp_arr.names, index=adata.obs_names
)
fp_df[CLUSTER_COL] = adata.obs[CLUSTER_COL].values
if CONDITION_COL in adata.obs.columns:
    fp_df[CONDITION_COL] = adata.obs[CONDITION_COL].values

fp_df.to_csv(TABLE_DIR / 'fate_probabilities.csv')
print(f"  fate_probabilities.csv  {fp_df.shape}")

# Per-cluster mean fate probability
cluster_fp = fp_df.drop(columns=[CONDITION_COL], errors='ignore') \
                   .groupby(CLUSTER_COL).mean(numeric_only=True)
cluster_fp.to_csv(TABLE_DIR / 'fate_prob_per_cluster.csv')
print(f"  fate_prob_per_cluster.csv  {cluster_fp.shape}")
print(cluster_fp.round(3))

# Fate by condition box plot
fp_cols = [col for col in fp_df.columns
           if col not in [CLUSTER_COL, CONDITION_COL]]
if fp_cols and CONDITION_COL in fp_df.columns:
    n = len(fp_cols)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4))
    axes = [axes] if n == 1 else list(axes)
    for ax, col in zip(axes, fp_cols):
        fp_df.boxplot(column=col, by=CONDITION_COL, ax=ax, grid=False, rot=30)
        ax.set_title(col, fontsize=9)
        ax.set_xlabel('')
    plt.suptitle(f'Fate probability by condition ({RUN_TAG})', y=1.02)
    plt.tight_layout()
    save_fig('fate_by_condition.png')

# ── 8. Lineage drivers ────────────────────────────────────────────────────────
print(f"\n=== 8. Lineage drivers ===")
c.compute_lineage_drivers(cluster_key=CLUSTER_COL, return_drivers=True)
term_names = list(c.terminal_states.cat.categories)
print(f"  Terminal states: {term_names}")

for ts in term_names:
    # Save driver gene table
    try:
        drivers = c.lineage_drivers
        ts_clean = ts.replace(' ', '_').replace('/', '-')
        if drivers is not None and ts in drivers.columns.get_level_values(0):
            drv_df = drivers[ts].sort_values('corr', ascending=False)
            drv_df.to_csv(TABLE_DIR / f'lineage_drivers_{ts_clean}.csv')
            print(f"  lineage_drivers_{ts_clean}.csv  {drv_df.shape}")
    except Exception as e:
        print(f"  [SKIP] driver table {ts}: {e}")

    # Save driver plot
    try:
        c.plot_lineage_drivers(lineage=ts, n_genes=8)
        save_fig(f'lineage_drivers/{ts_clean}.png')
    except Exception as e:
        plt.close('all')
        print(f"  [SKIP] driver plot {ts}: {e}")

# ── 9. Gene trends ────────────────────────────────────────────────────────────
print(f"\n=== 9. Gene trends ===")
time_key = ('latent_time'         if 'latent_time'            in adata.obs.columns
            else 'dpt_pseudotime'  if 'dpt_pseudotime'         in adata.obs.columns
            else 'cytotrace2_pseudotime')
genes_ok = [g for g in GENES_OF_INTEREST if g in adata.var_names]
print(f"  Time key: {time_key!r}   Genes: {genes_ok}")
if genes_ok:
    model = cr.models.GAM(adata)
    cr.pl.gene_trends(
        adata, model=model, genes=genes_ok,
        time_key=time_key, data_key='X', ncols=3, same_plot=False,
    )
    save_fig('gene_trends.png')

# ── 10. Save final h5ad ───────────────────────────────────────────────────────
print(f"\n=== 10. Save ===")
# Attach fate probabilities to obs
for col in fp_cols:
    adata.obs[f'fate_{col}'] = fp_df[col].values
if c.terminal_states is not None:
    adata.obs['terminal_state'] = c.terminal_states

out_h5ad = OUT_DIR / f'hsc_cellrank_{RUN_TAG}.h5ad'
adata.obs = _fix_nullable_strings(adata.obs)
adata.var = _fix_nullable_strings(adata.var)
adata.write_h5ad(out_h5ad)
print(f"Final h5ad: {out_h5ad}")
print(f"\n=== Done — {RUN_TAG} ===")
