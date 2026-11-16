"""
velocity_scvelo.py — scVelo RNA velocity from velocyto loom files.

Env vars
--------
H5AD_PATH       Path to input hsc.h5ad (with X_pca, X_umap)
LOOM_DIR        Dir containing per-sample loom subdirs (Rx_*/<sample>.loom)
VELOCITY_MODE   stochastic | dynamical  (default: dynamical)
OUT_DIR         Root output dir         (default: output/objects/velocity)
N_JOBS          Parallel jobs           (default: 8)

Resume logic
------------
Checkpoint h5ads are written after each expensive step.  Re-running the
script skips completed steps automatically:
  checkpoints/merged.h5ad               — after loom merge & barcode alignment
  checkpoints/moments_{mode}.h5ad       — after moments
  checkpoints/dynamics_{mode}.h5ad      — after recover_dynamics (dynamical only)

Output layout
-------------
<OUT_DIR>/
  checkpoints/
    merged.h5ad
    moments_{mode}.h5ad
    dynamics_{mode}.h5ad          (dynamical only)
  figures/
    velocity_stream_{mode}.png
    latent_time_{mode}.png        (dynamical only)
    phase_portraits_{mode}.png
    proportions_{mode}.png
  tables/
    velocity_genes_{mode}.csv     — per-gene kinetic params + likelihood
    latent_time_{mode}.csv        — per-cell latent time (dynamical only)
  hsc_velocity_{mode}.h5ad       — final merged object
"""
import os, warnings
# Must be set before numpy/scipy are imported — limits OpenBLAS/MKL threads so
# that fork() in multiprocessing children doesn't inherit a multi-threaded state
# and deadlock.  forkserver creates a clean single-threaded helper process instead.
os.environ.setdefault('OMP_NUM_THREADS',        '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS',   '1')
os.environ.setdefault('MKL_NUM_THREADS',        '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
# Do NOT set multiprocessing start method — scvelo uses joblib whose loky backend
# is spawn-based and handles fork-safety itself.
warnings.filterwarnings('ignore')

import anndata as ad
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import sparse
import joblib

# ── Config ────────────────────────────────────────────────────────────────────
H5AD_PATH     = Path(os.environ.get('H5AD_PATH',  'output/objects/hsc.h5ad'))
LOOM_DIR      = Path(os.environ.get('LOOM_DIR',   'output/objects/velocyto'))
MODE          = os.environ.get('VELOCITY_MODE', 'dynamical').lower()
OUT_DIR       = Path(os.environ.get('OUT_DIR',    'output/objects/velocity'))
N_JOBS        = int(os.environ.get('N_JOBS', 8))
# Optional: cluster label to use as velocity root (e.g. '0').
# If unset, scvelo auto-detects root cells (can be wrong — check stream plot).
ROOT_CLUSTER    = os.environ.get('ROOT_CLUSTER', None)
# Set NEGATE_VELOCITY=1 to flip all arrows (use when stream plot is inverted).
NEGATE_VELOCITY = os.environ.get('NEGATE_VELOCITY', '0').strip() not in ('', '0', 'false', 'False')

assert MODE in ('stochastic', 'dynamical'), \
    f"VELOCITY_MODE must be 'stochastic' or 'dynamical', got {MODE!r}"

CHKPT_DIR  = OUT_DIR / 'checkpoints'
FIG_DIR    = OUT_DIR / 'figures'
TABLE_DIR  = OUT_DIR / 'tables'
for d in (CHKPT_DIR, FIG_DIR, TABLE_DIR):
    d.mkdir(parents=True, exist_ok=True)

GENES_OF_INTEREST = ['Gata1', 'Pf4', 'Vwf', 'Csf1r', 'Ly6c2',
                     'Cd34', 'Kit', 'Slamf1', 'Mpl']

print(f"scvelo  {scv.__version__}")
print(f"Mode:   {MODE}")
print(f"Output: {OUT_DIR}")

def save_fig(name, dpi=150):
    path = FIG_DIR / name
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close('all')
    print(f"  → {path}")

# ── 1. Merge looms with h5ad ──────────────────────────────────────────────────
MERGED_CHKPT = CHKPT_DIR / 'merged.h5ad'

if MERGED_CHKPT.exists():
    print(f"\n=== 1. Load merged checkpoint ===")
    av = sc.read_h5ad(MERGED_CHKPT)
    print(f"  {av.shape}")
else:
    print(f"\n=== 1. Load h5ad + looms ===")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"  h5ad: {adata.shape}")

    loom_files = sorted(LOOM_DIR.glob('*/*.loom'))
    assert loom_files, f"No loom files found under {LOOM_DIR}"
    print(f"  Loom files: {len(loom_files)}")

    def fix_loom_barcode(bc, sample_id):
        """Convert velocyto barcode to Seurat format.

        velocyto --sampleid writes: Rx_2:AAACTTGCAACTTAGAx
          - prefixes with sampleid + ':'
          - replaces CellRanger gem-group suffix (-1) with 'x'
        Seurat merge format:        Rx_2_AAACTTGCAACTTAGA-1
        """
        core = bc.split(':')[-1]        # strip "Rx_2:" prefix
        core = core.rstrip('x') + '-1'  # "...x" → "...-1"
        return f"{sample_id}_{core}"

    looms = []
    for f in loom_files:
        lm = ad.read_loom(str(f))
        lm.var_names_make_unique()
        sample_id = f.stem              # "Rx_2"
        lm.obs_names = [fix_loom_barcode(bc, sample_id) for bc in lm.obs_names]
        print(f"    {f.name}: {lm.shape}  e.g. {lm.obs_names[0]}")
        looms.append(lm)

    loom_all = ad.concat(looms, merge='same') if len(looms) > 1 else looms[0]
    print(f"  h5ad barcode example:  {adata.obs_names[0]}")
    print(f"  loom barcode example:  {loom_all.obs_names[0]}")

    av = scv.utils.merge(adata.copy(), loom_all)
    print(f"  After merge: {av.shape[0]} cells retained")

    av.write_h5ad(MERGED_CHKPT)
    print(f"  Checkpoint: {MERGED_CHKPT}")

# ── 2. Filter, normalise, moments ─────────────────────────────────────────────
MOMENTS_CHKPT = CHKPT_DIR / f'moments_{MODE}.h5ad'

if MOMENTS_CHKPT.exists():
    print(f"\n=== 2. Load moments checkpoint ===")
    av = sc.read_h5ad(MOMENTS_CHKPT)
    print(f"  {av.shape}")
else:
    print(f"\n=== 2. Filter, normalise, moments ===")
    REPORTERS = ['mCherry', 'DRE']
    reporters_found = [g for g in REPORTERS if g in av.var_names]
    if reporters_found:
        av = av[:, ~av.var_names.isin(reporters_found)].copy()
        print(f"  Excluded reporters: {reporters_found}")

    # X is already SCT log-normalised (DietSeurat data layer) — skip normalisation.
    # spliced/unspliced are raw counts from velocyto looms; filter only on shared counts.
    # TODO: revisit — normalise spliced/unspliced with SCT for consistency.
    scv.pp.filter_genes(av, min_shared_counts=20)
    print(f"  After filter: {av.shape}")
    scv.pp.moments(av, n_pcs=30, n_neighbors=30)

    av.write_h5ad(MOMENTS_CHKPT)
    print(f"  Checkpoint: {MOMENTS_CHKPT}")

# ── 3. Velocity ───────────────────────────────────────────────────────────────
if MODE == 'dynamical':
    DYNAMICS_CHKPT = CHKPT_DIR / f'dynamics_{MODE}.h5ad'

    if DYNAMICS_CHKPT.exists():
        print(f"\n=== 3. Load dynamics checkpoint ===")
        av = sc.read_h5ad(DYNAMICS_CHKPT)
        print(f"  {av.shape}")
    else:
        print(f"\n=== 3. recover_dynamics (dynamical) ===")
        with joblib.parallel_backend('loky', n_jobs=N_JOBS):
            scv.tl.recover_dynamics(av, n_jobs=N_JOBS)
        av.write_h5ad(DYNAMICS_CHKPT)
        print(f"  Checkpoint: {DYNAMICS_CHKPT}")

    print(f"\n=== 4. velocity + velocity_graph (dynamical) ===")
    scv.tl.velocity(av, mode='dynamical')
    if NEGATE_VELOCITY:
        av.layers['velocity'] = -av.layers['velocity']
        print("  Negated velocity layer (NEGATE_VELOCITY=1)")
    # n_jobs=1 to avoid OOM: parallel workers each hold a copy of the graph in RAM
    scv.tl.velocity_graph(av, n_jobs=1)

    if ROOT_CLUSTER is not None:
        root_mask = av.obs[CLUSTER_COL].astype(str) == str(ROOT_CLUSTER)
        av.obs['root_cells'] = root_mask.astype(float)
        print(f"  Root cells: cluster {ROOT_CLUSTER} ({root_mask.sum()} cells)")
        scv.tl.latent_time(av, root_key='root_cells')
    else:
        scv.tl.latent_time(av)

    # Save per-gene kinetics table
    gene_cols = ['fit_alpha', 'fit_beta', 'fit_gamma',
                 'fit_likelihood', 'fit_r2', 'velocity_genes']
    gene_tbl = av.var[[c for c in gene_cols if c in av.var.columns]].copy()
    gene_tbl.to_csv(TABLE_DIR / f'velocity_genes_{MODE}.csv')
    print(f"  Saved: tables/velocity_genes_{MODE}.csv  ({len(gene_tbl)} genes)")

    # Save per-cell latent time
    lt_df = av.obs[['latent_time']].copy()
    if 'sample.name' in av.obs.columns:
        lt_df['sample_name'] = av.obs['sample.name']
    lt_df.to_csv(TABLE_DIR / f'latent_time_{MODE}.csv')
    print(f"  Saved: tables/latent_time_{MODE}.csv")

else:  # stochastic
    print(f"\n=== 3-4. velocity + velocity_graph (stochastic) ===")
    scv.tl.velocity(av, mode='stochastic')
    scv.tl.velocity_graph(av, n_jobs=1)

    gene_cols = ['velocity_genes']
    gene_tbl = av.var[[c for c in gene_cols if c in av.var.columns]].copy()
    gene_tbl.to_csv(TABLE_DIR / f'velocity_genes_{MODE}.csv')
    print(f"  Saved: tables/velocity_genes_{MODE}.csv")

# ── Detect cluster column (needed for root cell assignment + figures) ──────────
CLUSTER_COL = next(
    (c for c in av.obs.columns if 'res.0.9' in c),
    next((c for c in av.obs.columns if 'res.' in c), 'seurat_clusters')
)
print(f"  Cluster column: {CLUSTER_COL!r}")

# ── 5. Figures ────────────────────────────────────────────────────────────────
print(f"\n=== 5. Figures ===")

# Velocity stream
scv.pl.velocity_embedding_stream(
    av, basis='umap', color=CLUSTER_COL, show=False,
    title=f'RNA velocity ({MODE})')
save_fig(f'velocity_stream_{MODE}.png')

# Proportions of splicing states
scv.pl.proportions(av, groupby=CLUSTER_COL, show=False)
save_fig(f'proportions_{MODE}.png')

if MODE == 'dynamical':
    # Latent time UMAP
    scv.pl.scatter(av, color='latent_time', color_map='gnuplot', show=False,
                   title='scVelo latent time')
    save_fig(f'latent_time_{MODE}.png')

    # Phase portraits for selected genes
    genes_ok = [g for g in GENES_OF_INTEREST if g in av.var_names]
    if genes_ok:
        scv.pl.scatter(av, genes_ok[:6], ncols=3, show=False,
                       color='latent_time', color_map='gnuplot')
        save_fig(f'phase_portraits_{MODE}.png')

# ── 6. Save final h5ad ────────────────────────────────────────────────────────
out_h5ad = OUT_DIR / f'hsc_velocity_{MODE}.h5ad'
av.write_h5ad(out_h5ad)
print(f"\n=== Done ===")
print(f"Final h5ad: {out_h5ad}")
