"""
run_paga.py — Parameterised PAGA computation for the HSC compartment.

Runs PAGA globally (all cells) and per-condition (WT / mCh− / mCh+) for a
given cluster resolution, then exports the two CSVs consumed by
paga-traj-analysis-{res}.qmd.

Env vars
--------
RESOLUTION   Seurat cluster resolution string, e.g. "0.9", "1", "1.5"
             Default: "0.9"
H5AD_PATH    Path to input h5ad (default: output/objects/hsc.h5ad)
OUT_PAGA     Directory for global PAGA outputs
             Default: output/objects/paga
OUT_COND     Directory for per-condition PAGA outputs
             Default: output/objects/paga_conditions
N_PCS        Number of PCs to use for neighbour graph (default: 50)
N_NEIGH      Number of neighbours (default: 30)

Output files (per resolution)
-------------------------------
<OUT_PAGA>/paga_edges_{CLUST_COL}.csv
    cluster_A, cluster_B, weight   — global PAGA, all edge pairs > 0

<OUT_COND>/paga_edge_comparison_{CLUST_COL}.csv
    pair, WT, mCh-, mCh+, max_weight — per-condition PAGA, union of all pairs
"""
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

sc.settings.verbosity = 1

# ── Config ────────────────────────────────────────────────────────────────────
RESOLUTION = os.environ.get('RESOLUTION', '0.9')
N_PCS      = int(os.environ.get('N_PCS',   50))
N_NEIGH    = int(os.environ.get('N_NEIGH', 30))

H5AD_PATH = Path(os.environ.get('H5AD_PATH',
                                 'output/objects/hsc.h5ad'))
OUT_PAGA  = Path(os.environ.get('OUT_PAGA',
                                 'output/objects/paga'))
OUT_COND  = Path(os.environ.get('OUT_COND',
                                 'output/objects/paga_conditions'))

OUT_PAGA.mkdir(parents=True, exist_ok=True)
OUT_COND.mkdir(parents=True, exist_ok=True)

CLUST_COL   = f'sNN_pca.SCT_30_res.{RESOLUTION}'
COND_MAP    = {
    'HSC_WT_CD41+'      : 'WT',
    'HSC_Het_CD41+mCh-' : 'mCh-',
    'HSC_Het_CD41+mCh+' : 'mCh+',
}
COND_ORDER  = ['WT', 'mCh-', 'mCh+']

print(f"Resolution : {RESOLUTION}")
print(f"Cluster col: {CLUST_COL}")
print(f"H5AD       : {H5AD_PATH}")
print(f"Out (global): {OUT_PAGA}")
print(f"Out (cond)  : {OUT_COND}")

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("\n=== 1. Load h5ad ===")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.shape}")

assert CLUST_COL in adata.obs.columns, \
    f"Cluster column {CLUST_COL!r} not found in obs. " \
    f"Available: {[c for c in adata.obs.columns if 'res.' in c]}"

adata.obs[CLUST_COL] = adata.obs[CLUST_COL].astype(str).astype('category')
adata.obs['condition'] = adata.obs['sample.name'].map(COND_MAP).astype('category')

n_clusters = adata.obs[CLUST_COL].nunique()
print(f"  Clusters ({CLUST_COL}): {n_clusters}")
print(f"  Conditions: {adata.obs['condition'].value_counts().to_dict()}")

# ── 2. Global neighbours + PAGA ───────────────────────────────────────────────
print(f"\n=== 2. Global PAGA ===")
sc.pp.neighbors(adata, n_neighbors=N_NEIGH, n_pcs=N_PCS,
                use_rep='X_pca', random_state=42)
sc.tl.paga(adata, groups=CLUST_COL)

conn  = adata.uns['paga']['connectivities']
cats  = adata.obs[CLUST_COL].cat.categories
n     = len(cats)

rows = []
for i in range(n):
    for j in range(i + 1, n):
        w = float(conn[i, j])
        if w > 0:
            rows.append({
                'cluster_A': str(cats[i]),
                'cluster_B': str(cats[j]),
                'weight'   : round(w, 4)
            })

edge_df = pd.DataFrame(rows).sort_values('weight', ascending=False).reset_index(drop=True)
out_global = OUT_PAGA / f'paga_edges_{CLUST_COL}.csv'
edge_df.to_csv(out_global, index=False)
print(f"  {len(edge_df)} edges → {out_global}")

# ── 3. Per-condition neighbours + PAGA ───────────────────────────────────────
print(f"\n=== 3. Per-condition PAGA ===")

def build_paga(adata_full, condition, clust_col):
    sub = adata_full[adata_full.obs['condition'] == condition].copy()
    sub.obs[clust_col] = sub.obs[clust_col].cat.remove_unused_categories()
    sc.pp.neighbors(sub, n_neighbors=N_NEIGH, n_pcs=N_PCS,
                    use_rep='X_pca', random_state=42)
    sc.tl.paga(sub, groups=clust_col)
    n_cl = sub.obs[clust_col].nunique()
    print(f"  {condition}: {sub.n_obs:,} cells, {n_cl} clusters")
    return sub

def extract_edges(sub, clust_col, cond_label):
    conn = sub.uns['paga']['connectivities']
    cats = sub.obs[clust_col].cat.categories
    n    = len(cats)
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            w = float(conn[i, j])
            rows.append({
                'pair'     : f'{cats[i]}-{cats[j]}',
                'condition': cond_label,
                'weight'   : round(w, 4)
            })
    return pd.DataFrame(rows)

cond_dfs = []
for cond in COND_ORDER:
    sub = build_paga(adata, cond, CLUST_COL)
    cond_dfs.append(extract_edges(sub, CLUST_COL, cond))

all_edges = pd.concat(cond_dfs, ignore_index=True)

# Pivot: rows = cluster pairs, cols = conditions
edge_wide = (
    all_edges
    .pivot_table(index='pair', columns='condition',
                 values='weight', fill_value=0.0)
    [COND_ORDER]
    .reset_index()
)
edge_wide.columns.name = None
edge_wide['max_weight'] = edge_wide[COND_ORDER].max(axis=1)
edge_wide = edge_wide.sort_values('max_weight', ascending=False).reset_index(drop=True)

out_cond = OUT_COND / f'paga_edge_comparison_{CLUST_COL}.csv'
edge_wide.to_csv(out_cond, index=False)
print(f"  {len(edge_wide)} edge pairs → {out_cond}")

print(f"\n=== Done — res {RESOLUTION} ===")
print(f"  {out_global}")
print(f"  {out_cond}")
