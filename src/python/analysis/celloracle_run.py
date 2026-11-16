"""
celloracle_run.py — In silico TF perturbation of HSC lineages via CellOracle.

CellOracle uses the pySCENIC GRN (adjacencies.csv filtered by validated regulons)
to fit per-gene linear models and simulate the transcriptional effect of setting a
transcription factor's expression to zero (KO). The resulting cell-state shift
vectors are projected onto the UMAP embedding and aggregated per cluster / lineage.

Inputs (all produced by the pipeline or scenic_run.py)
------------------------------------------------------
  output/objects/hsc.h5ad           — SCT-normalized AnnData (cells × genes)
                                    obs: cluster cols, CytoTRACE2 scores, sample.name
                                    obsm: X_umap, X_pca
  output/objects/scenic/adjacencies_full.csv   — pySCENIC GRN: TF, target, importance
  output/objects/scenic/regulons_full.csv      — cisTarget-validated regulons (TF filter)
  output/objects/scenic/auc_mtx_raw_full.csv   — AUC matrix: cells × regulons (TF ranking)
  (all cells: the "full" subset of scenic_run.py; SCENIC_ADJ / SCENIC_REG / SCENIC_AUC override them)

Output layout
-------------
  output/objects/celloracle/
    checkpoints/
      oracle.pkl                    — serialized Oracle object (resume support)
    figures/
      00_grn_overview.png           — edges-per-TF distribution + top-TF bar
      01_umap_overview.png          — cluster + pseudotime UMAP
      02_perturbation_<TF>.png      — shift vectors on UMAP per TF
      03_summary_heatmap.png        — TFs × clusters z-scored shift magnitude
      04_pseudotime_inner_product.png — TF KO effect on differentiation direction
    tables/
      grn_stats.csv                 — per-TF edge count in GRN
      perturbation_shift_scores.csv — cells × TFs absolute shift magnitude
      cluster_shift_summary.csv     — TFs × clusters mean shift (z-scored)
      pseudotime_inner_products.csv — per-TF inner product with pseudotime gradient

Parameters
----------
  TOP_N_TFS    : int  — number of TFs to perturb (ranked by mean AUC, default 20)
  N_PROP       : int  — GRN propagation steps in simulate_shift (default 3)
  SIGMA_CORR   : float — bandwidth for embedding shift smoothing (default 0.05)
  N_PCA_COMPS  : int  — PCA components for knn_imputation (default 30)
  CLUSTER_COL  : str  — obs column used for lineage grouping
  PSEUDOTIME_COL: str — obs column with differentiation pseudotime

Run from the project root, in a Python environment that has celloracle:
  python -m analysis.celloracle_run
celloracle pins pandas<=1.5.3 and gimmemotifs<=0.17.1, which don't install on the project's
Python 3.12, so it is not in pyproject.toml: it needs an environment of its own (Python 3.10).
"""

import sys
import os
import warnings
import pickle
import traceback

warnings.filterwarnings('ignore')

try:
    import celloracle as co
    print(f"celloracle {co.__version__}")
except ImportError:
    sys.exit("celloracle is not installed in this Python (see the note at the top of this file)")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from pathlib import Path
from sklearn.linear_model import LinearRegression

# ── Config ────────────────────────────────────────────────────────────────────
H5AD          = Path('output/objects/hsc.h5ad')
ADJ           = Path(os.environ.get('SCENIC_ADJ', 'output/objects/scenic/adjacencies_full.csv'))
REG           = Path(os.environ.get('SCENIC_REG', 'output/objects/scenic/regulons_full.csv'))
AUC_MTX       = Path(os.environ.get('SCENIC_AUC', 'output/objects/scenic/auc_mtx_raw_full.csv'))
OUT           = Path('output/objects/celloracle')
CHKPT_DIR     = OUT / 'checkpoints'
FIG_DIR       = OUT / 'figures'
TABLE_DIR     = OUT / 'tables'

for d in (CHKPT_DIR, FIG_DIR, TABLE_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Perturbation parameters
TOP_N_TFS     = int(os.environ.get('TOP_N_TFS',    20))
N_PROP        = int(os.environ.get('N_PROP',        3))
N_NEIGHBORS   = int(os.environ.get('N_NEIGHBORS', 200))
N_PCA_COMPS   = int(os.environ.get('N_PCA_COMPS',  30))
SIGMA_CORR    = float(os.environ.get('SIGMA_CORR', 0.05))
CLUSTER_COL   = os.environ.get('CLUSTER_COL',  'sNN_pca.SCT_30_res.1.5')
PSEUDOTIME_COL = os.environ.get('PSEUDOTIME_COL', 'cytotrace2_score')
CONDITION_COL  = 'sample.name'

print("=" * 60)
print("CellOracle TF Perturbation — HSC")
print("=" * 60)
print(f"  H5AD:          {H5AD}")
print(f"  Adjacencies:   {ADJ}")
print(f"  Regulons:      {REG}")
print(f"  AUC matrix:    {AUC_MTX}")
print(f"  Cluster col:   {CLUSTER_COL}")
print(f"  Pseudotime:    {PSEUDOTIME_COL}")
print(f"  Top N TFs:     {TOP_N_TFS}")
print(f"  n_propagation: {N_PROP}")
print(f"  n_pca_comps:   {N_PCA_COMPS}")
print(f"  sigma_corr:    {SIGMA_CORR}")
print()

def save_fig(name, dpi=150):
    path = FIG_DIR / name
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close('all')
    print(f"  → {path}")

# ── 1. Load AnnData ───────────────────────────────────────────────────────────
print("=== 1. Load AnnData ===")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.shape}  ({adata.X.dtype})")

# Resolve cluster column — try exact match then prefix search
if CLUSTER_COL not in adata.obs.columns:
    candidates = [c for c in adata.obs.columns if 'res.' in c]
    if candidates:
        CLUSTER_COL = candidates[-1]
        print(f"  [WARN] CLUSTER_COL not found — using {CLUSTER_COL!r}")
    else:
        raise ValueError(f"No cluster column found. Available: {list(adata.obs.columns)}")

adata.obs[CLUSTER_COL] = adata.obs[CLUSTER_COL].astype(str).astype('category')
clusters  = adata.obs[CLUSTER_COL].cat.categories.tolist()
n_cells   = adata.n_obs
n_genes   = adata.n_vars
print(f"  Clusters ({len(clusters)}): {clusters}")
print(f"  Pseudotime col present: {PSEUDOTIME_COL in adata.obs.columns}")

# CellOracle needs X_umap and X_pca in obsm
assert 'X_umap' in adata.obsm, "X_umap missing from obsm"
assert 'X_pca'  in adata.obsm, "X_pca missing from obsm"

# SCT data in X is log-normalized — store as normalized_count layer so CellOracle
# skips its internal normalization (which would double-transform the data).
# CellOracle uses this layer for linear model fitting.
if sp.issparse(adata.X):
    adata.layers["normalized_count"] = adata.X.toarray()
else:
    adata.layers["normalized_count"] = adata.X.copy()

# ── 2. Build GRN from SCENIC outputs ─────────────────────────────────────────
print("\n=== 2. Build GRN (SCENIC adjacencies + validated regulons) ===")

adj = pd.read_csv(ADJ)   # columns: TF, target, importance
print(f"  Raw adjacencies: {len(adj):,} edges, {adj['TF'].nunique()} TFs")

# Parse regulon TF names from regulons.csv header (format: "TF(+)" or "TF(-)")
# Used to filter adjacencies to only cisTarget-validated TFs
validated_tfs = set()
try:
    with open(REG) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # First field is the regulon name like "Gata1(+)" or "Gata1"
            name = line.split(',')[0].strip('"')
            import re
            m = re.match(r'^([A-Za-z][A-Za-z0-9]+)', name)
            if m:
                validated_tfs.add(m.group(1))
    print(f"  Validated TFs from regulons.csv: {len(validated_tfs)}")
except Exception as e:
    print(f"  [WARN] Could not parse regulons.csv ({e}) — using all adjacency TFs")

if validated_tfs:
    adj_filt = adj[adj['TF'].isin(validated_tfs)].copy()
else:
    adj_filt = adj.copy()
print(f"  Filtered adjacencies: {len(adj_filt):,} edges, {adj_filt['TF'].nunique()} TFs")

# Keep only genes present in the AnnData
adj_filt = adj_filt[
    adj_filt['TF'].isin(adata.var_names) &
    adj_filt['target'].isin(adata.var_names)
].copy()
print(f"  After gene-overlap filter: {len(adj_filt):,} edges")

# Build TFdict: {target_gene: [TF1, TF2, ...]} required by CellOracle
tfdict = (adj_filt
          .groupby('target')['TF']
          .apply(list)
          .to_dict())
print(f"  TFdict: {len(tfdict):,} target genes with regulators")

# GRN stats table
grn_stats = (adj_filt
             .groupby('TF')
             .agg(n_targets=('target', 'count'),
                  mean_importance=('importance', 'mean'))
             .sort_values('n_targets', ascending=False)
             .reset_index())
grn_stats.to_csv(TABLE_DIR / 'grn_stats.csv', index=False)
print(f"  grn_stats.csv saved ({len(grn_stats)} TFs)")

# Figure: GRN overview
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
# Distribution of targets per TF
axes[0].hist(grn_stats['n_targets'], bins=40, color='#4393c3', edgecolor='white')
axes[0].set_xlabel('Targets per TF', fontsize=10)
axes[0].set_ylabel('Count', fontsize=10)
axes[0].set_title('GRN: targets per TF distribution', fontsize=11)
# Top 20 TFs
top20_grn = grn_stats.head(20)
axes[1].barh(top20_grn['TF'][::-1], top20_grn['n_targets'][::-1], color='#4393c3')
axes[1].set_xlabel('Number of target genes', fontsize=10)
axes[1].set_title('Top 20 TFs by target count', fontsize=11)
plt.tight_layout()
save_fig('00_grn_overview.png')

# ── 3. Select top TFs from SCENIC AUC matrix ─────────────────────────────────
print(f"\n=== 3. Select top {TOP_N_TFS} TFs by mean AUC ===")

auc_raw = pd.read_csv(AUC_MTX, index_col=0)
# Columns are like "Gata1(+)" — extract TF names
import re as _re
auc_tfs = {col: _re.sub(r'\([+-]\)$', '', col) for col in auc_raw.columns}
mean_auc = auc_raw.mean(axis=0).rename(index=auc_tfs)
mean_auc = mean_auc.groupby(level=0).max()   # if duplicates, keep highest

# Restrict to TFs present in the filtered GRN and in the AnnData
available = set(adj_filt['TF'].unique()) & set(adata.var_names)
mean_auc = mean_auc[mean_auc.index.isin(available)]
top_tfs = mean_auc.sort_values(ascending=False).head(TOP_N_TFS).index.tolist()

print(f"  TFs selected ({len(top_tfs)}): {top_tfs}")

# ── 4. Build Oracle object ────────────────────────────────────────────────────
ORACLE_PKL = CHKPT_DIR / 'oracle.pkl'

if ORACLE_PKL.exists():
    print(f"\n=== 4. Load Oracle checkpoint ===")
    with open(ORACLE_PKL, 'rb') as fh:
        oracle = pickle.load(fh)
    print(f"  Loaded from {ORACLE_PKL}")
else:
    print(f"\n=== 4. Build Oracle object ===")
    oracle = co.Oracle()

    # Import AnnData — pass the obs metadata and embeddings.
    # normalized_count layer contains SCT log-normalized values.
    oracle.import_anndata_as_raw_count(
        adata            = adata,
        cluster_column_name = CLUSTER_COL,
        embedding_name   = 'X_umap',
    )

    # Import GRN
    oracle.import_TF_data(TFdict=tfdict)
    print(f"  GRN imported: {len(tfdict):,} regulated genes")

    # PCA + KNN imputation (smoothing)
    oracle.perform_PCA()
    n_comps_avail = oracle.pca.components_.shape[0]
    n_comps = min(N_PCA_COMPS, n_comps_avail)
    print(f"  KNN imputation with {n_comps} PCA components")
    oracle.knn_imputation(n_pca_components=n_comps)

    # Fit linear regression GRN models (α = L2 regularisation)
    print("  Fitting GRN linear models (alpha=10) ...")
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=False)

    with open(ORACLE_PKL, 'wb') as fh:
        pickle.dump(oracle, fh)
    print(f"  Oracle saved: {ORACLE_PKL}")

# ── 5. Overview UMAPs ─────────────────────────────────────────────────────────
print(f"\n=== 5. Overview UMAPs ===")
panels = [CLUSTER_COL]
if PSEUDOTIME_COL in adata.obs.columns:
    panels.append(PSEUDOTIME_COL)
if CONDITION_COL in adata.obs.columns:
    panels.append(CONDITION_COL)

fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 5))
axes = [axes] if len(panels) == 1 else list(axes)
cmaps = {PSEUDOTIME_COL: 'RdYlBu_r'}
for ax, col in zip(axes, panels):
    sc.pl.umap(adata, color=col, ax=ax, show=False,
               color_map=cmaps.get(col), title=col, legend_loc='right margin')
plt.tight_layout()
save_fig('01_umap_overview.png')

# ── 6. Pseudotime gradient in UMAP space ─────────────────────────────────────
# Fit a linear model: UMAP ~ pseudotime to get the differentiation direction
# vector. Used later to compute the inner product (does KO advance/block diff?).
pt_gradient = None
if PSEUDOTIME_COL in adata.obs.columns:
    print(f"\n=== 6. Pseudotime gradient ===")
    # CytoTRACE2: high score = stem-like. Invert so high pt = differentiated.
    pt = 1.0 - adata.obs[PSEUDOTIME_COL].fillna(0.5).values
    umap = adata.obsm['X_umap']
    lr = LinearRegression(fit_intercept=True).fit(pt.reshape(-1, 1), umap)
    pt_gradient = lr.coef_[0]
    pt_gradient = pt_gradient / np.linalg.norm(pt_gradient)
    print(f"  Pseudotime gradient vector (UMAP): {pt_gradient}")

# ── 7. Perturbation simulation loop ──────────────────────────────────────────
print(f"\n=== 7. Perturbation loop ({len(top_tfs)} TFs) ===")

cell_ids  = oracle.adata.obs_names.tolist()
clust_lab = oracle.adata.obs[CLUSTER_COL].astype(str).values

# Storage
shift_scores  = pd.DataFrame(index=cell_ids)   # cells × TFs
cluster_means = {}                              # TF → {cluster: mean_shift}
pt_inner_prod = {}                              # TF → scalar inner product

for i, tf in enumerate(top_tfs):
    print(f"  [{i+1:2d}/{len(top_tfs)}] {tf} ...", end=' ', flush=True)

    try:
        oracle.simulate_shift(
            perturb_condition = {tf: 0.0},    # KO: set TF to 0
            n_propagation     = N_PROP,
        )
        oracle.estimate_transition_prob(
            n_neighbors      = N_NEIGHBORS,
            knn_random       = True,
            sampled_fraction = 1.0,
        )
        oracle.calculate_embedding_shift(sigma_corr=SIGMA_CORR)

        # Shift vectors in UMAP space — cells × 2
        delta = oracle.delta_embedding            # (n_cells, 2)
        shift_mag = np.linalg.norm(delta, axis=1)

        shift_scores[tf] = shift_mag

        # Per-cluster mean shift magnitude
        cmeans = {}
        for cl in np.unique(clust_lab):
            mask = clust_lab == cl
            cmeans[cl] = float(shift_mag[mask].mean())
        cluster_means[tf] = cmeans

        # Inner product with pseudotime gradient (positive = pro-differentiation)
        if pt_gradient is not None:
            ip = float(np.mean(delta @ pt_gradient))
            pt_inner_prod[tf] = ip

        # Per-TF perturbation UMAP figure
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        # Left: shift magnitude
        sc.pl.umap(oracle.adata, color=CLUSTER_COL, ax=axes[0], show=False,
                   title='Clusters', legend_loc='right margin')
        axes[1].scatter(
            oracle.adata.obsm['X_umap'][:, 0],
            oracle.adata.obsm['X_umap'][:, 1],
            c=shift_mag, cmap='magma', s=2, alpha=0.6,
            vmin=np.percentile(shift_mag, 2),
            vmax=np.percentile(shift_mag, 98),
        )
        sm = plt.cm.ScalarMappable(
            cmap='magma',
            norm=mcolors.Normalize(
                vmin=np.percentile(shift_mag, 2),
                vmax=np.percentile(shift_mag, 98),
            )
        )
        plt.colorbar(sm, ax=axes[1], label='Shift magnitude')
        axes[1].set_title(f'{tf} KO — shift magnitude')
        axes[1].axis('off')
        plt.suptitle(f'TF KO: {tf}', fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig(f'02_perturbation_{tf}.png')

        print("done")

    except Exception as exc:
        print(f"FAILED: {exc}")
        traceback.print_exc()
        cluster_means[tf] = {}
        pt_inner_prod[tf] = np.nan

# ── 8. Save per-cell shift score table ───────────────────────────────────────
print(f"\n=== 8. Save tables ===")

shift_scores.to_csv(TABLE_DIR / 'perturbation_shift_scores.csv')
print(f"  perturbation_shift_scores.csv  {shift_scores.shape}")

# Cluster × TF summary (raw mean shift)
cl_df = pd.DataFrame(cluster_means).T    # TFs × clusters
cl_df.index.name = 'TF'
cl_df.to_csv(TABLE_DIR / 'cluster_shift_summary.csv')
print(f"  cluster_shift_summary.csv  {cl_df.shape}")

# Pseudotime inner products
if pt_inner_prod:
    ip_df = (pd.Series(pt_inner_prod, name='pt_inner_product')
             .reset_index()
             .rename(columns={'index': 'TF'})
             .sort_values('pt_inner_product'))
    ip_df.to_csv(TABLE_DIR / 'pseudotime_inner_products.csv', index=False)
    print(f"  pseudotime_inner_products.csv  {ip_df.shape}")

# ── 9. Summary heatmap — TFs × clusters ──────────────────────────────────────
print(f"\n=== 9. Summary heatmap ===")

valid_tfs = [tf for tf in top_tfs if tf in cl_df.index and cl_df.loc[tf].notna().any()]
if valid_tfs:
    hm_data = cl_df.loc[valid_tfs].fillna(0)

    # Z-score within each TF row (how much is each cluster shifted relative to others?)
    from scipy.stats import zscore
    hm_z = hm_data.apply(lambda row: zscore(row, nan_policy='omit'), axis=1)
    hm_z = hm_z.fillna(0)

    # Sort clusters by mean z-score (most-shifted first)
    col_order = hm_z.mean(axis=0).sort_values(ascending=False).index.tolist()
    # Sort TFs by maximum cluster shift
    row_order = hm_z.max(axis=1).sort_values(ascending=False).index.tolist()

    hm_plot = hm_z.loc[row_order, col_order]

    fig_h = max(5, len(valid_tfs) * 0.45)
    fig_w = max(8, len(clusters) * 0.65)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        hm_plot,
        ax        = ax,
        cmap      = 'RdBu_r',
        center    = 0,
        vmin      = -3, vmax = 3,
        linewidths = 0.4,
        linecolor  = 'white',
        cbar_kws   = dict(label='Z-score (shift magnitude)', shrink=0.6),
        xticklabels = True,
        yticklabels = True,
    )
    ax.set_xlabel('Cluster', fontsize=11)
    ax.set_ylabel('TF (KO)', fontsize=11)
    ax.set_title('TF KO effect per cluster\n(z-scored mean shift magnitude)',
                 fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=9)
    plt.tight_layout()
    save_fig('03_summary_heatmap.png')

# ── 10. Pseudotime inner product bar plot ─────────────────────────────────────
if pt_inner_prod and any(np.isfinite(v) for v in pt_inner_prod.values()):
    print(f"\n=== 10. Pseudotime inner product ===")

    ip_plot = (pd.Series(pt_inner_prod)
               .dropna()
               .sort_values())

    colours = ['#b2182b' if v > 0 else '#2166ac' for v in ip_plot.values]
    fig, ax = plt.subplots(figsize=(8, max(4, len(ip_plot) * 0.45)))
    ax.barh(ip_plot.index, ip_plot.values, color=colours, edgecolor='white')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel('Inner product with pseudotime gradient\n'
                  '(+) pro-differentiation  |  (−) blocks differentiation',
                  fontsize=10)
    ax.set_title('TF KO: effect on differentiation direction\n'
                 '(CytoTRACE2 pseudotime gradient)', fontsize=11, fontweight='bold')
    ax.tick_params(axis='y', labelsize=9)
    plt.tight_layout()
    save_fig('04_pseudotime_inner_product.png')

# ── Done ─────────────────────────────────────────────────────────────────────
print(f"\n=== Done ===")
print(f"Outputs in {OUT}:")
for d in (FIG_DIR, TABLE_DIR):
    for f in sorted(d.iterdir()):
        print(f"  {str(f.relative_to(OUT)):55s}  {f.stat().st_size / 1e3:8.1f} KB")
