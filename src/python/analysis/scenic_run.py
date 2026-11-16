"""
scenic_run.py — pySCENIC for HSC compartment (mm10), run per condition.

Runs the full SCENIC pipeline independently for each condition subset
(mch_neg, mch_pos, wt_ctrl) defined by the sample.name metadata column,
and once on all cells (full). All output files carry a per-condition suffix.

Steps per subset:
  1. subset h5ad → hsc_{suffix}_expr.loom
  2. GRN          → adjacencies_{suffix}.csv
  3. cisTarget    → regulons_{suffix}.csv          (RAW)
  4. filter       → regulons_clean_{suffix}.dat
  5. tidy tables  → regulons_tidy_raw_{suffix}.csv
                    regulons_tidy_clean_{suffix}.csv
  6. AUCell       → auc_mtx_raw_{suffix}.csv
  7. AUCell       → auc_mtx_clean_{suffix}.csv

Each step is skipped if its output already exists.

Two modes, because the steps need different Python environments (run from the project root):
  prepare [subset ...]   step 1, in the project's environment (it reads the h5ad with anndata):
                           uv run --no-sync python -m analysis.scenic_run prepare
  infer [subset ...]     steps 2-7, inside the pySCENIC container (scripts/containers/pyscenic.def):
                           python src/python/analysis/scenic_run.py infer
Without subsets, every subset is processed. scripts/slurm/analysis/scenic.sh runs both modes,
one array task per subset. N_JOBS sets the worker count (the job passes its CPU count).
"""
import gc
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config ────────────────────────────────────────────────────────────────────
H5AD   = Path('output/objects/hsc.h5ad')
DB_DIR = Path('data/scenic')
OUT    = Path('output/objects/scenic')

TF_LIST = DB_DIR / 'mm_mgi_tfs.txt'
MOTIFS  = DB_DIR / 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
DB_500  = DB_DIR / 'mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
DB_10K  = DB_DIR / 'mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'

SUBSETS = {
    'mch_neg': 'HSC_Het_CD41+mCh-',
    'mch_pos': 'HSC_Het_CD41+mCh+',
    'wt_ctrl': 'HSC_WT_CD41+',
    'full':    None,   # all cells
}

# Filtering thresholds
NES_THRESHOLD   = 3.0
MIN_GENES       = 10
ACTIVATING_ONLY = True

N_JOBS   = int(os.environ.get('N_JOBS', 6))
PYSCENIC = ['pyscenic']


def subset_paths(suffix):
    return dict(
        LOOM           = OUT / f'hsc_{suffix}_expr.loom',
        ADJ            = OUT / f'adjacencies_{suffix}.csv',
        REG            = OUT / f'regulons_{suffix}.csv',
        REG_CLEAN      = OUT / f'regulons_clean_{suffix}.dat',
        REG_TIDY_RAW   = OUT / f'regulons_tidy_raw_{suffix}.csv',
        REG_TIDY_CLEAN = OUT / f'regulons_tidy_clean_{suffix}.csv',
        AUC_RAW        = OUT / f'auc_mtx_raw_{suffix}.csv',
        AUC_CLEAN      = OUT / f'auc_mtx_clean_{suffix}.csv',
    )


def done(path):
    return path.exists() and path.stat().st_size > 0


# ── 1. loom (mode: prepare) ───────────────────────────────────────────────────
def prepare(subsets):
    import anndata as ad
    import loompy
    import scipy.sparse as sp

    print("Loading h5ad ...")
    adata_full = ad.read_h5ad(H5AD)
    print(adata_full)

    for suffix in subsets:
        sample_name = SUBSETS[suffix]
        LOOM = subset_paths(suffix)['LOOM']
        print(f"\n--- 1. loom: {suffix} ({sample_name}) ---")
        if sample_name is None:
            adata = adata_full.copy()
        else:
            adata = adata_full[adata_full.obs['sample.name'] == sample_name].copy()
        print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        if done(LOOM):
            with loompy.connect(str(LOOM), mode='r', validate=False) as _ds:
                loom_shape = (_ds.shape[1], _ds.shape[0])
            if loom_shape != (adata.n_obs, adata.n_vars):
                print(f"  [STALE] loom {loom_shape} != h5ad ({adata.n_obs}, {adata.n_vars}) — regenerating")
                LOOM.unlink()
            else:
                print(f"  [SKIP] {LOOM.name} ({loom_shape[0]:,} cells × {loom_shape[1]:,} genes)")
        if not LOOM.exists():
            X = adata.X
            if sp.issparse(X):
                X = X.toarray()
            X = X.astype(np.float32).T
            loompy.create(
                str(LOOM),
                layers    = {'': X},
                row_attrs = {'Gene':   np.array(adata.var_names)},
                col_attrs = {'CellID': np.array(adata.obs_names)},
            )
            print(f"  Written: {LOOM.name}  ({LOOM.stat().st_size / 1e9:.2f} GB)")
        del adata; gc.collect()


# ── 2-7. GRN, cisTarget, filtering, AUCell (mode: infer) ─────────────────────
def infer(subsets):
    import pickle
    import loompy as lp

    for suffix in subsets:
        print(f"\n{'='*60}")
        print(f"  Subset: {suffix}  ({SUBSETS[suffix]})")
        print(f"{'='*60}")
        p = subset_paths(suffix)
        if not done(p['LOOM']):
            sys.exit(f"{p['LOOM']} missing: run the prepare mode first")

        # ── 2. GRN (GRNBoost2 through the pyscenic CLI, seeded) ──────────────────
        print("\n--- 2. GRN ---")
        if done(p['ADJ']):
            print(f"  [SKIP] {p['ADJ'].name} exists")
        else:
            subprocess.run(PYSCENIC + [
                'grn', str(p['LOOM']), str(TF_LIST),
                '--output', str(p['ADJ']),
                '--num_workers', str(N_JOBS),
                '--seed', '42',
            ], check=True)
            print(f"  Done: {p['ADJ'].name}")

        # ── 3. cisTarget ──────────────────────────────────────────────────────────
        print("\n--- 3. cisTarget ---")
        if done(p['REG']):
            print(f"  [SKIP] {p['REG'].name} exists")
        else:
            subprocess.run(PYSCENIC + [
                'ctx',
                str(p['ADJ']), str(DB_500), str(DB_10K),
                '--annotations_fname',    str(MOTIFS),
                '--expression_mtx_fname', str(p['LOOM']),
                '--output', str(p['REG']),
                '--num_workers', str(N_JOBS),
                '--mask_dropouts',
            ], check=True)
            print(f"  Done: {p['REG'].name}")

        # ── 4. Filter regulons ────────────────────────────────────────────────────
        print("\n--- 4. Filter regulons ---")
        if done(p['REG_CLEAN']):
            print(f"  [SKIP] {p['REG_CLEAN'].name} exists")
        else:
            from pyscenic.utils import load_motifs
            from pyscenic.prune import df2regulons

            motifs_df = load_motifs(str(p['REG']))

            def derive_regulons(motifs):
                E   = 'Enrichment'
                ctx = motifs[(E, 'Context')]
                keep = (~ctx.apply(lambda c: 'weight>50.0%' in c)
                        & (motifs[(E, 'NES')] >= NES_THRESHOLD))
                if ACTIVATING_ONLY:
                    keep &= ctx.apply(lambda c: 'activating' in c)
                kept = motifs[keep]
                assert not kept.empty, "No motifs survived filtering"
                regs = [r for r in df2regulons(kept) if len(r) >= MIN_GENES]
                if ACTIVATING_ONLY:
                    regs = [r.rename(r.transcription_factor) for r in regs]
                return regs

            regulons = derive_regulons(motifs_df)
            print(f"  {len(motifs_df):,} motif rows → {len(regulons)} clean regulons")
            with open(p['REG_CLEAN'], 'wb') as f:
                pickle.dump(regulons, f)
            print(f"  Written: {p['REG_CLEAN'].name}")
            del motifs_df, regulons; gc.collect()

        # ── 5. Tidy tables ────────────────────────────────────────────────────────
        print("\n--- 5. Tidy tables ---")
        if done(p['REG_TIDY_RAW']):
            print(f"  [SKIP] {p['REG_TIDY_RAW'].name} exists")
        else:
            from pyscenic.utils import load_motifs
            E = 'Enrichment'
            motifs_df = load_motifs(str(p['REG']))
            rows = []
            for (tf, motif_id), row in motifs_df.iterrows():
                nes, auc = row[(E, 'NES')], row[(E, 'AUC')]
                for gene, weight in row[(E, 'TargetGenes')]:
                    rows.append((tf, motif_id, nes, auc, gene, float(weight)))
            pd.DataFrame(rows, columns=['TF', 'motif_id', 'NES', 'AUC', 'target', 'weight']
                         ).to_csv(p['REG_TIDY_RAW'], index=False)
            print(f"  {len(rows):,} rows → {p['REG_TIDY_RAW'].name}")
            del motifs_df, rows; gc.collect()

        if done(p['REG_TIDY_CLEAN']):
            print(f"  [SKIP] {p['REG_TIDY_CLEAN'].name} exists")
        else:
            with open(p['REG_CLEAN'], 'rb') as f:
                regulons = pickle.load(f)
            rows = [
                (r.transcription_factor, gene, weight)
                for r in regulons
                for gene, weight in r.gene2weight.items()
            ]
            clean_tidy = pd.DataFrame(rows, columns=['TF', 'target', 'weight'])
            clean_tidy.to_csv(p['REG_TIDY_CLEAN'], index=False)
            print(f"  {len(clean_tidy):,} rows ({clean_tidy['TF'].nunique()} TFs) → {p['REG_TIDY_CLEAN'].name}")
            del regulons, rows, clean_tidy; gc.collect()

        # ── 6. AUCell — RAW ───────────────────────────────────────────────────────
        print("\n--- 6. AUCell (raw) ---")
        if done(p['AUC_RAW']):
            print(f"  [SKIP] {p['AUC_RAW'].name} exists")
        else:
            subprocess.run(PYSCENIC + [
                'aucell', str(p['LOOM']), str(p['REG']),
                '--output', str(p['AUC_RAW']),
                '--num_workers', str(N_JOBS),
            ], check=True)
            print(f"  Done: {p['AUC_RAW'].name}")

        # ── 7. AUCell — CLEAN (Python API — CLI can't load pickle) ───────────────
        print("\n--- 7. AUCell (clean) ---")
        if done(p['AUC_CLEAN']):
            print(f"  [SKIP] {p['AUC_CLEAN'].name} exists")
        else:
            from pyscenic.aucell import aucell

            with open(p['REG_CLEAN'], 'rb') as f:
                regulons = pickle.load(f)
            with lp.connect(str(p['LOOM']), mode='r', validate=False) as ds:
                expr_df = pd.DataFrame(
                    ds[:, :].T,
                    columns=ds.ra['Gene'],
                    index=ds.ca['CellID'],
                )
            auc_mtx = aucell(expr_df, regulons, num_workers=N_JOBS)
            auc_mtx.to_csv(p['AUC_CLEAN'])
            print(f"  Done: {p['AUC_CLEAN'].name} ({len(regulons)} regulons)")
            del regulons, expr_df, auc_mtx; gc.collect()


if __name__ == '__main__':
    mode, *chosen = sys.argv[1:] or ['']
    unknown = [s for s in chosen if s not in SUBSETS]
    if mode not in ('prepare', 'infer') or unknown:
        sys.exit(f"usage: scenic_run.py prepare|infer [{'|'.join(SUBSETS)} ...]"
                 + (f"  (unknown subset: {', '.join(unknown)})" if unknown else ''))
    OUT.mkdir(parents=True, exist_ok=True)
    (prepare if mode == 'prepare' else infer)(chosen or list(SUBSETS))

    print("\n=== Done ===")
    print(f"Outputs in {OUT}:")
    for f in sorted(OUT.iterdir()):
        if f.is_file():
            print(f"  {f.name:50s}  {f.stat().st_size / 1e6:8.1f} MB")
