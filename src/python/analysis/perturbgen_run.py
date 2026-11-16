"""
perturbgen_run.py — In silico TF perturbation of HSC lineages using PerturbGen.

PerturbGen (Lotfollahi-lab/Perturbgen) is a Geneformer-based generative Seq2Seq
transformer pre-trained on >95 M single cells.  Given a source cell state (HSC
progenitors) and a trajectory defined by pseudotime, it predicts the perturbed
transcriptome at later timepoints when a TF is knocked out.

Architecture summary
--------------------
  Source tokens  (progenitor cells, pseudotime bin 0)
       ↓ masked token generation
  Target tokens  (differentiated cells, bins 1–2)
       ↓ count decoder
  Predicted expression profiles  →  compare with unperturbed ground truth

Tokenisation uses Geneformer-style rank-based encoding: genes are ranked by
expression (median-normalised), generating a ranked list of Ensembl gene tokens.
Because the model vocabulary is human, mouse TFs are mapped to human 1:1
orthologues before specifying perturbation tokens.

Inputs
------
  output/objects/hsc.h5ad                    — SCT-normalised AnnData
  output/objects/scenic/auc_mtx_raw_full.csv      — SCENIC AUC for TF ranking (all cells; SCENIC_AUC overrides)
  pretraining_cohort/…ckpt                — downloaded PerturbGen encoder
  data/perturbgen/                     — Geneformer vocabulary (pkl files)

Pipeline steps (each skipped if checkpoint exists)
---------------------------------------------------
  0. Install PerturbGen from GitHub
  1. Download pre-trained encoder checkpoint
  2. Geneformer vocabulary files (gene medians, token dict, Ensembl map)
  3. Preprocess h5ad: mouse→human orthologs, Ensembl IDs, pseudotime bins
  4. perturbgen tokenise  (rank-based, stratified pairing by pseudotime bin)
  5. perturbgen train-decoder  (count decoder on top of frozen encoder)
  6. Select top TFs from SCENIC; build perturbation YAML configs
  7. perturbgen generate (val.py) per TF → per-TF AnnData
  8. Post-perturbation analysis: DEG counts, UMAP, comparison with CellOracle

Output layout
-------------
  output/objects/perturbgen/
    tokenized/
      dataset_src/                           — tokenised source dataset
      dataset_tgt_t1/ dataset_tgt_t2/        — target timepoint datasets
      h5ad_src.h5ad  h5ad_tgt/               — paired h5ads
      token_id_to_genename.pkl
      tokenid_to_rowid.pkl
    checkpoints/
      encoder/                               — pre-trained PerturbGen encoder
      count_decoder/best_model.ckpt
    configs/
      perturbation_<TF>.yaml                — one config per TF
    perturbation_results/
      adata_g<ensembl_id>_ssrc.h5ad         — predicted perturbed h5ads
    figures/
      00_preprocessing_umap.png
      01_perturbation_<TF>.png
      02_summary_heatmap.png
      03_pseudotime_inner_product.png
      04_comparison_celloracle_vs_perturbgen.png
    tables/
      ortho_used.csv
      perturbation_deg_summary.csv          — TF × cluster DEG count
      cluster_effect_summary.csv
      pseudotime_inner_products.csv

Comparison with CellOracle
--------------------------
  CellOracle: GRN linear propagation — fast, interpretable, no training needed.
  PerturbGen: generative transformer — captures higher-order co-expression
              patterns, predicts expression at future timepoints, slower but
              does not rely on a pre-specified GRN.

Run
---
  sbatch scripts/slurm/analysis/perturbgen.sh               (from the project root)
  TOP_N_TFS=10 sbatch scripts/slurm/analysis/perturbgen.sh
Dependencies: the project's `perturbgen` dependency group (uv sync --group perturbgen, once).
"""

import sys, os, subprocess, warnings, yaml, json, pickle, traceback, re
import shutil
warnings.filterwarnings('ignore')

# Disable wandb before any import touches it
os.environ.setdefault('WANDB_MODE', 'disabled')
os.environ.setdefault('WANDB_SILENT', 'true')

# ── Config ────────────────────────────────────────────────────────────────────
H5AD          = 'output/objects/hsc.h5ad'
AUC_MTX       = os.environ.get('SCENIC_AUC', 'output/objects/scenic/auc_mtx_raw_full.csv')
REPO_DIR      = 'data/perturbgen/Perturbgen'

# Put the cloned PerturbGen repository on the path of this process and of its subprocesses
# (tokenise, train-decoder, val.py), so `perturbgen`, `scmaskgit` and the other repo-local
# packages import from the clone; their dependencies come from the `perturbgen` group.
_existing_pythonpath = os.environ.get('PYTHONPATH', '')
_new_pythonpath = ':'.join(filter(None, [REPO_DIR, _existing_pythonpath]))
os.environ['PYTHONPATH'] = _new_pythonpath
sys.path.insert(0, REPO_DIR)
VOCAB_DIR     = 'data/perturbgen/vocab'      # Geneformer pkl files
ENCODER_CKPT  = 'data/perturbgen/encoder.ckpt'
OUT           = 'output/objects/perturbgen'
TOK_DIR       = f'{OUT}/tokenized'
CKPT_DIR      = f'{OUT}/checkpoints'
CFG_DIR       = f'{OUT}/configs'
PERT_DIR      = f'{OUT}/perturbation_results'
FIG_DIR       = f'{OUT}/figures'
TABLE_DIR     = f'{OUT}/tables'

TOP_N_TFS     = int(os.environ.get('TOP_N_TFS',    15))
N_HVG         = int(os.environ.get('N_HVG',      2000))
BATCH_SIZE    = int(os.environ.get('BATCH_SIZE',    32))
EPOCHS_DEC    = int(os.environ.get('EPOCHS_DEC',    10))   # count decoder epochs
CLUSTER_COL   = os.environ.get('CLUSTER_COL',  'sNN_pca.SCT_30_res.1.5')
PSEUDOTIME_COL = os.environ.get('PSEUDOTIME_COL', 'cytotrace2_score')
CONDITION_COL  = 'sample.name'
N_WORKERS      = int(os.environ.get('N_WORKERS', 4))

# Pre-training cohort checkpoint — available via Git LFS in the PerturbGen repo
ENCODER_REPO_PATH = (
    'pretraining_cohort/'
    '20250709_1223_cellgen_train_masking_lr_5e-05_wd_1e-06_'
    'batch_64_ptime_pos_sin_m_pow_tp_1-2-3_s_42-epoch=00.ckpt'
)

for d in [VOCAB_DIR, TOK_DIR, CKPT_DIR, CFG_DIR, PERT_DIR, FIG_DIR, TABLE_DIR,
          f'{CKPT_DIR}/encoder', f'{CKPT_DIR}/count_decoder']:
    os.makedirs(d, exist_ok=True)

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
from scipy.stats import zscore

CELLORACLE_TABLE = Path('output/objects/celloracle/tables/cluster_shift_summary.csv')

print("=" * 60)
print("PerturbGen TF Perturbation — HSC")
print("=" * 60)
print(f"  H5AD:        {H5AD}")
print(f"  Repo:        {REPO_DIR}")
print(f"  Cluster col: {CLUSTER_COL}")
print(f"  Top N TFs:   {TOP_N_TFS}")
print(f"  N HVG:       {N_HVG}")
print(f"  Decoder epochs: {EPOCHS_DEC}")
print()

def save_fig(name, dpi=150):
    path = f'{FIG_DIR}/{name}'
    plt.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close('all')
    print(f"  → {path}")

def run(cmd, **kwargs):
    """Run subprocess with live output."""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, check=True, **kwargs)
    return result

# ── 0. PerturbGen code ────────────────────────────────────────────────────────
# A clone of the PerturbGen repository in data/perturbgen/ (a link to a folder outside the
# project, like every input); cloned there on the first run if missing.
print("=== 0. PerturbGen code ===")
if not Path('data/perturbgen').is_dir():
    sys.exit("data/perturbgen missing: link the PerturbGen folder first (tools/link-input.sh <path>/perturbgen)")
if Path(f'{REPO_DIR}/perturbgen/__init__.py').exists():
    print(f"  [SKIP] Already cloned at {REPO_DIR}")
else:
    run(['git', 'clone', '--depth=1',
         'https://github.com/Lotfollahi-lab/Perturbgen.git', REPO_DIR])
    print("  Cloned repository")
import perturbgen
print("  perturbgen importable from the clone")

# ── 1. Download pre-trained encoder ──────────────────────────────────────────
print(f"\n=== 1. Download pre-trained encoder ===")

if Path(ENCODER_CKPT).exists() and Path(ENCODER_CKPT).stat().st_size > 1e8:
    print(f"  [SKIP] {ENCODER_CKPT} exists ({Path(ENCODER_CKPT).stat().st_size / 1e9:.2f} GB)")
else:
    # Download from GitHub via the Git LFS media URL
    ENCODER_URL = (
        'https://media.githubusercontent.com/media/Lotfollahi-lab/Perturbgen/main/'
        + ENCODER_REPO_PATH.replace('=', '%3D')
    )
    print(f"  Downloading encoder (~574 MB) from GitHub LFS ...")
    run(['wget', '-q', '--show-progress', '-O', ENCODER_CKPT, ENCODER_URL])
    print(f"  Saved: {ENCODER_CKPT}")

# ── 2. Geneformer vocabulary files ────────────────────────────────────────────
# PerturbGen bundles the 95M Geneformer vocabulary inside perturbgen/pp/
print(f"\n=== 2. Locate Geneformer vocabulary ===")

# These are bundled in the cloned repo; link into VOCAB_DIR for easy access
_pp_dir = Path(REPO_DIR) / 'perturbgen' / 'pp'
ENSEMBL_MAP   = str(_pp_dir / 'ensembl_mapping_dict_gc95M.pkl')
GENE_MEDIANS  = str(_pp_dir / 'gene_median_dict_gftokens_gc95M.pkl')
TOKEN_DICT    = str(_pp_dir / 'token_dict_gftokens_gc95M.pkl')

for p in [ENSEMBL_MAP, GENE_MEDIANS, TOKEN_DICT]:
    assert Path(p).exists(), f"Missing vocabulary file: {p}"
    print(f"  {Path(p).name}  ({Path(p).stat().st_size / 1e6:.1f} MB)")

# Load Ensembl → gene name mapping (used to convert TF names to token IDs)
with open(ENSEMBL_MAP, 'rb') as fh:
    ensembl_map = pickle.load(fh)   # {ensembl_id: gene_symbol}  or reverse

# Determine direction of the mapping
_first_key, _first_val = next(iter(ensembl_map.items()))
if _first_key.startswith('ENSG'):
    ensembl_to_sym = ensembl_map                                 # {ENSG: symbol}
    sym_to_ensembl = {v: k for k, v in ensembl_map.items()}     # {symbol: ENSG}
else:
    sym_to_ensembl = ensembl_map
    ensembl_to_sym = {v: k for k, v in ensembl_map.items()}

print(f"  Vocabulary: {len(sym_to_ensembl):,} human genes with Ensembl IDs")

# ── 3. Preprocess AnnData ─────────────────────────────────────────────────────
print(f"\n=== 3. Preprocess AnnData ===")

PROCESSED_H5AD = f'{OUT}/tokenized/hsc_preprocessed.h5ad'

if Path(PROCESSED_H5AD).exists():
    adata = sc.read_h5ad(PROCESSED_H5AD)
    print(f"  [SKIP] Loaded preprocessed h5ad: {adata.shape}")
else:
    adata = sc.read_h5ad(H5AD)
    print(f"  Raw AnnData: {adata.shape}")

    # 3a. Resolve cluster column
    if CLUSTER_COL not in adata.obs.columns:
        candidates = [c for c in adata.obs.columns if 'res.' in c]
        CLUSTER_COL = candidates[-1] if candidates else 'seurat_clusters'
        print(f"  [WARN] CLUSTER_COL fallback: {CLUSTER_COL!r}")
    adata.obs[CLUSTER_COL] = adata.obs[CLUSTER_COL].astype(str).astype('category')

    # 3b. Mouse→Human ortholog mapping (1:1 strict)
    #
    # Strategy (no network needed — container has no DNS):
    #   1. Load mm_hs_ortho_long exported from the R targets store
    #      (mm_symbol → hs_symbol; columns: Mouse.gene.name, human_ortho)
    #   2. Derive hs_ensembl from the PerturbGen vocab pkl
    #      (ensembl_to_sym is already loaded in step 2 above)
    #
    # Written by the pipeline (target mm_hs_ortho_long_path).

    ORTHO_CACHE   = f'{OUT}/tables/mm_to_hs_ensembl.csv'
    R_ORTHO_CSV   = 'output/tables/mm_hs_ortho_long.csv'

    if Path(ORTHO_CACHE).exists():
        ortho = pd.read_csv(ORTHO_CACHE)
        print(f"  Loaded cached orthologs: {len(ortho):,} pairs")

    elif Path(R_ORTHO_CSV).exists():
        print(f"  Building ortholog table from R export ({R_ORTHO_CSV}) ...")
        raw = pd.read_csv(R_ORTHO_CSV)

        # Normalise column names — R target uses Mouse.gene.name / human_ortho
        col_map = {}
        for c in raw.columns:
            cl = c.lower().replace('.', '_').replace(' ', '_')
            if 'mouse' in cl or 'mm' in cl:
                col_map[c] = 'mm_symbol'
            elif 'human' in cl or 'hs' in cl or 'ortho' in cl:
                col_map[c] = 'hs_symbol'
        raw = raw.rename(columns=col_map)
        assert 'mm_symbol' in raw.columns and 'hs_symbol' in raw.columns, \
            f"Cannot find mm/hs columns in {R_ORTHO_CSV}. Columns: {list(raw.columns)}"

        # Strict 1:1 filter
        raw = (raw[['mm_symbol', 'hs_symbol']]
               .dropna()
               .query("mm_symbol != '' and hs_symbol != ''")
               .drop_duplicates())
        raw = (raw
               .groupby('mm_symbol').filter(lambda x: len(x) == 1)
               .groupby('hs_symbol').filter(lambda x: len(x) == 1)
               .reset_index(drop=True))

        # Derive human Ensembl IDs from PerturbGen vocab pkl (sym_to_ensembl
        # was built in step 2 from ensembl_mapping_dict_gc95M.pkl — no network).
        raw['hs_ensembl'] = raw['hs_symbol'].map(sym_to_ensembl)
        ortho = raw.dropna(subset=['hs_ensembl']).reset_index(drop=True)
        ortho.to_csv(ORTHO_CACHE, index=False)
        print(f"  {len(ortho):,} 1:1 pairs with Ensembl IDs → {ORTHO_CACHE}")

    else:
        raise FileNotFoundError(
            f"Ortholog table not found: {R_ORTHO_CSV}. The pipeline writes it "
            f"(target mm_hs_ortho_long_path): run tools/run-pipeline.sh first."
        )

    mm_to_hs_sym     = dict(zip(ortho['mm_symbol'], ortho['hs_symbol']))
    mm_to_hs_ensembl = dict(zip(ortho['mm_symbol'], ortho['hs_ensembl']))

    # 3c. Subset to genes with valid human Ensembl ID in PerturbGen vocabulary
    keep_mask = np.array([
        g in mm_to_hs_ensembl and mm_to_hs_ensembl[g] in ensembl_to_sym
        for g in adata.var_names
    ])
    adata = adata[:, keep_mask].copy()
    print(f"  After ortholog filter: {adata.shape}")

    # Rename var_names to human Ensembl IDs (required by PerturbGen tokeniser)
    adata.var['mm_symbol']   = adata.var_names.tolist()
    adata.var['hs_symbol']   = [mm_to_hs_sym.get(g, g) for g in adata.var['mm_symbol']]
    adata.var['ensembl_id']  = [mm_to_hs_ensembl[g] for g in adata.var['mm_symbol']]
    adata.var_names = adata.var['ensembl_id'].tolist()
    adata.var_names_make_unique()
    # gene_name must be human symbol — tokeniser will set gene_name=var_names (Ensembl IDs)
    # if we don't pre-populate it, which breaks MT/ribo gene detection.
    adata.var['gene_name'] = adata.var['hs_symbol']

    # 3d. n_counts (required by Geneformer tokeniser)
    # SCT data is log1p-normalized; expm1 to get pseudo-counts, then sum
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    adata.obs['n_counts'] = np.expm1(X).sum(axis=1).astype(np.float32)

    # 3e. HVG selection (on the human-Ensembl-renamed data)
    # Cap N_HVG at the number of available genes (ortholog filter may reduce
    # the gene set below the requested N_HVG — use all genes in that case).
    n_avail = adata.n_vars
    n_hvg_actual = min(N_HVG, n_avail)
    if n_hvg_actual < n_avail:
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg_actual,
                                        flavor='seurat_v3', span=0.5)
        except ModuleNotFoundError:
            # skmisc not installed — fall back to dispersion-based selection
            print("  [WARN] skmisc missing, using flavor='seurat' instead")
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg_actual,
                                        flavor='seurat')
        adata = adata[:, adata.var['highly_variable']].copy()
        print(f"  After HVG ({n_hvg_actual}) filter: {adata.shape}")
    else:
        print(f"  Skipping HVG filter — only {n_avail} genes available (< {N_HVG}), keeping all")
    N_HVG = adata.n_vars   # update so tokenise + config use the actual count

    # 3f. Pseudotime bins → "timepoints" for PerturbGen pairing
    if PSEUDOTIME_COL in adata.obs.columns:
        # CytoTRACE2: high = stem-like. Invert so high pt = differentiated.
        pt = 1.0 - adata.obs[PSEUDOTIME_COL].fillna(0.5).values
        # Tertile bins: 0 = progenitor (source), 1/2 = differentiated (targets)
        bins = pd.qcut(pt, q=3, labels=['0', '1', '2'])
        adata.obs['pseudotime_bin'] = bins.astype(str)
    else:
        # Fallback: use cluster-level pseudotime (rough proxy)
        print("  [WARN] No pseudotime column — assigning uniform bins by cluster rank")
        cl = adata.obs[CLUSTER_COL].astype(int)
        thirds = pd.qcut(cl, q=3, labels=['0', '1', '2'], duplicates='drop')
        adata.obs['pseudotime_bin'] = thirds.fillna('1').astype(str)

    bin_counts = adata.obs['pseudotime_bin'].value_counts()
    print(f"  Pseudotime bins: {dict(bin_counts)}")

    # PerturbGen needs raw (or pseudo-raw) counts in X — use expm1 of SCT.
    # adata.X at this point is still log1p-normalised (SCT); expm1 converts
    # back to pseudo-count scale.  adata already contains only the genes to
    # keep (HVG-filtered or all), so index directly from adata.X.
    X_cur = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    X_out = np.expm1(X_cur).astype(np.float32)
    adata.X = sp.csr_matrix(X_out)

    # Convert nullable string columns (pd.StringDtype) to plain object dtype.
    # pandas ≥ 2.0 returns StringDtype for string columns; anndata's write_h5ad
    # raises RuntimeError unless allow_write_nullable_strings is True.
    # Converting to 'object' avoids the issue in both our write AND in the
    # PerturbGen tokeniser subprocess (which cannot inherit the setting).
    for _df, _name in [(adata.obs, 'obs'), (adata.var, 'var')]:
        for _col in _df.columns:
            if hasattr(_df[_col], 'dtype') and hasattr(_df[_col].dtype, 'name'):
                if 'string' in str(_df[_col].dtype).lower():
                    _df[_col] = _df[_col].astype(object)
    # Also ensure the obs index (cell barcodes) is plain object
    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)

    import anndata as _ad
    _ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(PROCESSED_H5AD)
    # Verify the write was complete — a corrupt write leaves no 'var' group
    with __import__('h5py').File(PROCESSED_H5AD, 'r') as _f:
        assert 'var' in _f, "write_h5ad produced a corrupt file (no 'var' group)"
        assert 'obs' in _f, "write_h5ad produced a corrupt file (no 'obs' group)"
    print(f"  Saved preprocessed h5ad: {PROCESSED_H5AD}")

# ── 4. Tokenise ───────────────────────────────────────────────────────────────
print(f"\n=== 4. Tokenise (perturbgen tokenise) ===")

DATASET_NAME    = 'hsc_perturbgen'

# Actual PerturbGen tokenised output root — derived from perturbgen/configs/paths.py:
#   ROOT = Path(__file__).resolve().parents[3]   ← data/perturbgen
#   TOKENIZED_DIR = ROOT / "T_perturb" / "tokenized_data"
# There is NO --output_dir CLI flag; the tokeniser always writes here.
PERTURBGEN_TOK_DIR = 'data/perturbgen/T_perturb/tokenized_data'
DATASET_SUBDIR     = f'{PERTURBGEN_TOK_DIR}/{DATASET_NAME}'

# Tokenise completion flag (placed after a successful run)
TOK_FLAG = Path(f'{DATASET_SUBDIR}/.tokenised')

if TOK_FLAG.exists():
    print(f"  [SKIP] Tokenised datasets already exist")
else:
    tok_cmd = [
        sys.executable, '-m', 'perturbgen', 'tokenise',
        '--h5ad_path',          PROCESSED_H5AD,
        '--dataset',            DATASET_NAME,
        '--gene_filtering_mode','all',          # HVG selection already done
        '--var_list',           CLUSTER_COL, 'pseudotime_bin',
        '--pairing_mode',       'stratified',
        '--time_obs',           'pseudotime_bin',
        '--main_pairing_obs',   CLUSTER_COL,
        '--nproc',              str(N_WORKERS),
        '--n_hvg',              str(N_HVG),
        '--reference_time',     '0',            # source = most progenitor-like
        '--time_point_order',   '0', '1', '2',
        '--gene_median_path',   GENE_MEDIANS,
        '--token_dict_path',    TOKEN_DICT,
        '--gene_mapping_path',  ENSEMBL_MAP,
        # NOTE: --output_dir is NOT a valid flag for perturbgen tokenise.
        # Output always goes to TOKENIZED_DIR/{dataset}/ (see perturbgen/configs/paths.py).
    ]
    run(tok_cmd, cwd=REPO_DIR)
    os.makedirs(DATASET_SUBDIR, exist_ok=True)
    TOK_FLAG.touch()
    print(f"  Tokenisation complete → {DATASET_SUBDIR}")

# Tokenised output paths (from perturbgen/configs/paths.py + GF_tokenisation.py):
# All paths are derived from TOKENIZED_DIR/{dataset}/ with suffix 'all' for --gene_filtering_mode=all.
#
#   {DATASET_SUBDIR}/
#     token_id_to_genename_all.pkl      — token_id → gene_name mapping
#     tokenid_to_rowid_all.pkl          — token_id → row_id mapping
#     h5ad_pairing_all_GF_genes/        — adata subset to GF vocabulary genes
#     h5ad_pairing_all_src/             — src adata per reference time point
#     h5ad_pairing_all_tgt/             — tgt adata per time point
#     dataset_all_src/                  — tokenised src HuggingFace datasets
#     dataset_all_tgt/                  — tokenised tgt HuggingFace datasets
#
MAPPING_PKL      = f'{DATASET_SUBDIR}/token_id_to_genename_all.pkl'
ROWID_PKL        = f'{DATASET_SUBDIR}/tokenid_to_rowid_all.pkl'
TGT_FOLDER       = f'{DATASET_SUBDIR}/dataset_all_tgt'   # folder; read_dataset_files lists *.dataset inside
TGT_ADATA_FOLDER = f'{DATASET_SUBDIR}/h5ad_pairing_all_tgt'  # folder; read_dataset_files lists *.h5ad inside

# load_from_disk needs the actual HF Dataset dir (contains dataset_info.json), not the parent folder.
# sc.read_h5ad needs the actual .h5ad file path, not the parent folder.
_src_ds_candidates = sorted(Path(f'{DATASET_SUBDIR}/dataset_all_src').glob('*.dataset'))
_src_h5_candidates = sorted(Path(f'{DATASET_SUBDIR}/h5ad_pairing_all_src').glob('*.h5ad'))
if not _src_ds_candidates:
    raise FileNotFoundError(f"No *.dataset found in {DATASET_SUBDIR}/dataset_all_src — re-run tokenise step")
if not _src_h5_candidates:
    raise FileNotFoundError(f"No *.h5ad found in {DATASET_SUBDIR}/h5ad_pairing_all_src — re-run tokenise step")
SRC_DATASET = str(_src_ds_candidates[0])
SRC_ADATA   = str(_src_h5_candidates[0])

print(f"  src dataset:  {SRC_DATASET}")
print(f"  tgt folder:   {TGT_FOLDER}")
print(f"  mapping pkl:  {MAPPING_PKL}")

# ── 5. Train count decoder ────────────────────────────────────────────────────
print(f"\n=== 5. Train count decoder ===")

DECODER_CKPT_DIR = f'{CKPT_DIR}/count_decoder'
DECODER_FLAG     = Path(f'{DECODER_CKPT_DIR}/.trained')

if DECODER_FLAG.exists():
    print(f"  [SKIP] Count decoder already trained")
else:
    decoder_cmd = [
        sys.executable, '-m', 'perturbgen', 'train-decoder',
        '--train_mode',          'count',
        '--split',               'False',
        '--splitting_mode',      'stratified',
        '--split_obs',           CLUSTER_COL,
        '--encoder',             'scmaskgit',
        '--encoder_path',        ENCODER_CKPT,   # scmaskgit weights; overrides HPC-hardcoded default
        '--output_dir',          DECODER_CKPT_DIR,
        '--pred_tps',            '1', '2',       # our tgt has only 1_1.dataset and 2_2.dataset
        '--var_list',            'pseudotime_bin',  # metadata cols present in our datasets; cluster col key differs by run so omitted
        '--src_dataset',         SRC_DATASET,
        '--tgt_dataset_folder',  TGT_FOLDER,
        '--src_adata',           SRC_ADATA,
        '--tgt_adata_folder',    TGT_ADATA_FOLDER,
        '--mapping_dict_path',   MAPPING_PKL,
        '--batch_size',          str(BATCH_SIZE),
        '--epochs',              str(EPOCHS_DEC),
        '--count_lr',            '5e-5',
        '--cellgen_lr',          '1e-5',
        '--n_workers',           str(N_WORKERS),
    ]
    run(decoder_cmd, cwd=REPO_DIR)
    DECODER_FLAG.touch()
    print(f"  Count decoder trained: {DECODER_CKPT_DIR}")

# Find decoder checkpoint
decoder_ckpts = sorted(Path(DECODER_CKPT_DIR).rglob('*.ckpt'))
DECODER_CKPT  = str(decoder_ckpts[-1]) if decoder_ckpts else ''
print(f"  Using decoder ckpt: {DECODER_CKPT}")

# ── 6. Select top TFs and map to Ensembl IDs ─────────────────────────────────
print(f"\n=== 6. Select top {TOP_N_TFS} TFs from SCENIC ===")

auc_raw = pd.read_csv(AUC_MTX, index_col=0)
auc_tfs = {col: re.sub(r'\([+-]\)$', '', col) for col in auc_raw.columns}
mean_auc = auc_raw.mean(axis=0).rename(index=auc_tfs).groupby(level=0).max()

# Load ortho table
ortho = pd.read_csv(f'{OUT}/tables/mm_to_hs_ensembl.csv')
mm_to_hs_ensembl = dict(zip(ortho['mm_symbol'], ortho['hs_ensembl']))
mm_to_hs_sym     = dict(zip(ortho['mm_symbol'], ortho['hs_symbol']))

# Load token dictionary to verify TF tokens exist
with open(TOKEN_DICT, 'rb') as fh:
    token_dict = pickle.load(fh)  # {gene_symbol: token_id} or {ensembl_id: token_id}

# Determine token dict key type
_tk, _ = next(iter(token_dict.items()))
token_key_is_ensembl = _tk.startswith('ENSG')

def get_token_id(mm_sym):
    """Get PerturbGen token ID for a mouse gene symbol."""
    hs_ensembl = mm_to_hs_ensembl.get(mm_sym)
    if hs_ensembl is None:
        return None
    key = hs_ensembl if token_key_is_ensembl else mm_to_hs_sym.get(mm_sym)
    return token_dict.get(key)

top_tfs = []
for tf in mean_auc.sort_values(ascending=False).index:
    hs_ens = mm_to_hs_ensembl.get(tf)
    if hs_ens and get_token_id(tf) is not None:
        top_tfs.append(tf)
    if len(top_tfs) >= TOP_N_TFS:
        break

print(f"  TFs selected ({len(top_tfs)}):")
for tf in top_tfs:
    print(f"    {tf:12s}  → {mm_to_hs_sym.get(tf,'?'):12s}  {mm_to_hs_ensembl.get(tf,'?')}")

# Load processed adata for analysis
adata_proc = sc.read_h5ad(PROCESSED_H5AD)
if CLUSTER_COL not in adata_proc.obs.columns:
    candidates = [c for c in adata_proc.obs.columns if 'res.' in c]
    CLUSTER_COL = candidates[-1] if candidates else 'seurat_clusters'
clusters    = adata_proc.obs[CLUSTER_COL].astype(str).unique().tolist()
cluster_lab = adata_proc.obs[CLUSTER_COL].astype(str).values

# Pseudotime gradient in UMAP space for inner product analysis
pt_gradient = None
if PSEUDOTIME_COL in adata_proc.obs.columns and 'X_umap' in adata_proc.obsm:
    # Reload original h5ad for UMAP coords (preprocessing discards original obsm)
    adata_orig = sc.read_h5ad(H5AD)
    if 'X_umap' in adata_orig.obsm:
        shared_idx = adata_orig.obs_names.isin(adata_proc.obs_names)
        pt = 1.0 - adata_orig.obs.loc[shared_idx, PSEUDOTIME_COL].fillna(0.5).values
        umap = adata_orig.obsm['X_umap'][shared_idx]
        lr = LinearRegression(fit_intercept=True).fit(pt.reshape(-1, 1), umap)
        pt_gradient = lr.coef_[0] / np.linalg.norm(lr.coef_[0])
        adata_umap = adata_orig[shared_idx].copy()
        del adata_orig
    print(f"  Pseudotime gradient: {pt_gradient}")

# ── 7. Perturbation inference (one config per TF) ─────────────────────────────
print(f"\n=== 7. Perturbation inference ===")

# Number of target timepoints
TGT_TIMEPOINTS = [1, 2]   # corresponds to pseudotime bins 1 and 2
PERT_SEQUENCE  = ['src']  # perturb source cells, predict effect at targets

def make_pert_config(tf_mm, ens_id, cfg_path):
    """Write a perturbation YAML config for one TF."""
    cfg = {
        'data': {
            'src_dataset_file':  SRC_DATASET,
            'tgt_dataset_folder': TGT_FOLDER,
            'src_adata':         SRC_ADATA,
            'tgt_adata_folder':  TGT_ADATA_FOLDER,
        },
        'trainer': {
            'use_count_decoder':     True,
            'n_samples':             3,
            'tgt_vocab_size':        N_HVG + 2,   # +2 for special tokens
            'd_model':               768,
            'generate':              False,
            'use_positional_encoding': False,
            'dropout':               0,
            'num_heads':             8,
            'num_layers':            6,
            'loss_mode':             'zinb',
            'd_ff':                  32,
            'max_seq_length':        N_HVG,
            'mask_scheduler':        'pow',
            'output_dir':            PERT_DIR,
            'mapping_dict_path':     MAPPING_PKL,
            'tokenid_to_rowid_path': ROWID_PKL,
            'encoder':               'scmaskgit',
            'encoder_path':          ENCODER_CKPT,
            'context_mode':          False,
            'temperature':           1.5,
            'iterations':            19,
            'pos_encoding_mode':     'time_pos_sin',
            'exclude_src':           False,
            'perturbation_mode':     'mask',
            'genes_to_perturb':      [ens_id],
            'validation_mode':       'inference',
            'perturbation_sequence': PERT_SEQUENCE,
            'pred_tps':              TGT_TIMEPOINTS,
            'ckpt_masking_path':     DECODER_CKPT,
        },
        'datamodule': {
            'batch_size':   BATCH_SIZE,
            'num_workers':  N_WORKERS,
            'shuffle':      False,
            'split':        False,
            'max_len':      N_HVG,
            'var_list':     [CLUSTER_COL, 'pseudotime_bin'],
            'pred_tps':     TGT_TIMEPOINTS,
        },
        'model': {
            'precision':         32,
            'ckpt_masking_path': DECODER_CKPT,
        },
    }
    with open(cfg_path, 'w') as fh:
        yaml.dump(cfg, fh, default_flow_style=False)

# Storage for post-hoc analysis
cluster_effects = {}    # TF → {cluster: mean_effect}
pt_inner_prod   = {}    # TF → scalar

val_script = str(Path(REPO_DIR) / 'perturbgen' / 'Perturb' / 'val.py')

for i, tf_mm in enumerate(top_tfs):
    ens_id  = mm_to_hs_ensembl[tf_mm]
    hs_sym  = mm_to_hs_sym[tf_mm]
    cfg_path = f'{CFG_DIR}/perturbation_{tf_mm}.yaml'

    print(f"\n  [{i+1:2d}/{len(top_tfs)}] {tf_mm} (→{hs_sym}, {ens_id})")

    # Expected output pattern from val.py
    result_pattern = f'adata_g{ens_id}_ssrc*.h5ad'
    existing = list(Path(PERT_DIR).glob(result_pattern))
    if existing:
        print(f"    [SKIP] Result exists: {existing[0].name}")
        adata_pert = sc.read_h5ad(str(existing[0]))
    else:
        make_pert_config(tf_mm, ens_id, cfg_path)
        try:
            run([sys.executable, val_script, '--config', cfg_path], cwd=REPO_DIR)
            new_results = list(Path(PERT_DIR).glob(result_pattern))
            if not new_results:
                print(f"    [WARN] No output h5ad found for {tf_mm}")
                cluster_effects[tf_mm] = {}
                pt_inner_prod[tf_mm] = np.nan
                continue
            adata_pert = sc.read_h5ad(str(new_results[-1]))
        except subprocess.CalledProcessError as exc:
            print(f"    [FAIL] {exc}")
            cluster_effects[tf_mm] = {}
            pt_inner_prod[tf_mm] = np.nan
            continue

    # ── Post-perturbation analysis ────────────────────────────────────────
    # PerturbGen outputs predicted expression at target timepoints.
    # Compare perturbed predictions against unperturbed predictions (ctrl)
    # by running DEG analysis within each cluster.

    # adata_pert.obs should have 'status' (perturbed/unperturbed) and cluster col
    if 'status' not in adata_pert.obs.columns:
        # Mark half as perturbed, half as ctrl based on convention
        n = len(adata_pert)
        adata_pert.obs['status'] = (
            ['perturbed'] * (n // 2) + ['unperturbed'] * (n - n // 2)
        )

    # Normalize for DE
    sc.pp.normalize_total(adata_pert, target_sum=1e4)
    sc.pp.log1p(adata_pert)

    cluster_deg_counts = {}
    clust_obs = CLUSTER_COL if CLUSTER_COL in adata_pert.obs.columns else None

    if clust_obs:
        for cl in adata_pert.obs[clust_obs].astype(str).unique():
            sub = adata_pert[adata_pert.obs[clust_obs].astype(str) == cl].copy()
            if len(sub['status'].unique()) < 2:
                cluster_deg_counts[cl] = 0
                continue
            try:
                sc.tl.rank_genes_groups(sub, groupby='status',
                                        reference='unperturbed',
                                        method='wilcoxon',
                                        n_genes=sub.n_vars)
                degs = sc.get.rank_genes_groups_df(sub, group='perturbed')
                n_sig = int((degs['pvals_adj'] < 0.05).sum())
                cluster_deg_counts[cl] = n_sig
            except Exception:
                cluster_deg_counts[cl] = 0
    else:
        # No cluster obs in output: use overall DEG count
        try:
            sc.tl.rank_genes_groups(adata_pert, groupby='status',
                                    reference='unperturbed', method='wilcoxon')
            degs = sc.get.rank_genes_groups_df(adata_pert, group='perturbed')
            n_sig = int((degs['pvals_adj'] < 0.05).sum())
            cluster_deg_counts = {'all': n_sig}
        except Exception:
            cluster_deg_counts = {'all': 0}

    cluster_effects[tf_mm] = cluster_deg_counts
    print(f"    DEG counts per cluster: {cluster_deg_counts}")

    # Pseudotime inner product — use mean UMAP shift of perturbed vs ctrl
    if pt_gradient is not None and 'X_umap' in adata_pert.obsm:
        try:
            pert_mask = adata_pert.obs['status'] == 'perturbed'
            ctrl_mask = adata_pert.obs['status'] == 'unperturbed'
            if pert_mask.sum() > 0 and ctrl_mask.sum() > 0:
                mean_pert = adata_pert.obsm['X_umap'][pert_mask].mean(axis=0)
                mean_ctrl = adata_pert.obsm['X_umap'][ctrl_mask].mean(axis=0)
                delta = mean_pert - mean_ctrl
                pt_inner_prod[tf_mm] = float(delta @ pt_gradient)
        except Exception:
            pt_inner_prod[tf_mm] = np.nan

    # Per-TF UMAP of DEG count per cluster
    if clust_obs and 'adata_umap' in dir() and cluster_deg_counts:
        try:
            umap_df = adata_umap.obs[[CLUSTER_COL]].copy()
            umap_df['deg_count'] = umap_df[CLUSTER_COL].astype(str).map(cluster_deg_counts).fillna(0)
            fig, ax = plt.subplots(figsize=(7, 6))
            sc.pl.umap(adata_umap, color=CLUSTER_COL, ax=ax, show=False,
                       title=f'{tf_mm} KO — clusters', legend_loc='right margin')
            plt.tight_layout()
            save_fig(f'01_perturbation_{tf_mm}.png')
        except Exception:
            plt.close('all')

# ── 8. Save summary tables ────────────────────────────────────────────────────
print(f"\n=== 8. Save tables ===")

# TF × cluster DEG count
cl_df = pd.DataFrame(cluster_effects).T   # TFs × clusters
cl_df.index.name = 'TF'
cl_df.to_csv(f'{TABLE_DIR}/cluster_effect_summary.csv')
print(f"  cluster_effect_summary.csv  {cl_df.shape}")

if pt_inner_prod:
    ip_df = (pd.Series(pt_inner_prod, name='pt_inner_product')
             .reset_index().rename(columns={'index': 'TF'})
             .dropna()
             .sort_values('pt_inner_product'))
    ip_df.to_csv(f'{TABLE_DIR}/pseudotime_inner_products.csv', index=False)
    print(f"  pseudotime_inner_products.csv  {ip_df.shape}")

# ── 9. Summary heatmap ────────────────────────────────────────────────────────
print(f"\n=== 9. Summary heatmap ===")

valid_tfs = [tf for tf in top_tfs if tf in cl_df.index and cl_df.loc[tf].notna().any()]
if valid_tfs and cl_df.shape[1] > 0:
    hm_data = cl_df.loc[valid_tfs].fillna(0).astype(float)
    # Z-score across clusters per TF (how does the TF affect clusters differently?)
    hm_z = hm_data.apply(lambda row: zscore(row, nan_policy='omit'), axis=1).fillna(0)
    col_order = hm_z.mean(axis=0).sort_values(ascending=False).index.tolist()
    row_order = hm_z.max(axis=1).sort_values(ascending=False).index.tolist()
    hm_plot = hm_z.loc[row_order, col_order]

    fig_h = max(5, len(valid_tfs) * 0.45)
    fig_w = max(8, hm_plot.shape[1] * 0.65)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        hm_plot, ax=ax, cmap='YlOrRd', vmin=0, vmax=3,
        linewidths=0.4, linecolor='white',
        cbar_kws=dict(label='Z-score (significant DEGs per cluster)', shrink=0.6),
    )
    ax.set_xlabel('Cluster', fontsize=11)
    ax.set_ylabel('TF (KO)', fontsize=11)
    ax.set_title('PerturbGen TF KO: DEGs per cluster\n'
                 '(z-scored count of FDR<0.05 genes)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    save_fig('02_summary_heatmap.png')

# ── 10. Pseudotime inner product bar ─────────────────────────────────────────
if pt_inner_prod and any(np.isfinite(v) for v in pt_inner_prod.values()):
    ip_plot = pd.Series(pt_inner_prod).dropna().sort_values()
    colours = ['#b2182b' if v > 0 else '#2166ac' for v in ip_plot.values]
    fig, ax = plt.subplots(figsize=(8, max(4, len(ip_plot) * 0.45)))
    ax.barh(ip_plot.index, ip_plot.values, color=colours, edgecolor='white')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel('UMAP shift along differentiation axis\n'
                  '(+) pro-differentiation  |  (−) blocks differentiation', fontsize=10)
    ax.set_title('PerturbGen TF KO: effect on differentiation direction',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    save_fig('03_pseudotime_inner_product.png')

# ── 11. Cross-tool comparison with CellOracle ─────────────────────────────────
if CELLORACLE_TABLE.exists() and valid_tfs and cl_df.shape[1] > 0:
    print(f"\n=== 11. Cross-tool comparison ===")
    co_df = pd.read_csv(str(CELLORACLE_TABLE), index_col=0)

    shared_tfs  = [tf for tf in valid_tfs if tf in co_df.index]
    shared_clus = [c for c in col_order if c in co_df.columns]

    if shared_tfs and shared_clus:
        co_z = co_df.loc[shared_tfs, shared_clus]
        co_z = co_z.apply(lambda row: zscore(row, nan_policy='omit'), axis=1).fillna(0)
        pg_z = hm_z.loc[shared_tfs, shared_clus]

        fig, axes = plt.subplots(1, 2, figsize=(max(10, len(shared_clus)) * 1.5,
                                                 max(5, len(shared_tfs) * 0.45)))
        for ax, (data, title) in zip(axes, [
            (co_z, 'CellOracle\n(GRN linear propagation)'),
            (pg_z, 'PerturbGen\n(generative transformer)'),
        ]):
            sns.heatmap(data, ax=ax, cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                        linewidths=0.3, linecolor='white',
                        cbar_kws=dict(label='Z-score', shrink=0.7))
            ax.set_title(title, fontsize=11, fontweight='bold')
            ax.set_xlabel('Cluster', fontsize=9)
            ax.set_ylabel('TF (KO)', fontsize=9)
            ax.tick_params(axis='y', labelsize=7)
        plt.suptitle('TF KO perturbation: CellOracle vs PerturbGen',
                     fontsize=13, y=1.02)
        plt.tight_layout()
        save_fig('04_comparison_celloracle_vs_perturbgen.png')
        print(f"  Comparison: {len(shared_tfs)} TFs × {len(shared_clus)} clusters")

# ── Done ─────────────────────────────────────────────────────────────────────
print(f"\n=== Done ===")
print(f"Outputs in {OUT}:")
for d in [FIG_DIR, TABLE_DIR]:
    for f in sorted(Path(d).iterdir()):
        print(f"  {str(Path(f).relative_to(OUT)):60s}  {f.stat().st_size/1e3:8.1f} KB")
