# DEVLOG — mpn_10x_camelia

---

## 2026-05-11 — Human reference processing (Suo/Jardine), Moran's I targets, container fixes

### What changed
- `data_hub/suo_2022/convert_h5ad_to_newformat.py` (new): pure h5py in-place fix for `nullable-string-array` groups — detects h5ad obs/var groups with `values`+`mask` children, extracts values, deletes group, creates flat `string-array` dataset; makes old-format h5ad files readable by `anndataR`.
- `data_hub/jardine_2021/convert_h5ad_to_newformat.py` (new): identical conversion logic for Jardine 2021 files (fig1b, prog, mnc).
- `data_hub/suo_2022/process_suo_refs.R` (new): loads converted Suo 2022 h5ad files via `anndataR`, builds Seurat reference objects (LogNormalize → FindVariableFeatures → ScaleData → PCA → Harmony for HSC_PROGENITORS only → UMAP with `return.model = TRUE`); saves RDS per dataset; MEM and MYELOID datasets have Harmony disabled; HSC_IMMUNE and STROMA commented out.
- `data_hub/jardine_2021/process_jardine_refs.R` (new): same processing pipeline for Jardine 2021 references; no Harmony for any dataset.
- Both processing R scripts updated: layer normalisation block now drops all `as_Seurat()` layers, converts the raw counts matrix from dense to sparse (`as(mat, "CsparseMatrix")`), assigns as `counts`, and drops `scale.data` before `saveRDS` — reduces RDS file size substantially.
- `_targets.R`:
  - Added `#jardine` section: three `tarchetypes::tar_map` blocks (jardine_fig1b, jardine_prog, jardine_mnc) × three compartments (stroma/hsc/mnc) = 18 new projection targets (`proj_jardine_21_hs_*` + `proj_path_jardine_21_hs_*`); mirrors existing `#suo` section format.
  - Added `moransi_hsc` `tar_map`: `Seurat::RunMoransI` per condition (`neg`/`pos`/`wt`) using `umap.nn` embeddings as spatial coordinates; returns gene × (morans_i, p_value, p_adj BH) data frame + CSV file target.
  - Added `moransi_pca_hsc` `tar_map`: identical structure but uses `pca.SCT[, 1:10]` as spatial coordinates instead of UMAP — captures higher-dimensional neighbourhood signal.
- `src/scrna-seq/scrna-seq_utilities.R`:
  - Added `moran_sparse()`: computes Moran's I for a single gene vector against a sparse weight matrix (`dgCMatrix`) without allocating a dense n² distance matrix; uses `Matrix::crossprod` throughout; suitable for ~50k-cell graphs.
  - Added `nn_to_sym_dist()`: converts a Seurat `Neighbor` object (asymmetric KNN) to a symmetric sparse distance matrix; three symmetrisation methods (`average`, `max`, `min`).
- `container/bioinfo.def`:
  - Added `libudunits2-dev libgdal-dev libgeos-dev libproj-dev` to system packages (required by `sf` → `spdep`).
  - Removed `celloracle` from Python pip install — pins `pandas<=1.5.3` and `gimmemotifs<=0.17.1`, both incompatible with Python 3.12 (`distutils` removed, `configparser.SafeConfigParser` removed).
- `container/renv_setup.R`:
  - Added `Rfast2` to CRAN packages (required by `Seurat::RunMoransI` for fast Moran's I computation).

### Why
Two parallel workstreams this session:

1. **Human reference processing**: Suo 2022 (fetal haematopoiesis) and Jardine 2021 (fetal BM) are key human references for label transfer onto the MPN HSC compartment. They are distributed as old-format h5ad files (anndata <0.9) which use `nullable-string-array` encoding for string metadata — a GROUP structure in HDF5 that `anndataR` cannot read (expects a flat DATASET). A pure h5py in-place rewrite was implemented to fix the encoding without re-writing the full file through anndata, which would re-introduce the problem via pandas type inference.

2. **Moran's I for spatially variable genes**: Identifying genes whose expression is locally autocorrelated within each condition (mCh+, mCh−, WT) provides complementary evidence to DE analysis — it highlights genes with coherent expression patches rather than bulk shifts. Two spatial coordinate systems are compared: UMAP (2D, captures global topology) and PCA 1–10 (higher-dimensional, captures more variance and avoids UMAP discontinuities).

### Bugs / issues
- **`anndataR` warns "No counts or data layer"**: `as_Seurat()` maps h5ad `X` to a layer but may name it `"X"` rather than `"counts"`. Old scripts did `LayerData(srt, "RNA", "counts") <- LayerData(srt, "RNA", "X")` but forgot to delete `"X"`, leaving both layers plus `"data"` and `"scale.data"` — all as dense matrices. Fixed: drop all layers first, convert the raw matrix to sparse, assign as `"counts"` only; drop `scale.data` before saving.
- **`configparser.SafeConfigParser` removed in Python 3.12**: `gimmemotifs 0.17.0` (pulled as `celloracle` dependency) fails at metadata generation stage. Attempted fix (install from GitHub first) was insufficient because `celloracle` pins `gimmemotifs<=0.17.1` and pip's resolver always downloads from PyPI. Root cause: `celloracle` is frozen at Python 3.10/3.11 era. Resolution: remove `celloracle` from container; install separately in a conda/mamba env with Python 3.10 if needed.
- **`pandas<=1.5.3` incompatible with Python 3.12**: `celloracle` forces a downgrade to pandas 1.5.3 which requires `distutils` (removed in 3.12). Confirmed unfixable without significant patching; removal of `celloracle` resolves the build failure.
- **`sf` fails to install without system libs**: `spdep` depends on `sf` which needs `libudunits2-dev`, `libgdal-dev`, `libgeos-dev`, `libproj-dev` — none were in the container. Added to `bioinfo.def` system packages block.
- **`Rfast2::moranI` integer overflow on large n**: `Rfast2::moranI` uses `seq.int` internally over the full n×n weights matrix; `n²` overflows `.Machine$integer.max` (~2.1B) when n > ~46k cells. Per-condition subsets (~15–20k cells each) fit within the limit, so `Seurat::RunMoransI` is called on subsets rather than the full object.
- **Trailing comma error in R list with commented entries**: `list(a, b, #,c)` causes "argument 4 is empty". Fixed by moving comma to start of commented line: `#,list(...)`.
- **RunHarmony "number of columns must match"**: NA values in grouping vars (`method`, `donor`, `sex`) cause Harmony to fail. Fixed by replacing NAs with `"unknown"` before `RunHarmony`.

### Decisions
- **No SCTransform for human references**: Label transfer via `project_mm_onto_hs_seurat()` accesses `ref@assays$RNA@counts` directly — using LogNormalize + PCA keeps the reference compatible and avoids SCTransform's Pearson residual scale (meaningless for cross-species projection).
- **No Harmony for Jardine datasets**: Jardine references are already batch-corrected in the published h5ad; re-running Harmony would destroy the corrected embedding.
- **Harmony only for Suo HSC_PROGENITORS**: MEM and MYELOID subsets have sufficient within-condition coherence that Harmony integration is not needed; HSC_PROGENITORS spans multiple donors/methods requiring correction.
- **`celloracle` excluded from container**: Fundamental Python 3.12 incompatibility. Installing in the container would require downgrading pandas globally, breaking scVelo/CellRank which require pandas ≥ 2.0. Decision: separate conda env.
- **Moran's I with UMAP and PCA coordinates in parallel**: Using both gives complementary spatial scales — UMAP preserves cluster topology, PCA 1–10 captures continuous variance without the non-linear distortions of UMAP.
- **`p.adjust(method = "BH")` per condition separately**: Each condition's Moran's I results are adjusted within that condition — the multiple-testing burden is the number of variable features tested, not all genes.

### Compartment status
- **Stroma** ✅: pipeline unchanged; Jardine projection targets added (stroma compartment included in `#jardine` tar_map).
- **HSC** ✅: Moran's I targets added (UMAP + PCA coordinate versions); Jardine + Suo projection targets in pipeline; `mpn_processed_hsc_fullcc` used as source throughout.
- **MNC** ✅: pipeline unchanged; Jardine MNC projection targets added.
- **HSC+MNC** ⏳: not touched this session.

### Next steps
- [ ] Delete stale Jardine RDS files and re-run `process_jardine_refs.R` (layer fix)
- [ ] Run `process_suo_refs.R` for MEM and MYELOID (HSC_PROGENITORS already processed)
- [ ] Rebuild container (`./scripts/build.sh --force`) to include `libudunits2-dev` + `Rfast2` + without `celloracle`
- [ ] Run Jardine projection tar_make targets
- [ ] Run `moransi_hsc_*` and `moransi_pca_hsc_*` targets; compare UMAP vs PCA rankings
- [ ] Install `spdep` in container (now that system libs are added) for p-value-aware sparse Moran's I

---

## 2026-04-25 — SuperCell metacell refinement + Milo per-cell scores + scatter-pie UMAP

### What changed
- `_targets.R`:
  - `supercells_hsc`: `gamma` raised from 50 → 50 (user reverted to 50 with `split_by = "sample.name"` and `fast.merge = TRUE`); `split_by` and `fast.merge` parameters now explicitly passed
  - Added `milo_hsc_cell_scores` target: aggregates neighbourhood-level logFC back to individual cells (mean across all nhoods each cell belongs to); returns named list of three data frames — one per contrast (`mChNeg_vs_WT`, `mChPos_vs_WT`, `mChPos_vs_mChNeg`)
- `src/scrna-seq/scrna-seq_utilities.R`:
  - `compute_supercells()`: added `split_by` (default `"sample.name"`) and `fast.merge` (default `FALSE`) parameters; `split_by = NULL` disables per-sample splitting for trajectory/global use; exposes `fast.merge` to `SCimplify_from_embedding`
- `notebooks/compositional-analysis.qmd`:
  - Added `plot_medoid_pie()` helper to the helpers chunk: places one `scatterpie` pie per group at the PCA-space medoid cell's UMAP coordinates; pie fill = proportion of a metadata variable; pie area ∝ sqrt(n_cells); optional background cell cloud
  - User added example usage in the HSC purity section (`res_20` clusters × `sample.name` fill)

### Why
Three related issues with the SuperCell and Milo workflows were addressed:

1. **SuperCell metacell count too high**: `gamma = 50` with `cell.split.condition = sample.name` produced 3110 supercells instead of the expected ~1074 because gamma is applied per sample group. The effective graining was ~17x, not 50x. Options discussed: increase gamma, coarsen the split condition, or disable splitting.

2. **SuperCell size distribution right-skewed**: Small per-sample groups force walktrap to produce many tiny communities that can't fill a full gamma-sized supercell. `fast.merge = TRUE` post-hoc merges communities smaller than gamma/2. Exposing `split_by` allows switching from per-sample to per-condition grouping or no grouping.

3. **Milo per-cell scores**: Milo assigns logFC to neighbourhoods, not cells. To visualise DA on a per-cell UMAP, neighbourhood logFC must be aggregated back. Mean across all neighbourhoods a cell belongs to (weighted equally) is the simplest and most interpretable approach. Cells not covered by any neighbourhood (n_nhoods == 0) receive NA.

### Bugs / issues
- CRAN mirror not set in fresh R environments: `install.packages()` fails with `trying to use CRAN without setting a mirror`. Fix: pass `repos = 'https://cloud.r-project.org'` explicitly, or add `options(repos = c(CRAN = "https://cloud.r-project.org"))` to `~/.Rprofile`.
- SuperCell effective gamma lower than set value: root cause is `cell.split.condition` applying gamma per group. With many small samples, effective cells-per-supercell = n_cells_per_sample / gamma, not n_cells_total / gamma.

### Decisions
- **`fast.merge` default kept `FALSE`** in `compute_supercells()`: it changes the membership structure and should be opt-in. The `_targets.R` target passes it explicitly as `TRUE`.
- **`split_by = NULL` for trajectory use**: Without sample splitting, supercells span conditions — acceptable for CellRank/velocity inputs but breaks pseudobulk DE. Decision: keep `split_by = "sample.name"` for the pipeline target; use `NULL` only for trajectory-specific targets.
- **Medoid approximation**: True medoid (minimising sum of pairwise distances) is O(n²) per group and impractical for thousands of cells. Cell nearest to the group centroid in PCA space is used as an approximation — standard practice and accurate enough for UMAP placement.
- **`scatterpie` for group composition UMAP**: Preferred over `ggforce::geom_arc_bar` for simplicity; radius scaled by sqrt(n_cells) so pie area is proportional to group size.
- **Milo per-cell mean logFC**: Cells at cluster boundaries with mixed-sign neighbourhood logFC will average toward 0 — a limitation to communicate in figure legends. Not used for inference, only visualisation.

### Compartment status
- **Stroma** ✅: pipeline complete; Milo stroma target unchanged.
- **HSC** ✅: `milo_hsc_cell_scores` added; `supercells_hsc` updated with new parameters; scatter-pie example in compositional notebook.
- **MNC** ✅: pipeline unchanged this session.
- **HSC+MNC** ⏳: not touched this session.

### Next steps
- [ ] Run `tar_make()` to validate `milo_hsc_cell_scores` and `supercells_hsc` with new parameters
- [ ] Plot `milo_hsc_cell_scores` per-cell logFC on HSC UMAP (three panels, one per contrast)
- [ ] Check `hist(supercells_hsc$supercell_size)` to confirm `fast.merge` reduced right-skew
- [ ] Evaluate whether `split_by = NULL` + higher gamma is needed for CellRank input supercells
- [ ] Add `scatterpie` to `container/renv_setup.R` if not present

---

## 2026-04-21 19:00 — pySCENIC cisTarget debugging: pandas 3.0 + numba compatibility fixes

### What changed
- `cache/python_pkgs/pyscenic/utils.py`:
  - Line 159: `.values` → `.values.copy()` — numba JIT `masked_rho4pairs` received a read-only array from newer pandas; `.copy()` forces a writable allocation
  - Lines 244–246: `groupby(target).apply(lambda grp: grp.nlargest(...))` → `sort_values + groupby.head(n)` — pandas 3.0 excludes the grouping column from `apply()` result by default (`include_groups` breaking change), causing `KeyError: 'target'` in `modules4top_factors`; `sort_values + head` avoids `apply` entirely and preserves all columns
- `notebooks/scenic_run.py`:
  - PYSCENIC binary detection changed from hardcoded `/project/cache/python_pkgs/bin/pyscenic` to dynamic: check bin path → `shutil.which` → Python entry-point fallback
  - Entry-point fallback now includes `-B` flag: `[sys.executable, '-B', '-c', 'from pyscenic.cli.pyscenic import main; main()']` — prevents stale `.pyc` in `cache/python_pkgs/pyscenic/` from masking source edits
- `scripts/run_scenic.sh`:
  - Added `-B` to outer Python invocation: `python3 -B ${SCRIPT}` — prevents `scenic_run.py` bytecode from being cached or read from stale cache
- `scripts/run_targets.sh`: minor cleanup (user-edited, already in diff)

### Why
The pySCENIC cisTarget step (step 3 of 4) failed repeatedly on the first real run after moving to the macOS Lima + Apptainer environment. The container uses Python 3.12 with `cache/python_pkgs` installed packages (pandas 3.0.2, numba). Two independent compatibility breaks prevented `pyscenic ctx` from completing:

1. **numba readonly array**: pandas 3.0 returns read-only views from `.values` on certain DataFrames. The numba-JIT function `masked_rho4pairs` was compiled expecting a writable `array(float64, 2d, C)` — it has no signature that accepts `readonly array(float64, 2d, C)`.
2. **pandas 3.0 `groupby.apply` column exclusion**: In pandas 3.0, `DataFrameGroupBy.apply()` no longer includes the groupby key column in the result DataFrame by default (`include_groups=False` is now the default). `modules4top_factors` grouped by `target`, then tried to access `df_grp['target']` from the result — which was no longer present.

Additionally, stale `.pyc` files in `cache/python_pkgs/pyscenic/` caused repeated confusion where source edits were not reflected in container runs — fixed by adding `-B` to both the outer Python call and the subprocess entry-point invocation.

### Bugs / issues
- **`FileNotFoundError: pyscenic`**: `bin/pyscenic` script doesn't exist in `cache/python_pkgs/bin/` (only dask/numpy tools). Fixed with dynamic detection + entry-point fallback.
- **Stale `.pyc` masking edits**: edits to `utils.py` not reflected in container runs because Python was loading `utils.cpython-312.pyc`. Fixed by: (a) deleting stale pyc files from host, (b) adding `-B` to prevent re-creation.
- **`TypeError: No matching definition` in numba** (masked_rho4pairs): read-only pandas `.values` array. Fixed: `.values.copy()`.
- **`KeyError: 'target'` in `modules4top_factors`**: pandas 3.0 breaking change — `groupby.apply` excludes groupby key from result. Fixed: `sort_values + groupby.head(n)`.
- **`lima` CLI defaulted to "default" VM**: Claude Code's Bash tool used `lima bash -c ...` which targets the "default" Lima instance (doesn't exist). Fixed: `limactl shell bioinfo bash -c ...`.
- **Bash tool runs on macOS host**: `./scripts/run_scenic.sh` requires `apptainer` which is inside the Lima VM, not on macOS. Fixed by using `limactl shell bioinfo` to reach the VM.

### Decisions
- **Patch `cache/python_pkgs/pyscenic/utils.py` directly** rather than pinning older pandas: the conda/pip environment inside the container is fixed at build time; downgrading pandas 3.0 → 2.x would break other packages. Minimal surgical patches to `utils.py` (2 lines) are safer and easier to document.
- **`sort_values + groupby.head` over `include_groups=True`**: `include_groups=True` is deprecated in pandas 3.0 (warning) and will be removed in 3.1. The `sort_values + head` idiom is semantically equivalent to `nlargest` per group and forwards-compatible.
- **`-B` flag permanently added**: rather than clearing pyc manually between edits, `-B` prevents the problem category entirely for this workflow. Acceptable performance trade-off (pyc provides <100 ms startup savings for a >30 min job).

### Compartment status
- **Stroma** ⏳: not the focus of this session
- **HSC** ✅: pySCENIC completed — `regulons.csv` (11.7 MB), `auc_mtx.csv` (145.1 MB) written to `output/scenic/`
- **MNC** ⏳: not the focus of this session
- **HSC+MNC** ⏳: not the focus of this session

### Next steps
- [ ] Load `auc_mtx.csv` in `notebooks/scenic-hsc.qmd` and render the results notebook
- [ ] Check regulon count and top regulons — expect 100–400 regulons for 53k HSC cells
- [ ] Velocyto SLURM array job: BAM files on remote cluster — need hostname, barcode paths, GTF path, output destination
- [ ] Render `cluster-annotation.qmd` and verify all 9 reference sections complete cleanly
- [ ] Select final cluster resolution per compartment from ARI curves

---

## 2026-04-21 — CellRank + pySCENIC Python workflows + notebook HTML fixes

### What changed
- `notebooks/scenic_run.py` (new file):
  - Full 3-step pySCENIC workflow: GRN (GRNBoost2) → cisTarget motif pruning → AUCell scoring
  - Stratified subsample of 5,000 cells (by cluster) for GRN; all cells for AUCell
  - Resume logic: each step checks for existing output before re-running (skip if file exists and non-empty)
  - Outputs: `adjacencies.csv`, `regulons.pkl`, `regulons_df.csv`, `auc_mtx.csv`, `regulon_summary.csv`, `rss_per_cluster.csv`, `hsc_scenic.h5ad`
  - Figures: regulon size histogram, top-9 UMAP, regulon × cluster heatmap (z-scored), regulon × condition heatmap, RSS heatmap
  - Packages from `/project/cache/python_pkgs` (pyscenic-0.12.1, arboreto-0.1.6, ctxcore, dask, loompy)
  - Correct aertslab cisTarget DB path: `mc_v10_clust/gene_based/` (not `mc_v10_nr_clust/`)
- `notebooks/scenic-hsc.qmd` (new file):
  - Sections: Overview, regulon overview (DT table), size distribution, top-9 UMAPs, regulon × cluster heatmap (pheatmap, z-scored), regulon × condition heatmap, RSS heatmap + table, regulon–target flat table (DT), key HSC TFs table, Methods
  - Methods covers GRN/cisTarget/AUCell algorithms, database provenance, parameters table, troubleshooting table
  - `show_fig()` helper: inlines PNG via `xfun::base64_uri()` in `results="asis"` context; graceful skip if file missing
  - `scenic_read()` helper: wraps `read_csv` with informative error if output not yet generated
  - `embed-resources: true, theme: cosmo` frontmatter for standalone sharing
- `notebooks/velocity_run.py` (previously completed):
  - CellRank 2.x workflow: PseudotimeKernel + ConnectivityKernel (0.8/0.2) → CFLARE → fate probabilities → lineage drivers → gene trends
  - CytoTRACE2 score inverted: `cytotrace2_pseudotime = 1 - cytotrace2_score` (HIGH=stem → LOW=stem for PseudotimeKernel)
  - CFLARE chosen over GPCCA to avoid densifying 53k×53k transition matrix (~23 GB OOM)
  - Fate probabilities: GMRES solver + ILU preconditioner (avoids negative values)
  - Per-state lineage driver PNGs (CellRank internal figure; `ax=` silently ignored)
- `notebooks/cellrank-hsc.qmd` (previously completed):
  - Sections: Overview, CytoTRACE2 UMAP, kernel projection, CFLARE spectrum, terminal states, fate probabilities (DT + pheatmap + condition boxplots), lineage drivers (per-state PNGs via `xfun::base64_uri()`), gene trends
  - Detailed Methods: Markov chain math, CFLARE vs GPCCA comparison table, linear system $(I-Q)\Pi=R$, pseudotime direction callout, alternative scores table, troubleshooting guide T1–T6
- `notebooks/_quarto.yml`:
  - Added `embed-resources: true` to `format: html:` block (fixes book chapter CSS embedding)
- All other notebooks (13 files): injected `format: html: embed-resources: true, theme: cosmo` frontmatter block (fixes broken HTML when rendering single notebooks outside book context)
- `scripts/render_notebooks.sh`:
  - Added `HOME="/Users/$(whoami | sed 's/\.guest$//')"` (Lima VM `.guest` suffix fix)
  - Changed default data path from `~/data/2026.01.10_MagicMouse10X` → `~/data_hub`
  - Single-chapter render: render in-place → `mv` to `OUTPUT_DIR` → `rm -rf *_files/` (Quarto `--output-dir` not supported for non-project single files)
- `scripts/run_jupyter.sh` (new file):
  - JupyterLab on port 8888; binds repo → `/project`, data_hub → `/project/data_hub`
  - `PYTHONPATH=/project/cache/python_pkgs`; notebook dir `/project/notebooks`
  - Same Lima HOME fix as other scripts
- `scripts/download_scenic_db.sh` (new file):
  - Correct aertslab URLs (probed from Lima VM): `mc_v10_clust/gene_based/`
  - `exec > >(tee -a LOG)` dual screen+file logging; FAILED counter; `-s` guard against 0-byte files
- `src/scrna-seq/scrna-seq_utilities.R`:
  - `anno_make_type_pal`: guard for empty types (`hue_pal()(0)` crash)
  - `extract_azimuth_umap_df`: standardises `cell_type` column (tries aliases: celltype, CellType, annotation, label, cluster_label)
  - `plot_ref_umap`: guard for missing `cell_type` column in ref_df

### Why
The two major Python analyses needed for regulatory and trajectory characterisation:

**CellRank**: RNA velocity is unavailable (no BAM files); CytoTRACE2 is used as pseudotime proxy. The key insight is that CytoTRACE2 scores HIGH for stem cells, requiring inversion before passing to `PseudotimeKernel` (which directs transitions LOW→HIGH). CFLARE avoids the ~23 GB memory requirement of GPCCA's Schur decomposition on the full 54k×54k transition matrix.

**pySCENIC**: Identifies transcription factor regulons supported by both co-expression (GRNBoost2) and cis-regulatory motif evidence (cisTarget). Stratified subsampling (5k cells) keeps GRN tractable; AUCell is applied to all cells for full resolution. The RSS (Regulon Specificity Score) identifies cluster-marker regulons analogous to FindAllMarkers for genes.

The HTML render fixes address a Quarto book-project behaviour: `document-css: false` is set by the book framework for all chapters, so individual notebooks rendered outside the book context lack Bootstrap/theme CSS. The fix is `embed-resources: true` + explicit theme in individual notebook frontmatter.

### Bugs / issues
- **`plot_spectrum(show=False)` AttributeError**: CellRank 2.x estimator plot methods do not accept `show=` → removed from all CellRank estimator calls.
- **Negative fate probability values**: GMRES default tolerance too loose → `tol=1e-7, preconditioner='ilu'` resolved.
- **`Lineage.to_pandas()` AttributeError**: Not present in CellRank 2.x; replaced with `pd.DataFrame(fp_arr.X, columns=fp_arr.names, index=adata.obs_names)`.
- **`gene_trends()` missing time_key**: CellRank 2.x requires explicit `time_key` argument.
- **`plot_lineage_drivers()` ignores `ax=`**: Creates internal matplotlib figure; cannot draw on shared axis. Fixed by saving one PNG per terminal state in a loop.
- **Empty lineage drivers figure**: Same issue as above — `ax=` silently ignored.
- **Noisy kernel arrows**: CytoTRACE2 direction inverted (HIGH=stem, kernel expects HIGH=differentiated) → `cytotrace2_pseudotime = 1 - cytotrace2_score`.
- **`render_notebooks.sh` mount error**: Lima VM appends `.guest` to username; `~/data_hub` path failed → HOME override fix.
- **HTML 91KB, no CSS**: Quarto book sets `document-css: false` for chapters → `embed-resources: true` + `theme: cosmo` in frontmatter.
- **`--output-dir` for single notebooks**: Not supported by Quarto for non-project renders → render in-place then mv.
- **pySCENIC download 404 (500bp feather)**: URL used `mc_v10_nr_clust/` (doesn't exist) → correct path `mc_v10_clust/gene_based/`.
- **0-byte feather left from failed download**: `-f` guard passed 0-byte file as "exists" → added `-s` (non-empty size) check.
- **`anno_make_type_pal` crash**: `hue_pal()(0)` when all reference types absent → guard for `length(all_types) == 0`.
- **`cell_type` column missing in `ref_df`**: `extract_azimuth_umap_df` didn't standardise varying column names across references → alias resolution + rename to `cell_type`.

### Decisions
- **CFLARE over GPCCA**: densifying the 53k×53k transition matrix for Schur decomposition would require ~23 GB RAM; CFLARE uses sparse eigenvectors (scales to any number of cells).
- **GRN on 5,000-cell subsample**: GRNBoost2 XGBoost trees scale ~O(n²) in cells; 5k cells yields a representative sample while keeping runtime under 30 min. AUCell applied to all cells for full single-cell resolution.
- **Resume logic in both scripts**: each step checks for existing non-empty output files before re-running; allows restart after partial failures without re-running the slow GRN/cisTarget steps.
- **`show_fig()` via `xfun::base64_uri()`**: produces fully self-contained HTML; no external asset paths that break when the notebook is moved or shared.
- **All notebooks patched for standalone rendering**: book-context CSS fix applied globally rather than only to new notebooks, ensuring all notebooks are shareable as standalone HTML.

### Compartment status
- **Stroma** ✅: fully processed; annotation, DE, ORA, Milo in pipeline.
- **HSC** ✅: CellRank trajectory + pySCENIC regulon workflows added; notebooks complete.
- **MNC** ✅: fully processed; DE and ORA in pipeline.
- **HSC+MNC** ⏳: joint clustering targets present; no dedicated trajectory/regulon analysis yet.

### Next steps
- [ ] Run `scenic_run.py` inside container: `apptainer exec ... python notebooks/scenic_run.py`
- [ ] Re-run `velocity_run.py` after re-checking CytoTRACE2 output quality on HSC subset
- [ ] Run `data_hub/hematlas/process_hematlas.R` to build HemAtlas azimuth references
- [ ] Rebuild container with scCustomize/dittoSeq/DropletUtils/Nebulosa additions
- [ ] Re-run `ferchen_26_mm_bm_umap` and other `*_umap` targets after `extract_azimuth_umap_df` fix

---

## 2026-04-20 — External geneset pipeline + HemAtlas processing script + container updates

### What changed
- `_targets.R`:
  - Added `geneset_study_table` (14-row reference index tribble) + `geneset_study_table_path`
  - Added 9 CSV file targets (`nestorowa_17_bl_mm_csv`, `kucinski_24_csc_mm_csv`, `isobe_23_cg_mm_csv`, `hematlas_hspc_mm_csv`, `hematlas_lineage_mm_csv`, `ferchen_26_mm_marker_csv`, `ferchen_26_mm_tf_csv`, `gavish_23_mm_marker_csv`, `gavish_23_mm_tf_csv`) — mouse-native genesets
  - Added 5 more CSV targets for human-sourced genesets: `psaila_20_hs_csv`, `pan_25_hs_csv`, `gavish_23_hs_scenic_csv`, `gavish_23_hs_meta_program_csv`, `gavish_23_hs_cancer_mp_csv`
  - Added 9 `_mm_` named-list geneset targets using `readr::read_csv(...) |> with(split(...))` tidy-pipe style
  - Added 5 `_hs_` named-list geneset targets using `inner_join(mm_hs_ortho_long, ...)` for ortholog mapping before split
  - Added 6 module score targets (`module_scores_nestorowa_17_bl_mm_hsc`, `module_scores_isobe_23_cg_mm_hsc`, `module_scores_psaila_20_hs_hsc`, `module_scores_pan_25_hs_hsc`, `module_scores_gavish_23_hs_meta_program_hsc`, `module_scores_gavish_23_hs_cancer_mp_hsc`)
  - Added HemAtlas Azimuth path targets (`hematlas_mm_bm_path`, `hematlas_mm_agm_path`, `hematlas_mm_fl_path`) and projection `tar_map` blocks for BM/AGM/FL × stroma/hsc/mnc (9 targets)
  - Added `ortholog_map_mm_to_hm` (orthogene gprofiler, mouse→human, for Zeng cross-species projection)
  - Added Zeng 2025 BoneMarrowMap `tar_map` (3 compartments)
- `src/scrna-seq/scrna-seq_utilities.R`:
  - Rewrote `extract_module_scores`: now uses `dplyr::rename(all_of(rename_map))` (where `rename_map = setNames(paste0("Cluster", seq_along(genesets)), names(genesets))`) to rename AddModuleScore columns back to geneset names; returns `cell_id` + score columns only (no UMAP coords — joined in notebook)
- `container/renv_setup.R`:
  - Added `scCustomize` to `cran_pkgs`
  - Added `dittoSeq`, `DropletUtils`, `Nebulosa` to `bioc_pkgs` (scCustomize Bioconductor dependencies)
  - Added `scCustomize` to the verification list
- `data_hub/hematlas/process_hematlas.R` (new file):
  - Extracts RDS from HemAtlas tar.gz archives (deep path stripped via temp dir)
  - Rebuilds Seurat from scratch: `GetAssayData(SCT, counts)` → `CreateSeuratObject` → `SCTransform` → `PCA` → `UMAP(return.model=TRUE)` → `DietSeurat` → `saveRDS` in `{tissue}/seurat/`
  - Builds Azimuth reference: `AzimuthReference` + `ValidateAzimuthReference` + `SaveAzimuthReference` in `{tissue}/azimuth/`
  - Structure: `build_archive()` function for steps 1-4 (safe inside function); Azimuth block intentionally left raw in the `for` loop (must run in global env)
- `data_hub/tusi_2018_nature/convert_tusi_to_10x.py`:
  - Fixed barcode source: was using `meta["barcode"]` (duplicated column) for `barcodes.tsv.gz`; corrected to `meta.index.tolist()` (cell_id, the unique identifier)
- `notebooks/geneset-exploration.qmd` (new file):
  - Setup chunk: `umap_df` (UMAP coords + sample.name from `anno_meta_hsc`) shared across all plotly sections; `geneset_dt()` helper (DT table from named gene list); `score_scatter()` helper (plotly UMAP coloured by module score)
  - Sections: Resource index → Nestorowa 2017 (plotly module scores) → Kucinski 2024 (FeaturePlot split.by) → Isobe 2023 (plotly) → HemAtlas HSPC subclusters (FeaturePlot) → HemAtlas mouse lineage (FeaturePlot) → Ferchen 2026 markers (FeaturePlot) → Ferchen 2026 TF (FeaturePlot) → Gavish 2023 mm markers (FeaturePlot, no split.by) → Gavish 2023 mm TF (FeaturePlot, no split.by) → Gavish 2023 SCENIC hs→mm (FeaturePlot, filter ≥2 orthologs, `rm(srt); gc()` here) → Psaila 2020 (plotly) → Pan 2025 (plotly) → Gavish meta programs (plotly) → Gavish cancer MP & NMF (plotly)
  - All FeaturePlot sections gate on `n_avail <- lapply(lst, function(x) sum(x %in% rownames(srt)))` before plotting

### Why
The geneset exploration notebook provides a systematic HSC-focused scan against 14 published reference datasets spanning lineage markers, trajectory programmes, TF activities, SCENIC regulons, and cancer meta-programs. Module scores (plotly, interactive) are used for whole-programme summaries; individual FeaturePlots are used for marker gene inspection. Human datasets are ortholog-mapped via the existing `mm_hs_ortho_long` target (Ensembl + MGI join) before scoring.

HemAtlas archives ship with only an incomplete SCT assay (no RNA assay, no PCA, no UMAP model), so the processing script must rebuild the entire reduction stack from `SCT@counts`. The Azimuth step must run in the global environment due to environment introspection inside `AzimuthReference`.

The `scCustomize` addition improves single-cell visualization options (dot plots, proportion plots) without requiring manual installation.

### Bugs / issues
- **HemAtlas: "Cannot delete default assay"**: initial approach renamed SCT → refAssay then tried to delete the old slot; Seurat blocked deletion of the default assay. Fixed by rebuilding from raw counts entirely (no rename needed).
- **HemAtlas: No PCA / SCT has no variable features**: HemAtlas objects had only the SCT assay with no PCA reduction and no variable features stored. Resolved by extracting `SCT@counts` → `CreateSeuratObject` → full SCTransform → PCA → UMAP pipeline.
- **HemAtlas: AzimuthReference inside function**: `ValidateAzimuthReference` failed when called inside `build_archive()` — Azimuth introspects the calling frame for an `AzimuthData` object. Fixed by splitting into `build_archive()` (steps 1–4) and running the Azimuth block raw in the `for` loop, mirroring the pattern in `10x_to_azimuth.R`.
- **`extract_module_scores` style**: initial rewrite used `match()`/`paste0("score_", ...)` column renaming. User required `dplyr::rename(all_of(rename_map))` style with `setNames(paste0("Cluster", seq_along(genesets)), names(genesets))`. Corrected.
- **`n_avail` filter dropped**: availability filter (`lapply(lst, function(x) sum(x %in% rownames(srt)))`) was omitted from some notebook sections. Restored consistently throughout all FeaturePlot sections.
- **Target naming missing species**: initial geneset targets lacked `_mm_` / `_hs_` species tags. Renamed throughout (e.g. `nestorowa_17_bl_genesets` → `nestorowa_17_bl_mm_genesets`).
- **`convert_tusi_to_10x.py` barcode bug**: `meta["barcode"]` column is duplicated (not unique); unique cell identifier is `meta.index` (`cell_id`). Fixed to `meta.index.tolist()`.

### Decisions
- **Tidy-pipe geneset targets**: `readr::read_csv(...) |> with(split(...))` — one expression per target, no intermediate variables. Human targets insert `dplyr::inner_join(mm_hs_ortho_long, ...)` in the same pipe before `with(split(...))`.
- **`umap_df` centralised in setup chunk**: all plotly score sections join to `umap_df` at render time rather than embedding UMAP coords in each target, keeping module score targets lightweight (cell_id + scores only).
- **`srt` lifetime in notebook**: loaded once in Kucinski chunk (first FeaturePlot section), freed with `rm(srt); gc()` after Gavish SCENIC (last FeaturePlot section); plotly sections load only the lightweight data frames.
- **Gavish SCENIC filter ≥ 2 orthologs**: plotting a single-gene FeaturePlot is uninformative and prone to noise; filter `lst_avail > 1` applied only for the SCENIC section.

### Compartment status
- **Stroma** ✅: fully processed; DE, ORA, Milo, annotation projection targets all in pipeline.
- **HSC** ✅: fully processed; geneset exploration notebook now covers 14 reference datasets with module scores and FeaturePlots.
- **MNC** ✅: fully processed; DE and ORA targets in pipeline.
- **HSC+MNC** ⏳: not yet implemented (explored and reverted this session).

### Next steps
- [ ] Run `data_hub/hematlas/process_hematlas.R` inside the container to produce BM/AGM/FL seurat/ + azimuth/ outputs
- [ ] Run `data_hub/tusi_2018_nature/convert_tusi_to_10x.py` (after barcode fix) then `10x_to_azimuth.R` for Tusi tissues
- [ ] Rebuild container with new scCustomize/dittoSeq/DropletUtils/Nebulosa additions
- [ ] Render `geneset-exploration.qmd` and review module score patterns across HSC clusters

---

## 2026-04-04 19:47 — DE notebook: ORA replacing GSEA, bubble plots, geneset membership tables, biased gene expression UMAPs

### What changed
- `notebooks/differential-expression.qmd`: extensive overhaul of the all_vs_all section for all three compartments
  - Replaced `fgsea::fgsea()` (ranked GSEA) with `fgsea::fora()` (ORA / Fisher's exact test) for HSC, MNC, and Stroma
  - Removed NES heatmaps; added ORA bubble plots (x = geneset, y = pathway; fill = log2 foldEnrichment; size = overlap/size ratio; black stroke = significant)
  - Added geneset membership tables split into summary (n_genes per geneset) and full membership (one row per gene) for all compartments
  - Added biased gene expression UMAP section: off-diagonal genes (mCh⁺-biased / mCh⁻-biased) plotted individually split by `sample.name`, grouped by geneset with `patchwork::plot_annotation` title
  - Added Stroma ORA + module score + biased expression sections (Stroma previously had no ORA or geneset visualisation)
  - Fixed DT widget rendering in `results="asis"` loops: replaced `print(DT::datatable())` with `print(htmltools::tagList(DT::datatable()))` 
  - Volcano y-axis: `expansion(mult = c(0, 0.15))` — no padding below zero
  - Geneset membership scatter: raised to `fig.height=9, fig.width=9` for symmetric 2×2 facet grid
  - Bubble plot: x-axis moved to top (`scale_x_discrete(position="top")`); fig sizes capped with `min()/max()` formula
  - Table columns: sorted by `foldEnrichment` descending; `overlap` shown as `overlap/size` ratio string
- `src/scrna-seq/scrna-seq_utilities.R`: `viridis_option` parameter added to `lean_dimplot`, `plot_umap_feature`, `plot_umap_feature_ext`

### Why
GSEA on a small per-geneset ranked list is methodologically flawed — ranking the few genes within a geneset defined by logFC position relative to the diagonal does not meaningfully test pathway enrichment. ORA via `fgsea::fora()` is the correct approach: it asks whether the geneset is over-represented in a MSigDB pathway relative to a consistent background (all significant scatter genes). The biased gene expression UMAPs directly show where mCh⁺- and mCh⁻-specific genes are expressed across cell states, complementing the module score visualisation with single-gene resolution.

### Bugs / issues
- **`htmltools::as.character()` not exported**: attempted `cat(htmltools::as.character(DT::datatable(...)))` in `results="asis"` — function not exported from htmltools namespace. Fixed by `print(htmltools::tagList(DT::datatable(...)))`.
- **Cairo surface size exceeded**: `plot_umap_feature_ext` stacks one row per feature (`ncol=1` patchwork); with 40 unique biased genes × 4 inch/row = 160 inch height, Cairo exceeded its ~32767 px limit. Fixed by iterating per geneset (≤20 genes each) rather than plotting all features in one call; `fig.height=30` is safe for 10-gene batches.
- **`fora()` `foldEnrichment` column**: initially computed fold enrichment manually as `(overlap/size) / (geneset_size/universe_size)`; discovered `fora()` already returns `foldEnrichment` column directly. Replaced manual calculation.

### Decisions
- **ORA over GSEA**: `fgsea::fora()` chosen over `clusterProfiler::enricher()` because fgsea is already loaded and `fora()` provides a consistent interface. Universe fixed to `unique(scatter_df$feature)` across all geneset tests within a compartment to ensure comparability.
- **Per-geneset iteration for expression UMAPs**: splitting by geneset (not arbitrary batches) naturally groups related genes, allows `patchwork::plot_annotation(title=gs_name)` to label each block, and keeps each figure within Cairo's size limit.
- **Stroma geneset definition**: single comparison (Het vs WT) → no scatter possible; genesets defined as `up_Het` / `down_Het` from the volcano (sig & |logFC| ≥ 1). Top 20 genes ranked by |logFC| for expression plots.
- **Bubble plot aesthetics**: significant hits filled + black stroke (0.3); non-significant grey (fill="grey88") with white stroke — avoids the heavy border from previous version while still distinguishing significance. Log2 fold enrichment chosen over raw to compress dynamic range.

### Compartment status
- **Stroma** ✅: full pipeline + ORA + geneset module scores + biased expression UMAPs now in DE notebook
- **HSC** ✅: full ORA pipeline; biased expression plots per geneset with geneset label
- **MNC** ✅: full ORA pipeline; biased expression plots per geneset with geneset label
- **HSC+MNC** ✅: full pipeline complete; no DE notebook section

### Next steps
- [ ] Run `tar_make()` to build all `mpn_subclust_*` targets
- [ ] Render DE notebook and review ORA bubble plots for biological signal
- [ ] Annotate clusters using canonical marker bubble plots
- [ ] Share updated notebook book with Camélia

---

## 2026-03-31 21:30 — Notebook render output restructured: versioned snapshots under analysis_book/

### What changed
- `notebooks/_quarto.yml`: removed `output-dir: _book` — output dir now passed dynamically at render time
- `scripts/render_notebooks.sh`: replaced fixed `notebooks/_book/` output with timestamped `--output-dir /project/analysis_book/YYYY-MM-DD_HH-MM`; removed `run_label` dependency (renders are independent of pipeline runs)
- `.gitignore`: replaced `notebooks/_book/` entry with `analysis_book/` to cover versioned snapshots
- `notebooks/.gitignore`: removed (fully redundant with root `.gitignore`; root already covers `notebooks/_freeze/`, `notebooks/*_cache/`, `notebooks/*_files/`)
- `notebooks/_targets.yaml`: removed (redundant — all `.qmd` notebooks explicitly call `tar_config_set(store = here::here("_targets"))` in setup chunks)
- `notebooks/_book/`, `notebooks/cell-cycle_files/`, `notebooks/qc-exploration_cache/`, `notebooks/qc-exploration_files/`: deleted render artifacts that had accumulated in the source directory

### Why
The `notebooks/` directory was accumulating large rendered outputs (`_book/` ~171 MB, `_freeze/` ~217 MB, `*_cache/`, `*_files/`) alongside the source `.qmd` files, making it hard to see what was source vs generated. The user also wanted to preserve rendered book snapshots at different project stages without overwriting them — each render should produce an independent, dated snapshot for sharing with collaborators.

### Bugs / issues
- `render_notebooks.sh` initially parsed `run_label` from config to name the output dir. Rejected — notebook renders are independent of pipeline runs and shouldn't require updating `project_config.yaml` before rendering. Switched to automatic `date '+%Y-%m-%d_%H-%M'` timestamp instead.
- Passing `--analysis_book` as a CLI argument to the script caused Quarto to fail with "No valid input files" — the script was treating it as a chapter filename. This was user error (testing the script); no fix needed.

### Decisions
- **Timestamp over run_label for output naming:** Pipeline runs are versioned by `run_label` (manually updated in config). Notebook renders are exploratory and should be self-dating — a timestamp is simpler and doesn't require any config update before rendering.
- **Single root `.gitignore` over nested ones:** The `notebooks/.gitignore` was fully covered by explicit `notebooks/*` patterns in the root `.gitignore`. Removed to reduce confusion about where ignore rules live.
- **Keep `_freeze/` and `.quarto/` in notebooks/:** These are Quarto's execution cache and metadata. They belong next to `_quarto.yml` and cannot be relocated without restructuring the project root. They are gitignored at the root level.
- **Remove `_targets.yaml` from notebooks/:** All `.qmd` notebooks already call `tar_config_set()` explicitly. The file was a redundant safety net.

### Compartment status
- **Stroma** ⏳: pipeline targets in place; no new changes this session
- **HSC** ⏳: pipeline targets in place; no new changes this session
- **MNC** ⏳: pipeline targets in place; no new changes this session
- **HSC+MNC** ⏳: pipeline targets in place; no new changes this session

### Next steps
- [ ] Run `./scripts/render_notebooks.sh` to verify book renders to `analysis_book/YYYY-MM-DD_HH-MM/`
- [ ] Confirm `notebooks/` stays clean after render (no `_book/` or `*_files/` appearing)
- [ ] Share `analysis_book/<timestamp>/index.html` with Camélia

---

## 2026-04-01 — Notebook QC/DE visualisations + PCA distribution targets

### What changed
- `notebooks/first-analysis-milestone-report.qmd`: added UMAP visualisations of `sample.name`, `scDblFinder.class`, and `scDblFinder.score` to both the pre-QC (section 3) and post-QC dimensionality reduction (section 4.1) sections for all three subsets; added `sample.name × dm2_class` UMAP grid (plain + highlight) for both pre- and post-QC; restructured `# Markers` into four subsections: `each_vs_rest`, `all_vs_all`, `mCh- vs mCh+`, `mCh- vs WT`; added explicit `tar_load` + `rm()`/`gc()` discipline throughout to prevent OOM kills during rendering.
- `src/scrna-seq/scrna-seq_utilities.R`: added `plot_umap_grid_sample_dm2()` (sample × dm2_class UMAP grid, manual highlight via stripped background layer); added `find_delegate_markers_all_vs_all()` (DElegate `findDE`, `compare = "all_vs_all"`); added `find_delegate_markers_mch()` (Het mCh- vs Het mCh+, HTO tag as replicate); added `find_delegate_markers_mch_vs_wt()` (Het mCh- vs WT, HTO tag as replicate); added `extract_pca_df()`, `plot_pc_parcoor()`, `plot_pc_umap_grid()` (PCA distribution analysis, moved from notebook).
- `_targets.R`: added `delegate_allvsall_markers`, `delegate_mch_markers`, `delegate_mch_vs_wt_markers` targets inside `tar_map`; added `pca_df` target (lightweight data.frame, no Seurat object in cache).
- `notebooks/pca-cell-distribution.qmd`: rewritten to load `pca_df_*` targets; explicit named chunks per subset; `N_PCS` from `proc_params`; plotting functions moved to utilities.
- `scripts/render_notebooks.sh`: timestamp changed from `YYYY-MM-DD_HH-MM` to `YYYY-MM-DD`; added optional `--label`/`-l` flag that concatenates to the timestamp (e.g. `2026-04-01_first-milestone`).

### Why
Camélia review preparation: the report needed richer QC visualisations (doublet scores on UMAP, sample distribution on embedding) and a more complete DE picture. The `each_vs_rest` DElegate comparison conflates the three genotype groups — pairwise and mCh-specific comparisons are needed to isolate the mCh signal from sex-driven differences. The PCA distribution notebook was refactored to integrate with the targets pipeline so that the heavy Seurat objects are not re-loaded during rendering.

### Bugs / issues
- **`gghighlight` incompatibility:** `gghighlight 0.5.0.9000` (dev) and even `0.4.1` (CRAN) both fail with `is.na()` on a language object inside their aes-param validation loop when used with the current ggplot2 version. Replaced with a manual two-layer approach: background `geom_point` uses `d_bg` (coordinate columns only, no faceting variables) so `facet_grid` does not subset it per panel; foreground layer inherits the full data and is subsetted normally.
- **`delegate_mch_markers_stroma` error:** Stroma samples are named `Stroma_Het` / `Stroma_WT` — no `mCh-`/`mCh+` pattern. Function errored with "No cells found". Fixed by adding an early-return guard returning `dplyr::tibble()` when `keep_cells` is empty.
- **`pca_df` target dependency not detected:** `colnames(mpn_processed_data@meta.data)` inside the target body — `@` slot access prevents static dependency detection (known targets gotcha). Fixed by assigning `srt <- mpn_processed_data` first.
- **OOM kill at chunk 102/131:** All three `mpn_processed_*` objects were simultaneously live inside a tabset. Fixed by adding `rm()`/`gc()` at the end of each grid chunk and adding explicit `tar_load` to all downstream chunks that assumed objects persisted.
- **`bubble-hsc` "object not found":** Same root cause — `mpn_processed_hsc` freed early but `bubble-*` chunks had no `tar_load`. Added `tar_load` + end-of-chunk `rm()`/`gc()` to all three canonical marker chunks.
- **`tar_map` hash invalidation:** Adding `pca_df` inside the `tar_map` changed the map's hash, triggering a full rebuild. Targets correctly skips unchanged upstream objects; only `pca_df_*` actually re-runs.

### Decisions
- **Manual highlight over `gghighlight`:** Rather than pinning to a specific `gghighlight` version or patching the package, implemented the highlight pattern directly: strip faceting columns from the background layer so ggplot2 doesn't subset it per facet. More transparent and dependency-free.
- **`dm2_class` as replicate column for mCh comparisons:** Within a Het tube, each HTO tag represents an independently demultiplexed biological replicate. Using it as `replicate_column` in `findDE` gives DElegate the replicate structure it needs for pseudobulk DE — much more statistically sound than treating all Het mCh- cells as one sample.
- **Empty-tibble guard instead of error:** Stroma has no mCh labelling by design. Silently returning an empty tibble lets the target succeed and the notebook render a blank table, making it clear the comparison is not applicable rather than crashing the pipeline.
- **`pca_df` as plain data.frame target:** Storing only the embedding + metadata (not the Seurat object) keeps the cache small and makes the notebook load instantaneous. `extract_pca_df` is the extraction boundary.
- **Explicit `tar_load` + `rm()`/`gc()` per chunk:** Tabset chunks all execute before Quarto moves to the next section, so all objects in a tabset coexist. Enforcing load/free within each chunk is the only reliable memory management strategy on a 16 GB machine.

### Compartment status
- **Stroma** ⏳: all targets built; `delegate_mch_markers_stroma` returns empty (no mCh labelling); `delegate_mch_vs_wt_markers_stroma` completed.
- **HSC** ⏳: all new DE and PCA targets built.
- **MNC** ⏳: all new DE and PCA targets built.
- **HSC+MNC** ⏳: no changes this session.

### Next steps
- [ ] Run `delegate_allvsall_markers_*` targets (not yet built)
- [ ] Render full notebook book and share with Camélia
- [ ] Review mCh- vs mCh+ DE results — check if sex-gene signal persists or drops out
- [ ] Annotate clusters using canonical markers bubble plots

---

## 2026-04-02 — Notebook visualisation overhaul: scatter/volcano DE plots, QC notebook rewrite, dim-reduction notebooks fixed

### What changed
- `notebooks/first-analysis-milestone-report.qmd`:
  - **Section 5.2 all_vs_all / HSC**: added scatter plot comparing log FC of Het mCh− vs WT (x) × Het mCh+ vs WT (y) for CD41+ HSC clusters; added per-cluster volcano plots (Het mCh− vs WT and Het mCh+ vs WT faceted); both use shape 21 / white border / alpha 0.5 / fill by cluster or comparison; point size = max(rate1, rate2); labels via ggrepel.
  - **Scatter plot**: significant genes only (`sig_neg | sig_pos`); labels top 20 per quadrant ranked by Euclidean distance from origin; all-clusters version precedes per-cluster tabset.
  - **Volcano plots**: per-cluster tabset; fill/colour by comparison (orange = mCh−, blue = mCh+); labels for `sig & (|log FC| ≥ 2 | −log10(padj) ≥ 20)`.
  - **Bug fix**: `inner_join` in scatter-prep was filtering `sig` on each side before joining, silently dropping genes significant in only one comparison. Fixed by joining full tables first then filtering `sig_neg | sig_pos` post-join.
- `notebooks/qc-exploration.qmd`: full rewrite — section 1 replaced with milestone-report jitter+violin+percentile-errorbar distributions (5 metrics, tabset); section 2 replaced with feature scatter plots copied from milestone report including MNC double-curve analysis; new section 2.5 doublet/multiplet overlap heatmap; section 3 `fig-height` increased to 12 for 3-stacked sweep plots; `hsc_mnc` removed from section 7 `PROCESSED_NAMES`; `str_detect` → `grepl` fix for missing `stringr`.
- `notebooks/dim-reduction-and-clustering-SCT.qmd`:
  - Added Section 1 loading cached `mpn_processed_data_*` targets for pipeline-identical UMAP.
  - **Core fix**: grid loop now calls `process_all()` directly (same function as `_targets.R`) instead of manually reimplementing `preprocess_seurat → find_neighbors_clusters → RunUMAP`. Confirmed visually identical to cached objects.
  - `umap-sample-tabs` fig size increased to 15×10 for 3×2 grid.
- `notebooks/dim-reduction-and-clustering-RNA.qmd`:
  - Grid loop restructured to match SCT pattern: `preprocess_rna()` helper defined, called fresh per grid combination from `raw_srt`; `raw_srt` kept in scope across the loop.
  - **Bug fix**: `preprocess_rna` was wrapping the pipe chain in `capture.output(...)` without `{}` braces + explicit assignment — function returned `NULL` (the captured output lines) instead of the Seurat object. Fixed with `out <- NULL; capture.output({ out <- ... }); out`.
  - `umap-sample-tabs` fig size increased to 15×10.
  - Stale explanatory text and inline comments removed from both dim-reduction notebooks.

### Why
The all_vs_all DE section needed richer visualisation: a scatter comparing the two WT comparisons directly reveals which genes are concordantly dysregulated in both Het mCh− and Het mCh+ vs WT, and which are specific to one genotype. Volcano plots per cluster give a per-cluster view of significance and effect size. The QC notebook was a simplified placeholder that didn't reflect the full QC exploration already in the milestone report. The dim-reduction notebooks had a critical workflow mismatch — they reimplemented `process_all()` manually and used PCA-based `RunUMAP` instead of NN-graph-based, producing different embeddings from the pipeline.

### Bugs / issues
- **Scatter `inner_join` pre-filtering bug**: filtering `sig` before joining retained only genes significant in *both* comparisons (implicit `&`). The intent was `|` (significant in at least one). Fixed by joining full tables, then filtering `sig_neg | sig_pos`.
- **`preprocess_rna` returns NULL**: `capture.output(pipe_chain)` returns the captured text lines, not the pipe result. Fixed by using `{}` block with explicit `out <-` assignment inside.
- **SCT notebook UMAP mismatch**: the notebook used `RunUMAP(..., dims = seq(n_dims), reduction = pca_name)` (reduction-based) while `process_all()` uses `RunUMAP(..., nn.name = Seurat::Neighbors(srt))` (NN-graph-based). Fixed by replacing the manual reimplementation with a direct call to `process_all()`.
- **`str_detect` not found in qc-exploration.qmd**: `tidyverse`/`stringr` not loaded. Replaced with base `grepl`.

### Decisions
- **`process_all()` as the single source of truth**: the grid loop in the SCT notebook now calls `process_all()` directly rather than reimplementing it. Any future changes to the pipeline function are automatically reflected in the notebook.
- **Quadrant-based labelling for scatter**: labelling top N genes per quadrant (by distance from origin) ensures labels are distributed across all four corners rather than concentrated around the largest absolute values in one direction.
- **Shape 21 + white border**: open-feeling points allow overlapping genes to remain distinguishable; alpha 0.5 further reduces overplotting while keeping colour legible.
- **Join-then-filter pattern**: always join full tables before applying significance filters to avoid silent data loss when conditions are applied independently per side.

### Compartment status
- **Stroma** ⏳: no pipeline changes; visualisation improvements in notebooks.
- **HSC** ✅: all_vs_all scatter + volcano plots complete; CD41+ HSC comparisons visualised.
- **MNC** ⏳: no pipeline changes; notebook QC section updated.
- **HSC+MNC** ⏳: no changes this session.

### Next steps
- [ ] Apply same scatter + volcano treatment to MNC and Stroma all_vs_all sections
- [ ] Render updated notebook book and share with Camélia
- [ ] Annotate clusters using canonical marker bubble plots
- [ ] Review concordant genes in HSC scatter (significant in both mCh− and mCh+ vs WT)

---

## 2026-04-03 — Compositional analysis notebook + sub-clustering pipeline

### What changed
- `notebooks/compositional-analysis.qmd`: new notebook added — cluster proportion bubble plots (sample + sample×HTO grid), row-scaled heatmap (pheatmap, % labels), and Milo differential abundance (miloR) for Stroma, HSC, MNC. Structured as two top-level sections: *Cluster proportions* and *Milo differential abundance*, each tabbed by compartment.
- `notebooks/subcluster-exploration.qmd`: new notebook — per-sub-cluster UMAP coloured by `sample.name` (re-processed embedding) + silhouette on parent compartment UMAP, pie chart of cell counts/%, per-sample highlight rows. Compartment sections with cluster tabs. Consistent `sample_pal` across all panels.
- `notebooks/_quarto.yml`: added `compositional-analysis.qmd` and `subcluster-exploration.qmd` to Supplementary section.
- `_targets.R`: added sub-clustering `tar_map` block — reads `subclusters:` from `project_config.yaml` at parse time via `subcluster_values`; produces `mpn_subclust_<name>` (subsetted + reprocessed Seurat) per entry. No markers or DE targets — kept to processing only.
- `src/scrna-seq/scrna-seq_utilities.R`: added `subset_clusters(srt, clusters, resolution)` — greps cluster column, subsets cells, adds `parent_cluster` metadata column.
- `project_config.yaml`: added `subclusters:` section with 17 entries — `hsc_c0/2/3/4/6`, `mnc_c0/1/2/3/4/6/7/8`, `stroma_c0/1/2/3`; moved `marker_min_pct`, `marker_logfc_threshold`, `marker_max_cells_per_ident` back to `processing:` (were accidentally nested under `subclusters:`).
- `container/renv_setup.R`: added `miloR` to `bioc_pkgs` and verification list.
- `container/bioinfo.def`: added `miloR` smoke test.

### Why
Camélia requested a compositional analysis to complement the cluster marker and DE results — understanding which clusters expand or contract between WT, Het mCh−, and Het mCh+ conditions is central to the MPN phenotype characterisation. Milo provides neighbourhood-level resolution beyond fixed clusters. Sub-clustering was requested to re-process individual clusters of interest (e.g. MK-like clusters in HSC, specific MNC populations) through the full SCTransform → PCA → UMAP pipeline to better visualise within-cluster heterogeneity by condition.

### Bugs / issues
- **Milo `testNhoods` design matrix mismatch**: rownames of design df were not set — fixed by `rownames(design_df) <- design_df$milo_sample`.
- **Milo NA dispersions**: Stroma has only 1 tube per condition → edgeR cannot estimate dispersion. Fixed by using `dm2_class` HTO tags as within-condition replicates (`milo_sample = paste(condition, dm2_class, sep = "__")`), consistent with the pseudobulk DE approach already used in this project.
- **`plotNhoodGraphDA` missing neighbourhood graph**: `buildNhoodGraph()` call was absent — added after `makeNhoods()`.
- **`plotNhoodGraphDA` UMAP name**: function expects `"UMAP"` in `reducedDim`; reduction is named `"umap.nn"`. Fixed by copying it under the canonical name inside `run_milo()`.
- **`plotNhoodGraphDA` `res` argument**: was passing `res = 3` (numeric) which was interpreted as `res_column` (string) argument — removed.
- **`JoinLayers` on SCT assay**: `subset_clusters()` called `JoinLayers()` on already-processed object; SCT assay doesn't support it. Removed the call entirely — subsetting a processed Seurat object does not split layers.
- **`subcluster_values` placement**: initially defined as a top-level R variable reading from yaml; then attempted inline `{...}` block in `tar_map(values = ...)` which is invalid. Resolved by keeping it as a top-level variable, later inlined properly into `tar_map(values = {...})` after user clarification — then reverted to yaml-reading pattern as user confirmed config is the right source of truth.

### Decisions
- **Milo replicates from HTO tags**: rather than treating each tube as one replicate (giving n=1 per condition for Stroma), HTO-demultiplexed `dm2_class` labels are used as pseudo-replicates within each condition. This is consistent with how `find_delegate_markers_mch()` works and gives edgeR enough degrees of freedom.
- **Sub-clustering targets produce only `mpn_subclust_<name>`**: no markers or DE targets added — user confirmed processing only is needed at this stage; markers/DE can be added later per specific cluster.
- **Sub-cluster naming**: `<compartment>_c<cluster>` pattern (e.g. `hsc_c0`, `stroma_c2`) chosen over plain `<compartment>` to support multiple sub-cluster runs per compartment. Config `name` field drives target names; notebook reads names dynamically from yaml.
- **Notebook loads config dynamically**: `subcluster-exploration.qmd` reads `subclusters:` from `project_config.yaml` at render time so adding new entries doesn't require editing the notebook.
- **Pie chart over bar**: user requested pie chart for the sample proportion overview; stacked bar replaced with `coord_polar(theta = "y")`.
- **UMAP legend top**: moved to `position = "top"` + `direction = "horizontal"` to free horizontal space for the UMAP panel.

### Compartment status
- **Stroma** ✅: full pipeline complete; `stroma_c0/1/2/3` sub-clusters defined and queued.
- **HSC** ✅: full pipeline complete; `hsc_c0/2/3/4/6` sub-clusters defined and queued.
- **MNC** ✅: full pipeline complete; `mnc_c0/1/2/3/4/6/7/8` sub-clusters defined and queued.
- **HSC+MNC** ✅: full pipeline complete; no sub-clusters defined yet.

### Next steps
- [ ] Run `tar_make()` to build all `mpn_subclust_*` targets
- [ ] Render `subcluster-exploration.qmd` and review UMAP patterns per cluster
- [ ] Annotate clusters using canonical marker bubble plots
- [ ] Apply scatter + volcano treatment to MNC and Stroma all_vs_all sections in milestone report

---

## 2026-04-03 17:00 — Differential expression notebook + Milo pairwise contrasts + scatter plot polish

### What changed
- `notebooks/differential-expression.qmd` (NEW): full markers section extracted from `first-analysis-milestone-report.qmd`; contains `each_vs_rest`, `all_vs_all`, `mCh- vs mCh+`, `mCh- vs WT` subsections for HSC; interactive plotly volcano (per cluster) and interactive plotly scatter (all clusters collectively); diagonal reference lines on all scatter plots.
- `notebooks/first-analysis-milestone-report.qmd`: entire `# Markers` section removed (~900 lines); stale freeze cache cleared (`rm -rf notebooks/_freeze/first-analysis-milestone-report`) after chunk numbering shifted.
- `notebooks/compositional-analysis.qmd`: Milo section overhauled — dropped `model.contrasts` in favour of two independent two-group runs (`.run_milo_pair()`); added `run_milo_contrasts()` wrapper; `plot_milo_contrasts()` renders 2-row layout (beeswarm | beeswarm / volcano | volcano); heatmap style improved (`#2166ac`/`#b2182b` palette, `border_color=NA`, `cellwidth=52`, `cellheight=22`, condition annotation bar, `"%.1f%%"` labels, `treeheight_row/col=0`).
- `notebooks/_quarto.yml`: added `differential-expression.qmd` as first chapter in Supplementary section.

### Why
The milestone report was growing unwieldy with ~900 lines of marker/DE code. Moving interactive plots and the full markers section to a dedicated notebook keeps the report focused on high-level interpretation while allowing deeper DE exploration in a separate document. The Milo `model.contrasts` approach was producing nonsensical results (identical ±40 logFC for both contrasts) because `testNhoods` exposes raw GLM coefficients on log-CPM scale, not fold-changes — two independent two-group runs are the correct approach for pairwise comparisons. The scatter plot diagonal lines (`y=x`, `y=x±1`) were requested to visually separate shared Het effects (along the diagonal) from mCh-specific effects (above/below the band).

### Bugs / issues
- **`model.contrasts` producing garbage logFC**: `testNhoods` with `model.contrasts` returns per-coefficient log-CPM columns (`logFC.1`, `logFC.0`, `logFC..1`), not pairwise fold-changes. Both contrasts produced identical ±40 values. Root cause: miloR's `model.contrasts` is not analogous to limma contrasts. Fix: dropped `model.contrasts` entirely; run two independent two-group Milo analyses via `.run_milo_pair()`, which subsets cells to just the two conditions of interest and fits a simple `~ condition` model.
- **Beeswarm empty for mCh+**: first attempt passed full three-condition object to `testNhoods` with a contrast vector; the `condition` factor still had three levels so the second run returned empty results. Fix: `.run_milo_pair()` subsets and re-factors to exactly two levels before building the Milo object.
- **Plotly widgets blank in `results="asis"` chunks**: knitr streams raw HTML text in `results="asis"` mode but does not register widget JS dependencies. Fix: use `htmltools::tagList(lapply(...))` as the last expression of a normal chunk (no `results="asis"`).
- **Stale freeze cache after markers removal**: after deleting ~900 lines from milestone report, Quarto's freeze metadata referenced wrong chunk line numbers, causing `load-raw-stroma` to error. Fix: `rm -rf notebooks/_freeze/first-analysis-milestone-report`.
- **`top_n_diagonal` defined inline**: function was defined inside the plot chunk rather than the shared helpers chunk, making it invisible to other chunks. Fix: moved to `allvsall-hsc-helpers` alongside `top_n_quadrant`, `repel_by_quadrant`, `repel_by_diagonal_half`.

### Decisions
- **Independent two-group Milo runs over `model.contrasts`**: miloR's GLM interface does not expose pairwise log-fold-changes when `model.contrasts` is used — it returns model coefficients. Running two separate analyses (mCh− vs WT, mCh+ vs WT) each with a binary `condition` factor is both simpler and produces correct, interpretable logFC values from `testNhoods`.
- **Separate DE notebook**: the milestone report should remain a readable summary for collaborators; all exploratory marker and DE content (static + interactive) moved to `differential-expression.qmd` in Supplementary. This also removes the `plotly`/`htmltools` dependencies from the milestone report.
- **`htmltools::tagList` for plotly loops**: the only reliable way to render multiple plotly widgets in a Quarto chunk is to return a `tagList` — not to print in a loop or use `results="asis"`. Adopted as the pattern for all future interactive widget sections.
- **Diagonal reference lines on scatter**: `y=x` (slope 1, intercept 0) marks genes with identical effect in both conditions (shared Het effect); `y=x+1` / `y=x-1` mark the ±1 logFC band — genes outside this band are specifically enriched in one condition. Colour `#4e9a8a` (teal) chosen to contrast with the orange/blue point palette without clashing.
- **Collective interactive scatter**: user requested one interactive scatter covering all clusters simultaneously (not per-cluster tabs), making it easy to compare cross-cluster gene behaviour in a single view.

### Compartment status
- **Stroma** ✅: full pipeline complete; compositional analysis (Milo) now includes two pairwise contrasts.
- **HSC** ✅: full pipeline + DE analysis in new dedicated notebook; scatter + volcano with interactive versions.
- **MNC** ✅: full pipeline complete; DE section (scatter/volcano) not yet added to DE notebook.
- **HSC+MNC** ✅: full pipeline complete; no DE notebook section yet.

### Next steps
- [ ] Add MNC and Stroma all_vs_all scatter + volcano to `differential-expression.qmd`
- [ ] Run `tar_make()` to build all `mpn_subclust_*` targets
- [ ] Render `subcluster-exploration.qmd` and review UMAP patterns per cluster
- [ ] Annotate clusters using canonical marker bubble plots
- [ ] Share updated notebook book with Camélia

---

## 2026-04-12 — Cluster annotation pipeline + Azimuth label transfer + DE/ORA targets moved to pipeline

### What changed
- `_targets.R`: major restructuring and large new annotation block
  - Removed `output_dir` stamped-run-key target; replaced with static `path_params` (reads `paths:` block from yaml — `output_table_dir`, `output_figure_dir`, `azimuth_dir`)
  - Removed all per-subset file-output targets from `tar_map` (umap_*, elbow_*, pca_*, cluster_markers_path_*, bubble_plot_*, marker_feature_plots_*) and manual hsc_mnc equivalents — output targets deleted to simplify pipeline; notebooks now load from targets store directly
  - Removed `output_dir` and `qc_violin_path` targets; QC CSV paths now reference `path_params$output_table_dir`
  - Removed `report` target (Quarto book render) and `analysis_src` tracker
  - Added `tar_config_set(store = here::here("_targets"))` at top of `_targets.R`
  - Added `miloR`, `SingleCellExperiment`, `fgsea`, `msigdbr` to `tar_option_set(packages = ...)`
  - Added `annotation_recluster` + `annotation_recluster_path` targets (`tar_map` over stroma/hsc/mnc) — `recluster_resolutions_df()` runs `FindClusters` at 8 resolutions (0.1–0.9) and returns a plain data.frame
  - Added reference loading targets: `tabula_muris_marrow/spleen/lung/kidney` (cellxgene census via `fetch_ref_census()`), `marrow_atlas` (pre-built syn53222728.Rds via `load_ref_marrow_atlas()`), `mca_pbmc/bm/bm_ckit/bm_mesench` (local MCA dirs via `load_ref_mca*()`)
  - Added 9 × 3 = 27 Azimuth projection `tar_map` blocks (`proj_tms_bm`, `proj_tms_sp`, `proj_tms_lu`, `proj_tms_ki`, `proj_bm_atlas`, `proj_mca_pb`, `proj_mca_bm`, `proj_mca_ckit`, `proj_mca_mes`) each over stroma/hsc/mnc subsets; outputs as CSV via `save_csv`
  - Added `anno_meta` / `anno_meta_path` `tar_map` targets — extract umap.nn coords + cluster cols + sample/HTO metadata to plain data.frame CSV
  - Added `anno_expr_stroma/hsc/mnc` — SCT expression matrix for all canonical marker genes
  - Added `anno_raw_meta_stroma/hsc/mnc` — pre-filter metadata (umap.nn + QC flags from `metadata_flagged_all`)
  - Added `ref_umap_df_*` per reference — tiny UMAP coord + cell_type data.frame for notebook plotting
  - Added `dim_grid_sct_pre` / `dim_grid_sct_runs` `tar_map` targets — pre-computed dim-reduction grid data (pre: elbow/PCA data frames; runs: UMAP + cluster cols per dims×k combination)
  - Milo targets moved into pipeline: `milo_cluster_col`, `milo_stroma`, `milo_hsc`, `milo_mnc` (format = "qs")
  - Added `anno_meta_cond` `tar_map` target — `add_condition_df()` + rename cluster col to `"cluster"`
  - DE/ORA/module-score targets moved into pipeline: `allvsall_genesets_hsc/mnc/stroma`, `ora_results_hsc/mnc/stroma`, `module_scores_hsc/mnc/stroma`, `de_feature_genes_*`, `de_feature_expr_*`, `biased_feature_expr_hsc/mnc/stroma`
- `src/scrna-seq/scrna-seq_utilities.R`: ~1800 lines of new helper functions added
  - Reference loaders: `fetch_ref_census()`, `seurat_ref_process()`, `load_ref_marrow_atlas()`, `load_ref_mca()` (generic) + four typed wrappers (`load_ref_mca_pbmc/bm/bm_ckit/bm_mesench`), `load_ref_scp1248()`
  - Azimuth label transfer: `build_azimuth_ref()`, `run_azimuth_df()`, `transfer_labels_df_azimuth()`
  - Annotation metadata extractors: `extract_raw_meta()`, `extract_anno_meta()`, `extract_expr_df()`, `extract_ref_umap_df()`
  - Data-frame UMAP helpers: `plot_umap_df()`, `plot_umap_feature_df()`, `plot_umap_grid_sample_dm2_df()`
  - Cluster annotation helpers: `anno_res_col()`, `anno_merge_recluster()`, `anno_merge_proj()`, `anno_load_proj_list()`, `anno_load_ref_df_list()`
  - Resolution and reference plotting: `plot_resolution_umap_grid()`, `plot_feature_umap_panels()`, `plot_ref_umap()`, `plot_ref_proj_umap()`, `plot_query_density_umap()`, `plot_query_label_umap()`, `plot_mapping_score_umap()`, `plot_prediction_score_umap()`, `print_projection_panel()`, `anno_make_type_pal()`
  - ARI + score diagnostics: `compute_ari_table()`, `plot_ari_line()`, `plot_mapping_score_violin()`, `plot_score_by_cluster()`, `plot_score_by_celltype()`, `plot_resolution_heatmap()`
  - Dim-grid helpers: `compute_dim_grid_sct_pre()`, `compute_dim_grid_sct_runs()`, `recluster_resolutions_df()`
  - Milo helpers: `add_condition_df()`, `add_condition_seurat()`, `run_milo()`, `run_milo_contrasts()`
  - DE/ORA helpers: `compute_allvsall_genesets()`, `compute_ora_results()`, `extract_module_scores()`, `extract_feature_expr_long()`
  - Dot plot helper: `plot_dot_df()`
- `notebooks/cluster-annotation.qmd` (NEW): full cluster annotation notebook with marker-based (dot plots + UMAP across resolutions) and reference-based (Azimuth label transfer from 9 references × 3 compartments, ARI line plots, mapping/prediction score violins, resolution heatmaps) sections
- `project_config.yaml`: added `output_table_dir`, `output_figure_dir`, `azimuth_dir`, `tms_marrow_h5ad`; removed `run_key`/`run_label`; added `dim_grid_dims`, `dim_grid_k`, `dim_grid_res_scan` to `processing:` block
- `notebooks/compositional-analysis.qmd`, `notebooks/differential-expression.qmd`, `notebooks/first-analysis-milestone-report.qmd`: updated to load from pipeline targets (data.frames, not Seurat objects) via the new `anno_meta`, `module_scores_*`, `ora_results_*`, `biased_feature_expr_*` targets
- Various notebooks: restructured to use `tar_load` from pipeline targets rather than inline heavy computation

### Why
The previous pipeline exposed only the final `mpn_processed_<subset>` Seurat objects to notebooks, which forced notebooks to either load heavy Seurat objects or rerun expensive computations inline. This session extracted all analysis-relevant quantities (UMAP coords, cluster labels, expression matrices, projection scores, Milo objects, geneset/ORA results) into lightweight targets (plain data.frames or qs-cached Seurat only when needed). The annotation notebook was the primary motivation: systematic label transfer from 9 curated references (TMS Marrow/Spleen/Lung/Kidney, Mouse BM Atlas, MCA PBMC/BM/BM-cKit/BM-Mesenchyme) provides principled cluster annotation support. Biologically: Lung and Kidney TMS references serve as negative controls — clusters labelled by lung/kidney cell types indicate annotation noise rather than true cell identity.

### Bugs / issues
- **`run_label`/`run_key` stamping**: the old `output_dir` target mutated `project_config.yaml` at runtime, which invalidated all downstream targets on every re-run (yaml write changed mtime). Fix: removed stamped run key entirely; `output_table_dir` and `output_figure_dir` are now static paths set manually in yaml.
- **`tar_config_set` missing from `_targets.R`**: targets was not finding the store when sourced from notebooks via `here::here("_targets")` in some contexts. Fix: added `tar_config_set(store = here::here("_targets"))` at top of `_targets.R`.
- **Reference UMAP naming**: `marrow_atlas` uses `refUMAP` embedding name while TMS and MCA references use `umap.nn`. `extract_ref_umap_df()` accepts `umap_col` parameter (defaults to `"refUMAP"`) to handle both.

### Decisions
- **9 references × 3 compartments = 27 projection `tar_map` blocks**: each reference gets its own `tar_map` block rather than a nested loop, which keeps targets names predictable (`proj_<ref>_<subset>`) and allows independent invalidation per reference.
- **Lung + Kidney as negative controls**: TMS Lung and Kidney are biologically unrelated to bone marrow. If a query cluster gets consistently labelled "Type II pneumocyte" or "proximal tubule cell", this signals weak label transfer, not a real cell type — the mapping score violin will show low confidence. Including them guards against over-interpretation of transfer scores.
- **Pre-filter raw metadata targets (`anno_raw_meta_*`)**: notebooks that visualise pre-QC UMAPs (milestone report section 3) were loading full `mpn_raw_processed_*` Seurat objects (>1 GB each). New `extract_raw_meta()` extracts only umap.nn coords + PC1 + QC flags into a lightweight data.frame, replacing the full Seurat load.
- **`allvsall_genesets` / `ora_results` / `module_scores` targets in pipeline**: ORA and module scoring are deterministic given the DE results; moving them into the pipeline means the notebook only reads pre-computed data.frames — faster render, lower memory, more reproducible.
- **`dim_grid_sct_runs` target**: runs SCTransform once per subset then sweeps dims × k combinations. Avoids 6 repeated SCTransform calls (one per grid cell) which would OOM on a 16 GB machine.

### Compartment status
- **Stroma** ✅: full pipeline; annotation targets (recluster + 9 projections + meta/expr extracts) defined; annotation notebook ready for render
- **HSC** ✅: full pipeline; annotation targets defined; notebook ready for render
- **MNC** ✅: full pipeline; annotation targets defined; notebook ready for render
- **HSC+MNC** ⏳: processing targets still present; annotation targets not defined (excluded from this session)

### Next steps
- [ ] Run `tar_make(names = tar_glob_raw("annotation_recluster_*|proj_*|anno_meta_*|anno_expr_*|anno_raw_meta_*|ref_umap_df_*"))` to build annotation targets
- [ ] Render `cluster-annotation.qmd` and review ARI lines + mapping score violins per compartment
- [ ] Select final cluster resolution per compartment based on ARI + biology
- [ ] Assign cell type labels to clusters and update `project_config.yaml → annotation:`
- [ ] Decide whether to keep or drop hsc_mnc combined subset

---

## 2026-04-14 — cluster-annotation.qmd full rewrite + density/mixing sections

### What changed
- `notebooks/cluster-annotation.qmd`: complete rewrite of the reference-based annotation section
  - Replaced single 170-line `results='asis'` loop (with base64/emit_png) with explicit per-reference markdown sections — one `##` heading and explicit `::: {.panel-tabset}` fences per reference, one named chunk per subsection
  - Each of 9 references has a `*-setup` chunk (loads `ref_umap_df_*` + projection data + builds palette) followed by static subsection chunks: Projected query: labels | Projected query: density | Mapping score | Summary
  - Removed "Predicted cell-type" section from all 9 references (redundant with labels tab)
  - Removed "Marker UMAPs" tab from all three compartments in the Canonical markers section
  - Removed `proj-load-all` / `ref-dfs-load` / `ref-palettes` upfront bulk-loading chunks; each reference setup chunk now loads its own data via `load_proj()` helper + `tar_read_raw("ref_umap_df_*")`
  - Violin chunks simplified: `ref_plot_violins()` now returns patchwork directly (no `$plot/$n` list); `fig.height=24` set statically in chunk header
  - x-axis label on cluster violin plots changed from `"Cluster"` to the full column name (e.g. `res_0_5`)
  - Compartment separator labels (bold ggplot void annotation) added between compartment blocks in violin stacks
  - x-axis text size in `plot_mapping_score_violin` increased 8 → 11
- `src/scrna-seq/scrna-seq_utilities.R`: added 8 per-reference section helpers
  - `ref_make_legend()`, `ref_plot_labels()`, `ref_plot_density()`, `ref_plot_celltype()`, `ref_plot_mapping()` — each draws one subsection for one compartment, returns NULL if data unavailable
  - `ref_plot_violins()` — stacks all compartment violin pairs with bold compartment label separators
  - `ref_draw_heatmaps()` — draws ComplexHeatmap for all compartments via `draw()`
  - `ref_plot_ari()` — builds ARI line plot across compartments
  - Updated `plot_score_by_cluster()`: x-label now passes `clust_col` (full column name) instead of `"Cluster"`
  - Updated `ref_plot_violins()`: returns `wrap_plots()` directly instead of `list(plot=, n=)`
- `notebooks/compositional-analysis.qmd`: added new `# Cell density & mixing` section before `# Cluster proportions`
  - Point density UMAP (ggpointdensity, magma palette) faceted by `sample.name`, one tab per compartment
  - Mixing metric UMAP (`Seurat::MixingMetric()` WT vs Het, turbo palette), one tab per compartment
- `notebooks/cluster-annotation.qmd`: added new `# Mixing metric & group-correlated genes` section after reference annotation
  - `Seurat::MixingMetric()` per compartment (WT vs Het on `umap.nn`)
  - `Seurat::GroupCorrelation(var = "mix")` to surface genes correlated with condition segregation; results shown as searchable DT table sorted by `|mix_r|`

### Why
The `results='asis'` loop approach (with `cat()` for all tabset fences and `emit_png`/base64 for plots) was unreadable, caused tab-leak bugs (last tab bleeding into next section), and made violin height control impossible. The rewrite makes the document structure explicit in markdown — each reference is a real section with real headings and fences, each plot is in a named chunk with `fig.width`/`fig.height` in the header. The section-level helper functions in utilities remove all the repetitive R logic from the notebook while keeping the notebook itself as a simple call-per-chunk document.

The mixing metric + group correlation analysis addresses a new biological question: beyond cluster boundaries, which genes track the degree of WT/Het intermixing in UMAP space? `MixingMetric` quantifies per-cell segregation; `GroupCorrelation` then correlates every gene's SCT expression with this score, prioritising genes that co-vary with condition separation.

### Bugs / issues
- **Tab leak (root cause)**: `cat(":::\n\n")` after a plot in `results='asis'` doesn't give Quarto a text node to anchor tabset closure — last tab bleeds into next section. Attempted `cat("&nbsp;\n\n")` workaround, then base64/temp-file approaches — all rejected. Fix: removed `results='asis'` loop entirely; explicit markdown structure has no tab leak.
- **Violin height** in `results='asis'` loop: `knitr::opts_chunk$set()` affects subsequent chunks only, not the current chunk. `knitr::opts_current$set()` was tried but unreliable. Fix: static `fig.height=24` in chunk header (max panels = 6 violins per reference × ~4 in each).
- **`ref_plot_violins` list vs plot**: originally returned `list(plot=, n=)` to communicate panel count for dynamic height — led to `if (!is.null(.v)) { knitr::opts_current$set(...); .v$plot }` in every violin chunk. Simplified to return the plot directly with static height.

### Decisions
- **One named chunk per subsection per reference**: 9 references × 4 subsections × 3 compartments = ~108 chunks. Verbose but completely explicit — no hidden loop logic, freeze cache works correctly per chunk, errors are localised.
- **`*-setup` chunk pattern**: each reference's `.ref`, `.name`, `.dfs`, `.ref_df`, `.pal` assigned in a single `include=FALSE` setup chunk. The 3–5 remaining chunks per subsection are single function calls with no local variables — maximum readability.
- **`load_proj()` helper retained**: loading 3 compartment projections per reference still requires a small loop (3 `tar_read_raw` calls). This is in a helper function defined once, not repeated 9 times in the notebook.
- **Compartment label in violin stack**: rather than splitting into 6 separate chunks (one per violin), a `ggplot() + annotate("text") + theme_void()` spacer panel is inserted as a patchwork element before each compartment pair — one function handles all compartments.

### Compartment status
- **Stroma** ✅: annotation notebook sections complete; density/mixing sections added
- **HSC** ✅: annotation notebook sections complete; density/mixing sections added
- **MNC** ✅: annotation notebook sections complete; density/mixing sections added
- **HSC+MNC** ⏳: excluded from annotation notebook (no projection targets defined)

### Next steps
- [ ] Render `cluster-annotation.qmd` with `-f` flag and verify all 9 reference sections render cleanly
- [ ] Render `compositional-analysis.qmd` and verify density + mixing sections
- [ ] Check `GroupCorrelation` column name — confirm `mix_r` is correct for `var = "mix"`
- [ ] Select final resolution per compartment from ARI curves
- [ ] Assign cell type labels and add to `project_config.yaml → annotation:`

---

## 2026-04-16 14:01 — Machine migration: WSL2 → macOS M5 Pro via Lima + Apptainer

### What changed
- `container/bioinfo.def`: fixed hardcoded `linux-amd64` Quarto download URL — replaced with `dpkg --print-architecture` so the correct `.deb` is fetched on arm64
- `container/python_setup.py`: created placeholder file (checked by `build.sh` but not referenced in `bioinfo.def`)
- `scripts/run_rstudio.sh`: fixed `$HOME` for Lima `.guest` username suffix; moved `RSTUDIO_RUN` to `/tmp` (virtiofs cannot create Unix sockets); added `--auth-minimum-user-id=500` for macOS UIDs below 1000; updated default data path to `data_hub/mpn_10x_camelia`
- `scripts/run_shell.sh`: fixed `$HOME` for Lima `.guest` username; updated default data path
- `scripts/run_targets.sh`: fixed `$HOME` for Lima `.guest` username; updated default data path
- `MIGRATION.md`: fixed Apptainer install instructions (PPA instead of arch-specific `.deb`; Ubuntu universe fallback removed); fixed Quarto arch detection note
- `_targets.yaml`: minor update (user-edited)

### Why
Project migrated from WSL2 Linux machine to macOS Apple M5 Pro. The existing container and scripts were written for x86_64 Linux and needed several fixes to run on arm64 macOS via Lima VM + Apptainer.

### Bugs / issues
- **Quarto `amd64.deb` on arm64**: `bioinfo.def` downloaded `quarto-*-linux-amd64.deb` unconditionally — `dpkg` rejected it with "package architecture (amd64) does not match system (arm64)". Fixed with `QUARTO_ARCH="$(dpkg --print-architecture)"`.
- **Apptainer not in Ubuntu 24.04 universe**: Lima VM uses Ubuntu 24.04; `apptainer` is not in the default universe repo for that release. Fixed by adding the Apptainer PPA (`ppa:apptainer/ppa`) — installs v1.4.5.
- **Lima mount read-only**: `limactl start` with default config mounts `~` as read-only. Required deleting and recreating the VM from a config file (`~/.lima/bioinfo.yaml`) with `writable: true` explicitly set. Lima does not pick up this change on stop/start of an existing instance.
- **Lima `.guest` username suffix**: inside the Lima VM, the Mac username becomes `mohamed_shoeb.guest` (UID 503 conflicts with Linux UID space). `$HOME` resolves to `/home/mohamed_shoeb.guest` — a path that doesn't contain the Mac files. Fixed by overriding `HOME` at the top of all scripts: `HOME="/Users/$(whoami | sed 's/\.guest$//')"`.
- **RStudio socket on virtiofs**: `RSTUDIO_RUN` was bound to a project subdirectory on the virtiofs mount; RStudio Server failed to `chmod` a Unix socket on that filesystem. Moved `RSTUDIO_RUN` to `/tmp/rstudio_run` (local VM fs).
- **RStudio UID < 1000**: Mac UIDs start at 501; RStudio Server's default `auth-minimum-user-id=1000` blocked login. Fixed with `--auth-minimum-user-id=500`.

### Decisions
- **Lima over Podman Desktop**: Lima with Apptainer chosen to preserve the existing Singularity-based workflow (`apptainer exec` / `apptainer shell`) without modification to the core pipeline scripts.
- **`/tmp` for RStudio runtime dir**: virtiofs (used by Lima's `vz` VM type on Apple Silicon) does not support Unix domain sockets. Runtime state (`/var/run/rstudio-server`) must live on the VM's local tmpfs. Session state and config (`RSTUDIO_VAR`, `RSTUDIO_SESSIONS`) remain on the project mount as they only need regular file I/O.
- **PPA for Apptainer**: `ppa:apptainer/ppa` provides v1.4.5 on Ubuntu 24.04 arm64 reliably. Ubuntu universe does not package apptainer for Noble (24.04); the PPA is the canonical install path recommended by the Apptainer project.
- **VM resources**: 14 CPUs / 40 GiB RAM / 200 GiB disk allocated from M5 Pro (18 cores / 48 GB total), leaving headroom for macOS.

### Compartment status
- **Stroma** ⏳: pipeline not yet run on new machine — container just built
- **HSC** ⏳: pipeline not yet run on new machine
- **MNC** ⏳: pipeline not yet run on new machine
- **HSC+MNC** ⏳: pipeline not yet run on new machine

### Next steps
- [ ] Verify RStudio Server launches cleanly at http://localhost:8787
- [ ] Transfer or symlink raw data to `~/data_hub/mpn_10x_camelia`
- [ ] Run `tar_make(names = "config_file")` smoke test inside container
- [ ] Run full pipeline on new machine and verify output matches WSL2 results
- [ ] Commit `bioinfo.yaml` Lima template to repo under `container/bioinfo.lima.yaml` for future reference

---

## 2026-04-25 18:58 — Container build fixes + GHCR push script

### What changed
- `container/bioinfo.def`:
  - Added `/root/.R/Makevars` creation step before `Rscript /opt/renv_setup.R`; sets `CXX11 = g++ -std=c++14`, `CXX11STD = -std=gnu++14`, `CXX14 = g++ -std=c++14`, `CXX14STD = -std=gnu++14` — fixes presto/RcppArmadillo C++14 compilation failure
  - Removed `cellxgene.census` smoke test
  - Removed `Azimuth` smoke test
- `container/renv_setup.R`:
  - Added standalone `miloR` install with `force = TRUE` before batch BiocManager run; retry block after batch
  - Added post-install `requireNamespace` checks for `presto` (with retry + early `stop()`)
  - Added post-install checks for `Signac` and `Azimuth` (with retry + early `stop()`) — then removed both packages entirely
  - Removed `cellxgene.census` install block (tiledb has no Linux aarch64 support)
  - Removed `Azimuth` + `Signac` install blocks (GitHub rate limits + R 4.5 fragility)
  - Added explanatory comments for all exclusions
- `scripts/push_container.sh` (new): pushes `bioinfo_latest.sif` to GHCR as OCI artifact via `apptainer push`; supports `--tag` and `--sif` flags; reads `GHCR_USER` + `GHCR_TOKEN` env vars

### Why
Container build was failing after ~49 min with "Missing packages: miloR, presto, Azimuth". Root causes were three independent issues: (1) presto's `src/Makevars` declares `CXX_STD = CXX11` but the installed RcppArmadillo requires C++14 — R's default `CXX11STD = -std=gnu++11` was too old; (2) miloR silently dropped during the batch BiocManager run; (3) Azimuth/Signac hit GitHub API rate limits (unauthenticated: 60 req/hr). The container is built on Lima (aarch64) but deployed on an x86_64 cluster via GHCR ORAS, so a push script was needed.

### Bugs / issues
- **C++14 compile error (presto):** `RcppArmadillo` ≥ 0.12 enforces C++14 at compile time; presto declares `CXX_STD = CXX11` in its `src/Makevars`, so R used `CXX11STD = -std=gnu++11`. Fixed by overriding `CXX11STD = -std=gnu++14` in `/root/.R/Makevars`.
- **miloR silently missing:** BiocManager batch installs can silently skip packages when transitive dependencies fail earlier in the run. Fixed by installing miloR standalone first with `force = TRUE`.
- **GitHub rate limit (Azimuth):** Without `GITHUB_PAT` properly forwarded, `remotes::install_github()` hits the 60 req/hr unauthenticated limit across the many GitHub installs. Root cause of previous failures was passing `$GITHUB_PAT` (empty) instead of `$APPTAINER_GITHUB_PAT`. Azimuth and Signac additionally dropped because they are fragile on R 4.5 and not critical for the pipeline.
- **cellxgene.census / tiledb:** tiledb has no Linux aarch64 binary — fails on Lima VM. Dropped from container; install via host renv on x86_64 if needed.

### Decisions
- **Drop Azimuth + Signac from container:** Azimuth label transfer results are already cached as CSV targets in `_targets.R`; the package itself is not needed at pipeline runtime. Can be installed via host renv on the x86_64 server if re-running annotation targets.
- **Drop cellxgene.census from container:** Used only for loading TMS reference data, which is already pre-processed and cached as local h5ad/RDS files. Architecture mismatch (aarch64 build, amd64 deploy) makes tiledb unbuildable in Lima.
- **ORAS push script over GitHub Actions:** Simpler for now — build locally in Lima, push manually. GHA approach deferred until the container definition stabilises.

### Compartment status
- **Stroma** ✅: pipeline targets complete; annotation notebooks functional
- **HSC** ✅: pipeline targets complete; CellRank + SCENIC done; Milo per-cell scores added
- **MNC** ✅: pipeline targets complete
- **HSC+MNC** ✅: combined targets complete

### Next steps
- [ ] Push `bioinfo_latest.sif` to GHCR: `GHCR_TOKEN=<pat> ./scripts/push_container.sh`
- [ ] Pull on cluster: `apptainer pull oras://ghcr.io/mohamedrefaat92/bioinfo:latest`
- [ ] Install Azimuth via host renv on cluster if annotation re-runs are needed
- [ ] Consider GitHub Actions build for future container versions (x86_64 native, no cross-arch issues)

---

## 2026-04-25 19:15 — run_velocyto.sh updated for CellRanger Multi layout

### What changed
- `scripts/run_velocyto.sh`:
  - Rewrote for CellRanger Multi output layout: `Rx_*/outs/per_sample_outs/*/count/` instead of `Rx_*/outs/possorted_genome_bam.bam`
  - BAM is now `sample_alignments.bam` (Multi) instead of `possorted_genome_bam.bam` (count)
  - Switched `velocyto run10x` → `velocyto run` with explicit `--bcfile` and `--sampleid` flags
  - Bind mounts changed: `BIND_DATA_HUB` → `BIND_DATA` (bound to `/data`) + `BIND_GENOMES` (bound to `/genomes`)
  - GTF auto-detection from CellRanger `_log` removed; `BIND_GENOMES` must point to a dir containing `GRCm39/genes/genes.gtf.gz`
  - Added gzip decompression of GTF to `/tmp/genes.gtf` (velocyto does not accept `.gz`)
  - Added explicit `BIND_DATA` / `BIND_GENOMES` required-env checks with clear error messages
  - Skip logic now checks both BAM and barcodes existence per sample
  - Removed Lima-specific `HOME` override (script now targets cluster directly)

### Why
The CellRanger run on the cluster used `cellranger multi` (not `cellranger count`), which produces a different output directory structure. The old script assumed `run10x` layout with a single BAM at `outs/possorted_genome_bam.bam`; the Multi layout has per-sample BAMs at `outs/per_sample_outs/<sample>/count/sample_alignments.bam`.

### Bugs / issues
None in this change — preventive rewrite ahead of running on cluster.

### Decisions
- **`velocyto run` over `run10x`:** `run10x` is a convenience wrapper that assumes the old count layout; `run` gives full control over BAM + barcodes paths, required for Multi output.
- **GTF decompression to /tmp:** velocyto's C extension does not accept gzipped GTF. Decompressing at runtime avoids storing a duplicate uncompressed copy.
- **Separate `BIND_GENOMES` mount:** Keeps genome references decoupled from sample data, matching the cluster's storage layout.

### Compartment status
- **Stroma** ✅: pipeline complete
- **HSC** ✅: pipeline complete; velocity inputs will be generated by this script
- **MNC** ✅: pipeline complete
- **HSC+MNC** ✅: combined targets complete

### Next steps
- [ ] Push `bioinfo_latest.sif` to GHCR: `GHCR_TOKEN=<pat> ./scripts/push_container.sh`
- [ ] Pull container on cluster: `apptainer pull oras://ghcr.io/mohamedrefaat92/bioinfo:latest`
- [ ] Run velocyto on cluster: `BIND_DATA=... BIND_GENOMES=... ./scripts/run_velocyto.sh`
- [ ] Load `.loom` files in `notebooks/velocity_run.py` and re-run CellRank with spliced/unspliced

---

## 2026-05-28 — Mpig6b expression notebook

### What changed
- `notebooks/mpig6b-expression.qmd` (new): exploratory notebook for *Mpig6b* expression across HSC, MNC, and Stroma compartments in raw and processed Seurat objects.

### Why
*Mpig6b* is a low-abundance gene (platelet CLEC-2 / ITAM signalling receptor) of biological interest in the MPN context. It is expected to be dropped by the pipeline's `min_cells_per_gene = 20` gene filter and potentially also by SCTransform's internal variance filter, so it cannot be assessed from the standard processed (`SCT`) objects alone. The notebook addresses this by querying the raw `RNA` assay (pre-QC, pre-gene-filter) to maximise sensitivity.

### Bugs / issues
- **`pass_cell` not found in processed `*_proc_df`:** `anno_meta_*` is extracted from `mpn_processed_*` via `extract_anno_meta()`, which does not carry `pass_cell` — that flag only lives in the raw-side metadata. Since all cells in the processed objects already passed QC, the fix was to apply the `pass_cell` filter only when the column is present in the data frame. Fixed with `{ if ("pass_cell" %in% colnames(.)) dplyr::filter(., pass_cell == TRUE) else . }` in `violin_condition()` and `violin_cluster()`.

### Decisions
- **Raw counts via `mpn_raw_*` + UMAP from `anno_raw_meta_*`:** `mpn_raw_*` contains the unfiltered RNA assay needed for sensitive detection; `anno_raw_meta_*` (extracted from `mpn_raw_processed_*`) provides the pre-QC UMAP without loading the full Seurat twice. This avoids OOM on the 16 GB machine.
- **Conditional `eval=` on SCT sections:** Chunk option `eval = <compartment>_in_proc` suppresses the processed-expression panels entirely (rather than printing an empty plot) when the gene is absent from SCT — cleaner output for the expected case.
- **Load → extract → `rm` / `gc()` pattern:** Follows the established low-memory notebook convention; Seurat objects are discarded immediately after expression extraction.
- **Separated zero-count cells on UMAP:** Grey fixed-colour layer for zero-count cells rendered first, coloured layer on top — avoids viridis scale being dominated by zeros while keeping spatial context.

### Compartment status
- **Stroma** ✅: no change
- **HSC** ✅: no change
- **MNC** ✅: no change
- **HSC+MNC** ✅: no change

### Next steps
- [ ] Render `notebooks/mpig6b-expression.qmd` inside container and verify detection tables
- [ ] If Mpig6b is detected above threshold in any compartment, consider adding it to `canonical_markers_vec` in `annotation_config` so `anno_expr_*` targets pick it up
- [ ] Pull container on cluster: `apptainer pull oras://ghcr.io/mohamedrefaat92/bioinfo:latest`
- [ ] Run velocyto on cluster: `BIND_DATA=... BIND_GENOMES=... ./scripts/run_velocyto.sh`

---

## 2026-10-27 — hdWGCNA pipeline integration, per-condition transition matrices, trajectory refinement, genomic annotation tools

### What changed
- `_targets.R`: added `sex_degs`, `hsc_pairwise_comparisons`, `mnc_pairwise_comparisons`, `stroma_pairwise_comparisons`, `cytotrace_hsc_wt`, `hem_markers_claude_df`; external data file targets `yoshida_2020_ocr_csv` (ImmGen ATAC), `hpa_secretome_tsv`, `naba_matrisome_mm_csv`; SCENIC targets `scenic_agg_column` + `scenic_meta_hsc`; merged multi-compartment objects `mpn_all`, `mpn_wt_all`, `mpn_mch_neg_all`, `mpn_mch_pos_all`; hdWGCNA targets (`hdwgcna_hsc_modules`, `hdwgcna_hsc_hvg`, `hdwgcna_hsc_soft_power`, `hdwgcna_hsc_metacells`, `hdwgcna_hsc_network`) — +981 lines overall
- `src/scrna-seq/scrna-seq_utilities.R`: added `find_delegate_pairwise()`, `plot_nhood_graph_groups()`, genomic range helpers (`build_protein_coding_genes_gr`, `build_protein_coding_promoters_gr`, `link_ranges_to_genes`, `parse_yoshida_orig_links`), pseudotime binning helpers (`bin_expr_along_pseudotime`, `plot_lineage_heatmap`, `expr_scale_pal_gradientn`), `compute_magic()` — +683 lines
- `notebooks/hdwgcna-hsc.qmd`: major refactor — removed all inline hdWGCNA computation (SetupForWGCNA, MetacellsByGroups, NormalizeMetacells, SetDatExpr, TestSoftPowers, ConstructNetwork, ScaleData, ModuleEigengenes, ModuleConnectivity); now loads from pipeline targets; TOM written to `output/hdwgcna/` by pipeline
- `notebooks/transition-matrices.qmd`: major expansion — per-condition (WT / mCh− / mCh+) cluster transition matrices for all 6 matrix types; CellRank cluster × cluster Markov transition matrix loaded from `output/cellrank/cytotrace/tables/cluster_transition_matrix.csv`; hdWGCNA metacell igraph analysis (majority-vote cluster and condition per metacell); clickable image lightbox overlay added
- `notebooks/trajectory-slingshot.qmd`: refactored to use `mpn_processed_hsc_fullcc` (full cell-cycle object); dedicated config block for `CLUS_COL`, `N_PCA_DIMS`, `EXCLUDE_CLUSTER`, `START_CLUSTER`, `END_CLUSTER`; cells filtered before PCA embedding passed to Slingshot; added `ggridges` dependency
- `notebooks/graph-analysis.qmd`: updated source object to `mpn_processed_hsc_fullcc`; simplified MST layout to UMAP centroid positions (removed Reingold-Tilford tree layout)
- `notebooks/scenic-hsc.qmd`: added Cell UMAP section (cluster + condition); renamed pipeline targets to `_full` suffix; alias shim preserves downstream code
- `notebooks/compositional-analysis.qmd`: Milo DA nhood graphs now use `scale_fill_gradient2` (blue–white–red, limits ±10); mCh neg vs mCh pos neighbourhood LFC scatter upgraded to `ggpointdensity`
- `notebooks/cellrank_run.py`: added `dpt` as third `INPUT_TYPE` (alongside `cytotrace` and `velocity`); configurable `CLUSTER_RES` env var for flexible cluster resolution in Python; nullable string dtype fix for pandas/anndata compatibility
- `notebooks/differential-expression.qmd`: major content trim (migrated redundant sections to dedicated notebooks)
- `container/bioinfo.def`, `container/renv_setup.R`: minor dependency updates
- All notebooks: `embed-resources: true` → `embed-resources: false` (external resources; reduces per-file size for notebooks rendered outside Quarto book context)
- `daily_notebooks/2026_06_01.qmd` – `2026_08_03.qmd`: new daily notebooks exploring p53 pathway, Wang 2021 developmental reference projection, ChIP-Atlas TF binding enrichment, STEAM mm10 enhancer predictions, JASPAR motif analysis, sccore graph diffusion label propagation

### Why
The hdWGCNA computation (metacell construction, TOM, eigengenes) is expensive (~20–30 min) and reproducibility required moving it into `targets` cache. The transition matrices notebook was extended to allow direct comparison of cluster connectivity profiles across conditions (WT vs Het), giving a quantitative view of how the MPN mutation alters cluster-to-cluster Markov flux and transcriptional similarity. The slingshot and graph notebooks needed to use `mpn_processed_hsc_fullcc` (the high-res cell-cycle-aware object) rather than the lower-res `mpn_processed_hsc` to match the cluster resolution used in all other HSC analyses. Exploratory daily notebooks tracked ongoing biological interest in TF binding, regulatory genomics, and signal pathway enrichment.

### Bugs / issues
- **`label_condition` setNames crash (hdwgcna-hsc):** Original code called `setNames(., names(mpn_processed_hsc_fullcc$sample.name))` outside of any mutation context, crashing when `mpn_processed_hsc_fullcc` was not loaded at the top level. Fixed by moving `setNames` call to after the pipeline object is renamed to `srt`.
- **CellRank nullable string dtype (cellrank_run.py):** `anndataR`-written h5ad files store string columns with pandas `StringDtype` / pyarrow backend; CellRank's internals expected plain `object` dtype, causing downstream attribute errors. Added `_fix_nullable_strings()` helper to cast all string columns to `object` at load time.
- **Transition matrices CELLRANK_DIR path:** Notebook was hardcoded to `cellrank/velocity_dynamical`; actual CytoTRACE run outputs land at `cellrank/cytotrace`. Fixed path constant.

### Decisions
- **hdWGCNA targets in pipeline:** Keeps the TOM file deterministic (reproducible soft power, fixed seeds in `tar_target`) and avoids re-running the 20-min construction on every notebook open. `output/hdwgcna/` is the TOM directory, created by `ConstructNetwork` inside the target.
- **`mpn_all` / `mpn_wt_all` / `mpn_mch_neg_all` / `mpn_mch_pos_all`:** Cross-compartment merged objects needed for any analysis requiring simultaneous access to HSC + MNC + Stroma cells (e.g. shared marker comparison, SCENIC cross-compartment regulon activity). Merged with `merge.data = FALSE` to avoid storing the full data layer twice. mCh+ Stroma special-cased to include all Het Stroma cells (no mCh+ sorting for Stroma).
- **`embed-resources: false`:** Notebook HTML files were approaching 50–100 MB each with embedded resources; external mode keeps them renderable while UMAP PNGs live separately.
- **DPT as third CellRank direction:** Diffusion pseudotime (from `dpt_hsc.ipynb`) provides a topology-agnostic direction signal as a counterpart to CytoTRACE and velocity; keeping all three modes in `cellrank_run.py` allows direct comparison without code duplication.

### Compartment status
- **Stroma** ✅: pipeline complete; `mpn_all` merged object available for cross-compartment analyses
- **HSC** ✅: hdWGCNA cached in pipeline; transition matrices, slingshot, SCENIC, CellRank all updated
- **MNC** ✅: pipeline complete; `mpn_all` merged object available
- **HSC+MNC** ✅: combined targets complete

### Next steps
- [ ] Run `tar_make()` to build new targets (`hdwgcna_hsc_*`, `mpn_all`, `scenic_meta_hsc`, pairwise comparisons)
- [ ] Render `notebooks/hdwgcna-hsc.qmd` to verify pipeline-loaded object matches prior interactive run
- [ ] Run `cellrank_run.py` with `INPUT_TYPE=dpt` to generate DPT-based fate probabilities
- [ ] Integrate ChIP-Atlas TF binding data into a genomic annotation notebook (promoter overlap → target gene lists)
- [ ] Add STEAM enhancer overlap to HSC differential gene loci for regulatory context
