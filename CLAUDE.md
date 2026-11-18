# CLAUDE.md — MPLW515L mouse bone marrow scRNA-seq

## Project overview

**Question:** How does heterozygous MPL mutation alter haematopoietic compartment composition and cell identity in a mouse model?
**Focus:** scRNA-seq characterisation of bone marrow compartments (Stroma, HSC, MNC) from WT, Het mCh-, and Het mCh+ mice sorted by 10x Chromium. Cluster marker discovery and differential expression between genotype groups per compartment.
**Organism:** Mouse. Reads aligned to GRCm39 + GENCODE vM33 with the mCherry and DRE reporters added (Cell Ranger 10.1.0; the first processing used 9.0.0); some annotation resources (EnsDb.Mmusculus.v79, BSgenome mm10, ImmGen ATAC peaks, cisTarget databases) are mm10.
**Collaborator:** Camélia (Psaila lab)

---

## Scientific background

MPN (myeloproliferative neoplasm) mouse models with heterozygous *Mpl* mutations show altered megakaryocyte and HSC compartments. This project characterises the transcriptional landscape of:

- **Stroma** — bone marrow stromal cells (MSC, fibroblast, osteolineage, pericyte, endothelial)
- **HSC** — haematopoietic stem and progenitor cells
- **MNC** — mononuclear cells (monocyte, macrophage, T/B cells, erythroid, MK lineage)

Samples: tubes Rx_1–Rx_8 (Stroma Rx_1–2, HSC Rx_3–5, MNC Rx_6–8), including WT, Het mCh-, and Het mCh+ conditions.
HTO demultiplexing (TotalSeq-B) tags samples within each tube; deMULTIplex2 assigns singlet/multiplet/negative classes.

---

## Project structure

The project follows the hyprid-proj template: `tools/README.md` explains the environment (conda `./env` for R and system libraries, uv `.venv` for Python, renv for R packages), `data/` + the data catalogue, the targets pipeline and its reports. `README.md` describes this project's data and the order of the steps.

```
_targets.R                         — the pipeline (~400 analysis targets + 64 reports)
metadata/project_config.yaml       — analysis parameters (QC, processing, subclusters, markers, colours)
metadata/samples.csv               — sample sheet: Tube → sample_name_short / sample_name_long
metadata/cellranger/               — Cell Ranger multi inputs per tube (libraries, hashtag feature references)
src/R/scrna-seq_utilities.R        — all analysis functions (sourced by the pipeline and by many reports)
src/python/analysis/               — Python analyses run as cluster jobs (scVelo, CellRank, PAGA, pySCENIC, PerturbGen, CellOracle)
notebooks/<report>/<report>.qmd    — the reports, rendered by the pipeline (tar_quarto) into output/reports/latest/
notebooks/exploration/             — notebooks run by hand (daily notes, Jupyter notebooks, metacells)
notebooks/archive/                 — reports that read targets the pipeline no longer has (README lists them)
scripts/slurm/                     — cluster jobs: Cell Ranger, setup, pipeline, velocity, analysis
scripts/sh/                        — linking inputs, downloading cisTarget resources
scripts/containers/                — Apptainer recipes (Cell Ranger, pySCENIC, velocyto)
data/                              — links to inputs outside the project (Cell Ranger results, references)
output/                            — everything generated: _targets/, reports/, objects/, tables/, containers/ (the .sif images), logs/
```

---

## Running

```bash
sbatch scripts/slurm/pipeline/run-pipeline.sh       # the pipeline and its reports, as a cluster job
env/bin/Rscript -e 'targets::tar_make(names = targets::tar_glob_raw("*_stroma"))'   # a subset
env/bin/Rscript -e 'targets::tar_outdated()'        # what would run
env/bin/Rscript -e 'targets::tar_errored()'         # what failed (error = "continue": others keep going)
```

The Python analyses run as jobs (`scripts/slurm/velocity/`, `scripts/slurm/analysis/`) after the pipeline has written `output/objects/hsc.h5ad`; the reports that read their results render on the next pipeline run.

---

## Pipeline targets — key groups

| Group | Targets |
|---|---|
| Config | `config` (read once at the top of `_targets.R`; targets use `!!config$<section>$<key>`); for reports: `path_params`, `qc_params`, `proc_params`, `grid_params`, `annotation_config`, `colors_vec` |
| Raw data | `tube_sample_map_file` (metadata/samples.csv), `cellranger_matrix_dirs` (data/cellranger_Rx_<n>), `mpn_per_sample`, `mpn_raw_all` |
| QC per compartment (`tar_map` stroma/hsc/mnc) | `demux_res_<x>`, `mpn_raw_<x>`, `mpn_raw_processed_<x>` |
| Processing per compartment | `pass_gene_<x>`, `mpn_clean_<x>`, `mpn_processed_<x>`, `mpn_processed_data_<x>`, `pca_df_<x>` |
| Cell-cycle-regressed HSC | `mpn_processed_hsc_fullcc` (the HSC object most analyses use), `_diffcc`, per condition `_wt/_mchneg/_mchpos_fullcc` |
| DE | `cluster_de_each_vs_rest_*`, `cluster_de_all_vs_all_*`, `condition_de_per_cluster_*`, `scatter_df_*`, `genesets_hsc`, `ora_results_hsc` |
| Annotation | `proj_*` (Azimuth and cross-species projections × compartment), `annotation_recluster_*`, `anno_meta_*` |
| Other analyses | Milo (`milo_*`), miloDE, hdWGCNA (`hdwgcna_hsc_*`), Moran's I, SCENT, CytoTRACE2, SuperCell, MAGIC, PHATE/FA2, SCENIC aggregation (`scenic_*`), gene sets (`*_genesets`, `ucell_scores_*`) |
| Exports | `hsc_h5ad_path` (output/objects/hsc.h5ad, for the Python jobs), `mm_hs_ortho_long_path`, CSV `*_path` targets in output/tables/ |
| Reports | `report_<name>` (one per notebooks/<name>/) |

---

## Conventions

- Parameters come from `metadata/project_config.yaml`, read once into `config` at the top of `_targets.R`. Targets inject them with `!!config$<section>$<key>`, so the value is part of the command and changing one parameter re-runs only the targets that use it. Reports use the section targets (`tar_load(qc_params)` etc.). Never hard-code paths.
- Inputs are file targets on `data/` (links) and `metadata/`; outputs go under `output/` (`path_params`: `output_table_dir`, `objects_dir`).
- Heavy Seurat objects use `format = "qs"`; plot outputs are files, not ggplot objects in the store.
- `DietSeurat` after full processing — downstream targets use `mpn_processed_data` (slim), not `mpn_processed`.
- Reports read results with `tar_read()` / `tar_load()` (tracked); files they read outside targets are listed with their `quarto_report()` line in `_targets.R`.
- R packages come from renv (`scripts/R/install-r-packages.R`), Python packages from uv (`pyproject.toml`); never install from inside a notebook.

---

## Key parameters (metadata/project_config.yaml)

| Parameter | Value |
|-----------|-------|
| QC: mito max | 7% |
| QC: HB max | 20% |
| QC: per compartment | nFeature and ribosomal % ranges in `compartment_qc` |
| QC: min cells per gene | 20 (with ≥ 3 counts) |
| Processing: ref assay | SCT |
| Processing: n features | 3000 |
| Processing: n PCA dims | 50 |
| Processing: cluster resolutions | 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.2, 1.3, 1.5 |
| Processing: k neighbors | 30 |
| Processing: n UMAP neighbors | 30 |
| Marker resolution | 0.1 |
| Marker min.pct / logfc.threshold / max cells per ident | 0.25 / 0.25 / 500 |
| Exclude reporters | mCherry, DRE |

---

## Cluster column naming convention

After `process_all()`, cluster columns follow the pattern:
`sNN_pca.<reduction>_res.<resolution>`

Grep pattern in targets: `grep(paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", marker_res), "$"), colnames(srt@meta.data), value = TRUE)[1]`

---

## Known technical issues / gotchas

### Memory
- The pipeline runs as a cluster job (`scripts/slurm/pipeline/run-pipeline.sh`, 256 GB). The load → extract → `rm()` / `gc()` pattern in the reports dates from a 16 GB laptop and still keeps them lean.
- Set `options(future.globals.maxSize = 8 * 1024^3)` before SCTransform / DElegate.

### presto for FindAllMarkers
- Without presto: ~47 min for 16 clusters. With presto: <1 min. Installed from CRAN by renv (`_dependencies.R` keeps it in renv.lock).

### deMULTIplex2 assignment
- Use `srt@meta.data$dm2_class <- assign_vec[match(colnames(srt), names(assign_vec))]` — direct slot assignment; `$<-` dispatch fails with NAs.
- `demultiplexTags()` writes `<plot.path>/_<timestamp>_assignment.pdf` whenever `plot.umap != "none"`, even with `plot.diagnostics = FALSE`; `plot.path` defaults to `getwd()` (the project root) and `NULL` writes to `/`. `demux_per_rx()` keeps `plot.umap = "umi"` (the result's UMAP needs it) and writes to a temporary folder it deletes; the plots are redrawn from `demux_res_<x>` in `notebooks/demultiplex-diagnostics/` (`dm2_summary_plot()`, `dm2_diagnostic_plots()`).

### Seurat v5 + JoinLayers
- After `merge()`, layers are split per-sample. Call `SeuratObject::JoinLayers()` before downstream steps.
- Use `SeuratObject::` prefix — not exported from `Seurat` namespace.

### DimPlot memory leak
- DimPlot/FeaturePlot ggplots capture calling environment (800 MB–2 GB).
- Mitigation: save to PDF (`format = "file"`); no ggplot objects in targets cache.

### HB gene pattern
- Mouse: `^Hb[^pes]` (excludes Hbp1 transcription factor, Hbs1l, Hbegf).

### HTO CITE-seq loading
- `Read10X()` returns a list when multiple feature types present.
- Extract with: `if (is.list(mat)) mat[["Gene Expression"]] else mat`

### Quarto reads every file in its project
- `_quarto.yml` limits the project to `notebooks/*/*.qmd` minus `exploration/` and `archive/`: a broken front matter (or a Python notebook whose kernel isn't installed) anywhere in the project would otherwise stop every report.

### Quarto cache on GPFS
- Quarto's SASS cache is a SQLite file under `$XDG_CACHE_HOME` (default `~/.cache`); SQLite can't open it on GPFS → `ERROR: unable to open database file: .../sass.kv`. The project's `_environment` sets `XDG_CACHE_HOME=/tmp/${USER}/xdg-cache` for every Quarto render (pipeline, Positron preview, CLI). renv and uv caches are unaffected (`RENV_PATHS_CACHE`, `UV_CACHE_DIR` set explicitly).

### Python from R
- phateR and Rmagic use reticulate; `_targets.R` sets `RETICULATE_PYTHON` to the project's `.venv`, which reaches the reports too. Notebooks with Python chunks call `reticulate::use_virtualenv(".venv")`.

---

## Project development history

### 2026-01 — Pipeline bootstrap

- Scaffold from Shoeb2025_CG-AMKL; adapted for MPN 10x multi-compartment design.
- Implemented `read_10x_sample()`, `merge_seurat_list()`, `flag_qc_metadata()`, `process_all()`.
- Added `compartment` column derived from `sample.name` for Stroma/HSC/MNC subsetting.
- deMULTIplex2 HTO demultiplexing with `dm2_class` metadata column.
- scDblFinder doublet scoring integrated into `flag_qc_metadata()`.
- Quarto book with 6 chapters; all figures saved as PDF/PNG (no Seurat objects in cache).

### 2026-03 — Restructure: harmony removed, per-compartment QC

- Removed harmony integration; pipeline now PCA-only.
- Refactored to per-compartment QC pipeline: `tar_map` over stroma/hsc/mnc for both QC and processing stages.
- Added `mpn_raw_processed_<x>` pre-filter dim-reduction targets for QC notebook visualisation.
- Added DElegate per-cluster DE (`delegate_samples_markers_<x>`) for genotype comparison.
- Added `canonical_markers_path`, `bubble_plot_path`, `marker_feature_plots_path` per subset.
- Added `mnc_qc_curves_clusters` for MNC complexity-curve characterisation.
- HSC+MNC combined subset added as manual targets (not tar_map row).
- `DietSeurat` step added: `mpn_processed_data_<x>` slim objects for downstream targets.
- Cluster column grep updated to `sNN_pca\..*_res\.` pattern (PCA-based, no harmony).
- `cluster_resolutions` changed to scalar `0.1` (single resolution for all subsets).

### 2026-03-31 — Notebook render output restructured

- Rendered book moved out of `notebooks/` — now lands at `analysis_book/<timestamp>/` (timestamped per render).
- `render_notebooks.sh` updated to pass `--output-dir` dynamically using `date '+%Y-%m-%d_%H-%M'`; no config update needed before rendering.
- `notebooks/.gitignore` and `notebooks/_targets.yaml` removed (redundant).
- Root `.gitignore` updated: `notebooks/_book/` → `analysis_book/`.
- Stale render artifacts (`_book/`, `*_files/`, `*_cache/`) deleted from `notebooks/`.

### 2026-04-01 — Notebook QC/DE visualisations + PCA distribution targets

- Added UMAP visualisations of `sample.name`, `scDblFinder.class`, `scDblFinder.score` to pre-QC (section 3) and post-QC dim-reduction (section 4.1) in `first-analysis-milestone-report.qmd`.
- Added `sample.name × dm2_class` UMAP grid (plain + manual highlight) for all subsets, pre- and post-QC.
- Added `plot_umap_grid_sample_dm2()` to utilities; replaced `gghighlight` with a manual two-layer facet approach (ggplot2 version incompatibility).
- Restructured `# Markers` into four subsections: `each_vs_rest`, `all_vs_all`, `mCh- vs mCh+`, `mCh- vs WT`.
- Added `find_delegate_markers_all_vs_all()`, `find_delegate_markers_mch()`, `find_delegate_markers_mch_vs_wt()` to utilities; `dm2_class` used as `replicate_column` for mCh comparisons.
- Added `delegate_allvsall_markers`, `delegate_mch_markers`, `delegate_mch_vs_wt_markers` targets in `tar_map`.
- Added `pca_df` target (plain data.frame, no Seurat object); added `extract_pca_df()`, `plot_pc_parcoor()`, `plot_pc_umap_grid()` to utilities; rewrote `pca-cell-distribution.qmd` to load from targets.
- `render_notebooks.sh`: timestamp changed to day-only (`YYYY-MM-DD`); optional `--label`/`-l` flag added.
- Enforced `tar_load` + `rm()`/`gc()` per chunk throughout the milestone report to prevent OOM kills on 16 GB machine.

### 2026-04-03 — Compositional analysis notebook + sub-clustering pipeline

- Added `notebooks/compositional-analysis.qmd`: cluster proportion bubble plots (sample + sample×HTO), row-scaled pheatmap with % labels, and Milo DA (miloR) for Stroma/HSC/MNC; two top-level sections (*Cluster proportions* / *Milo DA*), each tabbed by compartment.
- Added `notebooks/subcluster-exploration.qmd`: per-sub-cluster UMAP + pie chart + per-sample silhouettes on both sub-cluster and parent embeddings; consistent `sample_pal`; tabs by compartment.
- Added `subclusters:` to `project_config.yaml`: 17 entries (hsc_c0/2/3/4/6, mnc_c0/1/2/3/4/6/7/8, stroma_c0/1/2/3); `marker_*` params moved back to `processing:`.
- Added sub-clustering `tar_map` to `_targets.R`: reads `subclusters:` yaml at parse time; produces `mpn_subclust_<name>` (full reprocessed Seurat) per entry — no markers/DE targets.
- Added `subset_clusters()` to utilities: greps cluster column, subsets cells, adds `parent_cluster` metadata.
- Added `miloR` to `container/renv_setup.R` and `bioinfo.def` smoke tests.
- Milo uses `dm2_class` HTO tags as within-condition replicates (consistent with pseudobulk DE approach) to give edgeR enough degrees of freedom when only one tube exists per condition.

### 2026-04-02 — Notebook visualisation overhaul + dim-reduction workflow fix

- Added scatter plot (Het mCh− vs WT × Het mCh+ vs WT log FC) and per-cluster volcano plots to `first-analysis-milestone-report.qmd` section 5.2 all_vs_all HSC; shape 21 / white border / alpha 0.5; size = max(rate1, rate2); quadrant-based top-N labelling.
- Fixed `inner_join` scatter bug: was filtering `sig` before joining (implicit `&`); changed to join full tables then filter `sig_neg | sig_pos`.
- Rewrote `qc-exploration.qmd`: sections 1–2 now match milestone report (jitter+violin+errorbars, feature scatter, MNC double-curve); added section 2.5 doublet/multiplet heatmap; sweep `fig-height` → 12; removed `hsc_mnc` from UMAP section.
- `dim-reduction-and-clustering-SCT.qmd`: Section 1 loads cached `mpn_processed_data_*` (pipeline-identical); grid loop now calls `process_all()` directly — was manually reimplementing with PCA-based `RunUMAP` instead of NN-graph-based, causing UMAP mismatch.
- `dim-reduction-and-clustering-RNA.qmd`: same structural fix; `preprocess_rna()` bug fixed (`capture.output` returned NULL instead of Seurat object).
- Both dim-reduction notebooks: `umap-sample-tabs` fig size → 15×10 for 3×2 grid.

### 2026-04-03 — DE notebook + Milo pairwise contrasts + scatter diagonal lines

- Created `notebooks/differential-expression.qmd` (Supplementary): full markers section extracted from milestone report; `each_vs_rest`, `all_vs_all`, `mCh- vs mCh+`, `mCh- vs WT` subsections; static + plotly interactive volcano (per cluster) and collective scatter (all clusters).
- Removed entire `# Markers` section from `first-analysis-milestone-report.qmd` (~900 lines); cleared stale freeze cache.
- Fixed Milo contrasts: dropped `model.contrasts` (returned model coefficients, not fold-changes); replaced with two independent two-group runs via `.run_milo_pair()` + `run_milo_contrasts()` wrapper.
- `plot_milo_contrasts()`: 2-row layout (beeswarm | beeswarm / volcano | volcano); heatmap style updated (`#2166ac`/`#b2182b`, no borders, condition annotation bar).
- Diagonal reference lines on scatter: `y=x` (shared Het effect), `y=x±1` (±1 logFC band) in teal; replicated in plotly `shapes`.
- Interactive plotly via `htmltools::tagList(lapply(...))` — the only reliable pattern for multiple widgets in Quarto.
- `top_n_diagonal()` and `repel_by_diagonal_half()` helpers added to shared `allvsall-hsc-helpers` chunk.

### 2026-04-04 — ORA replacing GSEA, bubble plots, geneset membership tables, biased gene expression UMAPs

- Replaced `fgsea::fgsea()` with `fgsea::fora()` (ORA / Fisher's exact test) for HSC, MNC, Stroma; consistent background = all significant scatter genes per compartment.
- Removed NES heatmaps; added ORA bubble plots (fill = log2 foldEnrichment, size = overlap/size, black stroke = padj < 0.05).
- Added geneset membership tables: summary (n_genes per geneset) + full membership (one row per gene), sorted by foldEnrichment.
- Added biased gene expression UMAP section per compartment: off-diagonal genes plotted individually per geneset with `patchwork::plot_annotation` title, split by `sample.name`.
- Added Stroma ORA + module score + biased expression (previously missing); genesets defined as `up_Het`/`down_Het` from volcano.
- `viridis_option` parameter threaded through `lean_dimplot` → `plot_umap_feature` → `plot_umap_feature_ext` in utilities.
- DT widget fix in `results="asis"` loops: `print(htmltools::tagList(DT::datatable(...)))`.
- Cairo size limit fix: iterate per geneset (≤20 genes) rather than all biased features in one call.

### 2026-04-13/14 — cluster-annotation.qmd full rewrite + density/mixing sections

- Rewrote reference-based annotation section: replaced single `results='asis'` loop (170 lines, base64/emit_png, tab-leak bugs) with explicit per-reference `##` sections and named chunks — one chunk per subsection per compartment.
- Added 8 per-reference section helpers to `scrna-seq_utilities.R`: `ref_plot_labels/density/celltype/mapping/violins/heatmaps/ari` + `ref_make_legend`; each returns NULL if data unavailable.
- Removed "Predicted cell-type" section from all 9 references; removed "Marker UMAPs" tab from Canonical markers (HSC/MNC/Stroma).
- Decomposed upfront bulk-loading (`proj-load-all`) into per-reference `*-setup` chunks via `load_proj()` helper + `tar_read_raw()`.
- Violin plots: `ref_plot_violins()` now returns patchwork directly; bold compartment label separators inserted between compartment blocks; x-axis text size 8→11; cluster violin x-label shows full column name.
- Added `# Cell density & mixing` section to `compositional-analysis.qmd`: point density UMAP (ggpointdensity/magma) + `MixingMetric` UMAP (turbo), per compartment.
- Added `# Mixing metric & group-correlated genes` section to `cluster-annotation.qmd`: `MixingMetric` + `GroupCorrelation(var="mix")` per compartment; results as searchable DT table sorted by `|mix_r|`.

### 2026-04-12 — Cluster annotation pipeline + Azimuth label transfer + DE/ORA targets moved into pipeline

- Added `notebooks/cluster-annotation.qmd`: marker-based (dot plots + UMAP across 8 resolutions) and reference-based (Azimuth label transfer from 9 curated references × 3 compartments, ARI line plots, mapping score diagnostics, resolution heatmaps) annotation workflow.
- Added 27 Azimuth projection `tar_map` blocks in `_targets.R` (9 references × stroma/hsc/mnc): TMS Marrow/Spleen/Lung/Kidney, Mouse BM Atlas, MCA PBMC/BM/BM-cKit/BM-Mesenchyme; outputs as per-compartment CSV.
- Added `annotation_recluster` targets: `recluster_resolutions_df()` runs `FindClusters` at resolutions 0.1–0.9 and returns a plain data.frame CSV.
- Added lightweight annotation extract targets: `anno_meta`, `anno_expr`, `anno_raw_meta`, `ref_umap_df_*` — replace Seurat object loads in all notebooks.
- Moved Milo, ORA, module scores, and biased-expression targets into `_targets.R` pipeline: `milo_stroma/hsc/mnc`, `allvsall_genesets_*`, `ora_results_*`, `module_scores_*`, `biased_feature_expr_*`.
- Removed `output_dir` stamped-run-key target (mutated yaml at runtime causing false invalidation); replaced with static `path_params` from `paths:` yaml block.
- Removed all per-subset file-output targets (umap_*, elbow_*, bubble_*, marker_feature_*) and `report` target from `_targets.R`.
- Added `dim_grid_sct_pre/runs` targets: pre-computed dim-reduction grid (elbow/PCA/UMAP data.frames per dims×k sweep) to support dim-reduction notebook without inline recomputation.
- Added ~1800 lines of new helpers to `scrna-seq_utilities.R`: reference loaders, Azimuth wrappers, annotation plotting, ARI diagnostics, data-frame UMAP helpers, Milo/ORA/module-score helpers.
- `project_config.yaml`: added `output_table_dir`, `output_figure_dir`, `azimuth_dir`; removed `run_key/run_label`; added `dim_grid_dims/k/res_scan`.

### 2026-04-16 — Machine migration: WSL2 → macOS M5 Pro (Lima + Apptainer)

- Migrated project from WSL2 x86_64 to macOS Apple M5 Pro (arm64) using Lima VM + Apptainer.
- Fixed `bioinfo.def`: Quarto download URL now uses `dpkg --print-architecture` instead of hardcoded `amd64`.
- Added `container/python_setup.py` placeholder (required by `build.sh` pre-flight check).
- Fixed all run scripts (`run_rstudio.sh`, `run_shell.sh`, `run_targets.sh`): override `$HOME` to strip Lima's `.guest` username suffix; updated default data path to `data_hub/mpn_10x_camelia`.
- `run_rstudio.sh`: `RSTUDIO_RUN` moved to `/tmp` (virtiofs cannot create Unix sockets); added `--auth-minimum-user-id=500` for macOS UIDs < 1000.
- Lima VM created from explicit config (`~/.lima/bioinfo.yaml`) with `writable: true`, 14 CPUs / 40 GiB RAM / 200 GiB disk; Apptainer installed via `ppa:apptainer/ppa` (v1.4.5).

### 2026-04-20 — External geneset pipeline + HemAtlas processing script + container updates

- Added `notebooks/geneset-exploration.qmd`: 14-section HSC geneset browser with plotly module score UMAPs (Nestorowa, Isobe, Psaila, Pan, Gavish meta/cancer) and FeaturePlots (Kucinski, HemAtlas, Ferchen, Gavish mm); `umap_df` centralised in setup chunk; `geneset_dt()` + `score_scatter()` helpers.
- Added 14 geneset targets to `_targets.R`: 9 mouse-native (`_mm_`) via `readr::read_csv() |> with(split(...))` and 5 human-sourced (`_hs_`) with `inner_join(mm_hs_ortho_long, ...)` ortholog mapping; all named with species tag.
- Added 6 module score targets in pipeline (`module_scores_nestorowa_17_bl_mm_hsc`, `module_scores_isobe_23_cg_mm_hsc`, `module_scores_psaila_20_hs_hsc`, `module_scores_pan_25_hs_hsc`, `module_scores_gavish_23_hs_meta_program_hsc`, `module_scores_gavish_23_hs_cancer_mp_hsc`).
- Rewrote `extract_module_scores`: uses `dplyr::rename(all_of(rename_map))` to recover geneset names from AddModuleScore's `Cluster1..N` output; returns `cell_id` + scores only.
- Added `data_hub/hematlas/process_hematlas.R`: extracts RDS from tar.gz → rebuilds SCTransform/PCA/UMAP from `SCT@counts` → `DietSeurat` → `AzimuthReference`; Azimuth block runs raw in `for` loop (global env requirement).
- Fixed `data_hub/tusi_2018_nature/convert_tusi_to_10x.py`: barcode source changed from duplicated `meta["barcode"]` column to `meta.index` (unique `cell_id`).
- Added `scCustomize` (CRAN) + `dittoSeq`, `DropletUtils`, `Nebulosa` (Bioconductor) to `container/renv_setup.R`.

### 2026-04-21 — CellRank + pySCENIC Python workflows + notebook HTML fixes

- `notebooks/velocity_run.py`: CellRank 2.x workflow — PseudotimeKernel (CytoTRACE2 inverted: `1 - score`) + ConnectivityKernel → CFLARE → fate probabilities (GMRES+ILU) → per-state lineage driver PNGs → gene trends.
- `notebooks/cellrank-hsc.qmd`: presentation notebook; Methods covers Markov math, CFLARE vs GPCCA, inversion rationale, troubleshooting guide.
- `notebooks/scenic_run.py` (new): pySCENIC GRN → cisTarget → AUCell for HSC; 5k stratified-subsample for GRN, all cells for AUCell; resume logic per step; outputs adjacencies, regulons, auc_mtx, RSS, h5ad.
- `notebooks/scenic-hsc.qmd` (new): SCENIC results notebook — regulon DT table, size distribution, top-9 UMAPs, cluster/condition heatmaps (pheatmap), RSS table, Methods.
- HTML fix: added `embed-resources: true, theme: cosmo` to all 13+ notebook frontmatters; `_quarto.yml` gets `embed-resources: true`; fixes broken CSS when notebooks rendered outside book context.
- `scripts/render_notebooks.sh`: Lima `.guest` suffix fix; `~/data_hub` data path; in-place render + mv workaround.
- `scripts/run_jupyter.sh` (new): JupyterLab via Apptainer, port 8888, `PYTHONPATH=/project/cache/python_pkgs`.
- `scripts/download_scenic_db.sh` (new): correct aertslab URL (`mc_v10_clust/gene_based/`); dual screen+log; 0-byte guard.
- `src/scrna-seq/scrna-seq_utilities.R`: `anno_make_type_pal` empty-type guard; `extract_azimuth_umap_df` column-name alias resolution; `plot_ref_umap` missing-column guard.

### 2026-04-21 — pySCENIC cisTarget debugging: pandas 3.0 + numba compatibility

- Patched `cache/python_pkgs/pyscenic/utils.py`: `.values` → `.values.copy()` (numba readonly fix); `groupby.apply(nlargest)` → `sort_values + groupby.head` (pandas 3.0 `include_groups` breaking change, `KeyError: 'target'`).
- `notebooks/scenic_run.py`: dynamic pyscenic binary detection (bin path → `shutil.which` → Python entry-point fallback); added `-B` to entry-point subprocess to prevent stale `.pyc` masking source edits.
- `scripts/run_scenic.sh`: added `-B` to outer `python3` invocation.
- pySCENIC pipeline completed: `regulons.csv` (11.7 MB) + `auc_mtx.csv` (145.1 MB) in `output/scenic/`.
- Bash tool runs on macOS host (no apptainer); fixed by using `limactl shell bioinfo` to reach Lima VM.

### 2026-05-11 — Human reference processing (Suo/Jardine) + Moran's I targets + container fixes

- Added `data_hub/suo_2022/convert_h5ad_to_newformat.py` and `data_hub/jardine_2021/convert_h5ad_to_newformat.py`: pure h5py in-place fix for `nullable-string-array` HDF5 GROUP encoding → flat `string-array` DATASET; makes old-format h5ad readable by `anndataR`.
- Added `process_suo_refs.R` and `process_jardine_refs.R`: build Seurat reference objects (LogNormalize → PCA → UMAP with `return.model=TRUE`) for human label transfer; Harmony enabled only for Suo HSC_PROGENITORS; all layers converted dense→sparse before saving; `scale.data` dropped from RDS.
- `_targets.R`: added `#jardine` `tar_map` section (18 targets mirroring `#suo`); added `moransi_hsc` and `moransi_pca_hsc` `tar_map` (3 conditions × 2 coordinate systems = 12 targets); both include BH-adjusted p-values.
- `scrna-seq_utilities.R`: added `moran_sparse()` (sparse-matrix Moran's I, no n² allocation) and `nn_to_sym_dist()` (asymmetric KNN → symmetric sparse distance matrix).
- `container/bioinfo.def`: added `libudunits2-dev libgdal-dev libgeos-dev libproj-dev` (required by `sf`/`spdep`); removed `celloracle` (pins `pandas<=1.5.3` — incompatible with Python 3.12).
- `container/renv_setup.R`: added `Rfast2` (required by `Seurat::RunMoransI` for fast computation).

### 2026-04-25 — SuperCell refinement + Milo per-cell scores + scatter-pie UMAP

- `compute_supercells()` in utilities: added `split_by` and `fast.merge` parameters; `split_by = NULL` disables per-sample splitting (for trajectory use); `fast.merge = TRUE` merges small walktrap communities to reduce right-skew in metacell size distribution.
- `supercells_hsc` target updated: `fast.merge = TRUE`, `split_by = "sample.name"` now explicit.
- Added `milo_hsc_cell_scores` target: aggregates Milo neighbourhood logFC to per-cell mean logFC using the `cells × nhoods` membership matrix from `milo_hsc_obj`; returns named list of three data frames (one per contrast).
- Added `plot_medoid_pie()` to `compositional-analysis.qmd` helpers: `scatterpie`-based UMAP where each group is a pie placed at its PCA-space medoid cell; pie area ∝ sqrt(n_cells).
- Key insight documented: SuperCell `gamma` is applied per `cell.split.condition` group, not globally — effective graining = n_cells_per_group / gamma; many small samples cause right-skewed metacell sizes and lower-than-expected total metacell count.

### 2026-04-25 — Container build fixes + GHCR push script

- `container/bioinfo.def`: added `/root/.R/Makevars` with `CXX11STD = -std=gnu++14` and `CXX14STD = -std=gnu++14` before `Rscript /opt/renv_setup.R`; fixes presto/RcppArmadillo C++14 compile failure.
- `container/renv_setup.R`: standalone `miloR` install + retry; `requireNamespace` checks with early `stop()` for `presto`; removed `Azimuth`, `Signac`, and `cellxgene.census` install blocks.
- Azimuth/Signac dropped: fragile on R 4.5 + GitHub rate-limit failures; label-transfer results already cached as CSVs in `_targets.R`.
- cellxgene.census/tiledb dropped: tiledb has no Linux aarch64 support — fails on Lima VM (aarch64 build, amd64 deploy).
- `scripts/push_container.sh` (new): pushes `.sif` to GHCR as OCI artifact; reads `GHCR_USER` + `GHCR_TOKEN`; supports `--tag` and `--sif` flags.
- Pull on cluster: `apptainer pull oras://ghcr.io/mohamedrefaat92/bioinfo:latest`

### 2026-04-25 — run_velocyto.sh updated for CellRanger Multi layout

- Rewrote `scripts/run_velocyto.sh` for CellRanger Multi output: BAM is `per_sample_outs/*/count/sample_alignments.bam`; switched `velocyto run10x` → `velocyto run` with explicit `--bcfile` / `--sampleid`.
- Bind mounts: `BIND_DATA_HUB` replaced by `BIND_DATA` (→ `/data`) + `BIND_GENOMES` (→ `/genomes`); required-env guards added.
- GTF decompressed from `.gz` to `/tmp/genes.gtf` at runtime (velocyto does not accept gzipped GTF).
- Lima-specific `HOME` override removed; script now targets cluster environment directly.

### 2026-05-28 — Mpig6b expression notebook

- Added `notebooks/mpig6b-expression.qmd`: explores *Mpig6b* (platelet CLEC-2/ITAM) expression in raw (`RNA` assay, pre-gene-filter) and processed (`SCT`) Seurat objects for HSC, MNC, and Stroma.
- Raw counts sourced from `mpn_raw_*`; UMAP coords from `anno_raw_meta_*` (pre-QC); post-QC views use `mpn_processed_data_*` + `anno_meta_*`.
- Gene-filter check reports whether Mpig6b exceeds `min_cells_per_gene = 20` threshold per compartment; SCT sections suppressed via `eval=` when gene is absent.
- Fixed `pass_cell` not found error: `anno_meta_*` lacks this column (all processed cells already passed QC); fixed with conditional filter in `violin_condition()` / `violin_cluster()`.
- Load → extract → `rm` / `gc()` pattern maintained throughout to respect 16 GB memory limit.

### 2026-10-27 — hdWGCNA pipeline integration, per-condition transition matrices, trajectory refinement, genomic annotation tools

- Moved hdWGCNA computation into pipeline: new targets `hdwgcna_hsc_modules`, `hdwgcna_hsc_hvg`, `hdwgcna_hsc_soft_power`, `hdwgcna_hsc_metacells`, `hdwgcna_hsc_network`; TOM written to `output/hdwgcna/`; `notebooks/hdwgcna-hsc.qmd` now loads from cache (all inline SetupForWGCNA/ConstructNetwork/ModuleEigengenes calls removed).
- Added merged multi-compartment targets: `mpn_all`, `mpn_wt_all`, `mpn_mch_neg_all`, `mpn_mch_pos_all` (HSC + MNC + Stroma merged by condition; Stroma Het treated as mCh+).
- Added external data file targets: `yoshida_2020_ocr_csv` (ImmGen ATAC OCRs), `hpa_secretome_tsv`, `naba_matrisome_mm_csv`; added `sex_degs` gene vector, `hem_markers_claude_df`, `cytotrace_hsc_wt`, per-compartment pairwise comparison targets.
- `notebooks/transition-matrices.qmd`: per-condition (WT / mCh− / mCh+) cluster transition matrices for all 6 matrix types; CellRank Markov cluster transition matrix loaded from `cellrank/cytotrace/`; path constant updated from `velocity_dynamical` → `cytotrace`.
- `notebooks/trajectory-slingshot.qmd`: refactored to use `mpn_processed_hsc_fullcc`; dedicated config block (`CLUS_COL`, `EXCLUDE_CLUSTER`, `START_CLUSTER`, `END_CLUSTER`, `N_PCA_DIMS`); cells filtered before PCA passed to `slingshot()`.
- `notebooks/cellrank_run.py`: added `dpt` as third `INPUT_TYPE`; configurable `CLUSTER_RES` env var; nullable string dtype fix (`_fix_nullable_strings()`).
- `src/scrna-seq/scrna-seq_utilities.R`: added `find_delegate_pairwise()`, genomic range helpers (EnsDb → GRanges, ImmGen OCR parsing), pseudotime binning (`bin_expr_along_pseudotime`, `plot_lineage_heatmap`), `compute_magic()`, `plot_nhood_graph_groups()`.
- All notebooks: `embed-resources: true` → `embed-resources: false` to reduce HTML file size.

### 2026-11-02 — Moved into the hyprid-proj template (this repository)

- Recast from the container-based project into the template layout: `_targets.R` + `src/R/scrna-seq_utilities.R`; config in `metadata/project_config.yaml`; inputs as `data/` links; outputs in `output/` (`tables/`, `objects/` for h5ad and the Python tools' results).
- Environment rebuilt in three layers instead of the bioinfo container: C libraries in `environment.yml`, R packages via renv (`scripts/R/install-r-packages.R`: CRAN, Bioconductor 3.22, 13 from GitHub), Python via uv (`pyproject.toml`); pySCENIC and velocyto in pinned containers.
- Cell Ranger multi added to the project (`scripts/slurm/cellranger/`): 8 tubes from the FASTQ tar archives with Cell Ranger 10.1.0 in a private container (`scripts/containers/cellranger.def`, from the downloaded release), same reference and settings as the first processing (9.0.0); results in the data hub, linked as `data/cellranger_Rx_<n>`. Cell Ranger 10 writes each sample's matrix and BAM directly in `per_sample_outs/<id>/` (9.0 had a `count/` subfolder).
- 64 notebooks became pipeline reports (`notebooks/<name>/`, provenance + session info); 9 that read removed targets went to `notebooks/archive/`; scratch, daily and Jupyter notebooks to `notebooks/exploration/`.
- Fixes: `proj_path_hematlas_mm_agm` wrote over the MCA BM CSV; two notebooks had scrambled `include-after-body` YAML; `scatter_df_mnc`/`_stroma` (and their genesets / ORA / biased-expression targets) read targets that no longer exist and are commented out; scenic_run.py split into a prepare step (uv) and an infer step (container), GRN through the pyscenic CLI.
