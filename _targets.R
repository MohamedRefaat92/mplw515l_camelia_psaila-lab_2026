# The pipeline: MPN 10x scRNA-seq analysis (Stroma, HSC and MNC compartments).
# Run with tools/run-pipeline.sh, on the cluster with sbatch scripts/slurm/pipeline/run-pipeline.sh
# (or targets::tar_make() in R). The store lives in output/_targets (set in _targets.yaml).
# Adapted from Shoeb2025_CG-AMKL (https://github.com/cancerbits/Shoeb2025_CG-AMKL)
library(targets)
library(tarchetypes)

tar_source(c("src/R", "tools/R/provenance.R"))   # analysis functions + bookkeeping
# (not all of tools/R: some files run when sourced, e.g. register-r-kernel.R)

ggplot2::theme_set(cowplot::theme_cowplot())

# phateR and Rmagic run Python packages (phate, magic-impute) through reticulate: the project's
# .venv has them. Set here, it also reaches the reports, which targets renders in child processes.
Sys.setenv(RETICULATE_PYTHON = file.path(getwd(), ".venv", "bin", "python"))

#find available cores
n_cpus <- parallelly::availableCores()

# Analysis parameters, read once when the pipeline is defined. Targets take their values with
# !!config$<section>$<key>: the value is written into the target's command, so changing one
# parameter re-runs exactly the targets that use it (and what depends on them), nothing else.
config <- yaml::read_yaml("metadata/project_config.yaml")

# ---------------------------------------------------------------
# Global options
# ---------------------------------------------------------------
tar_option_set(
  # A failing target doesn't stop the others: several reports read outputs of the cluster jobs
  # in scripts/slurm/ (SCENIC, CellRank, PAGA), which may not exist yet. See targets::tar_errored().
  error = "continue",
  packages = c(
    "Seurat",
    "scDblFinder",
    "deMULTIplex2",
    "magrittr",
    "tidyverse",
    "patchwork",
    "cowplot",
    "ggrepel",
    "yaml",
    "quarto",
    "future",
    "rlang",
    "reshape2",
    "miloR",
    "SingleCellExperiment",
    "fgsea",
    "msigdbr",
    #"SCENT",
    "org.Mm.eg.db",
    "org.Hs.eg.db",
    "AnnotationDbi"
  ),
  format             = "qs",
  memory             = "transient",
  garbage_collection = 1,
  resources = tar_resources(qs = tar_resources_qs(nthreads = n_cpus))
)


# Reports: one folder per report under notebooks/, rendered to output/reports/latest/.
# extra_files: also re-render when the report style, the provenance helper, the reviewed data catalogue
# or renv.lock changes, or a file the report reads outside targets (passed in ...).
report_extra_files <- c("_quarto.yml", "_environment", "tools/R/provenance.R", "metadata/catalogue_data.tsv", "renv.lock")
# Like tar_quarto_raw(), without its `quarto inspect` call (~12 s per report, every time the pipeline
# is defined). The output path follows from _quarto.yml (output-dir mirrors notebooks/<name>/) and the
# included files are the shared ones in notebooks/; dependencies come from tar_knitr_deps() (R code only).
report_output_dir <- "output/reports/latest"
report_includes   <- c("notebooks/_provenance.qmd", "notebooks/_session-info.qmd")
# Usage, like tar_target(): quarto_report(report_qc_exploration, "qc-exploration", <extra files>...)
# — the target name (unquoted), then the report folder under notebooks/.
quarto_report <- function(target, name, ...) {
  target  <- deparse(substitute(target))
  qmd     <- file.path("notebooks", name, paste0(name, ".qmd"))
  sources <- c(qmd, report_includes)
  deps    <- sort(unique(unlist(lapply(sources, tarchetypes::tar_knitr_deps))))
  command <- substitute(
    tarchetypes:::tar_quarto_run(
      args = list(input = QMD, execute = TRUE, execute_dir = getwd(), execute_daemon = 0,
                  execute_daemon_restart = FALSE, execute_debug = FALSE, cache_refresh = FALSE,
                  debug = FALSE, quiet = TRUE, as_job = FALSE),
      deps = DEPS, sources = SOURCES, output = OUTPUT, input = INPUT
    ),
    list(
      QMD     = qmd,
      DEPS    = as.call(c(as.symbol("list"), lapply(deps, as.symbol))),
      SOURCES = sources,
      OUTPUT  = file.path(report_output_dir, "notebooks", name, paste0(name, ".html")),
      INPUT   = c(report_extra_files, ...)
    )
  )
  tar_target_raw(
    target,
    command    = as.expression(command),
    packages   = "quarto",
    format     = "file",
    repository = "local"
  )
}
utilities_r <- "src/R/scrna-seq_utilities.R"    # sourced by many reports
config_yaml <- "metadata/project_config.yaml"   # read by some

# ---------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------
list(

  # ── Config, for the reports ─────────────────────────────────────
  # The pipeline's own targets inject parameters with !!config$... (top of this file). These
  # targets hand sections of the config to reports (tar_load(qc_params) etc.), so a report
  # re-renders when a section it reads changes. proc_params and grid_params keep their old flat
  # names (ref_assay, marker_resolution, dim_grid_dims, ...), which the reports use.
  # (!!! splices a section's entries into list(...): a command must be a call, not a list.)
  tar_target(path_params,       list(!!!config$paths)),
  tar_target(qc_params,         list(!!!config$qc)),
  tar_target(annotation_config, list(!!!config$annotation)),
  tar_target(proc_params, list(
    !!!config$reduction,
    cluster_resolutions        = !!config$clustering$resolutions,
    marker_resolution          = !!config$markers$resolution,
    marker_min_pct             = !!config$markers$min_pct,
    marker_logfc_threshold     = !!config$markers$logfc_threshold,
    marker_max_cells_per_ident = !!config$markers$max_cells_per_ident,
    transgene_features         = !!config$features$transgene_features
  )),
  tar_target(grid_params, list(
    dim_grid_dims     = !!config$dim_grid$dims,
    dim_grid_k        = !!config$dim_grid$k,
    dim_grid_res_scan = !!config$dim_grid$res_scan
  )),
  # Cell-type colours (named vector)
  tar_target(colors_vec, c(!!!unlist(config$plotting$colours))),

  # ── Raw data ─────────────────────────────────────────────────────
  # The sample sheet: one row per tube (Tube, sample_name_short, sample_name_long)
  tar_target(
    tube_sample_map_file,
    "metadata/tube_sample_map.csv",
    format = "file"
  ),
  tar_target(
    tube_sample_map,
    readr::read_csv(tube_sample_map_file, show_col_types = FALSE)
  ),
  tar_target(
    sample_id_vec,
    sort(tube_sample_map$Tube)
  ),
  # Each tube's Cell Ranger multi matrix (scripts/slurm/cellranger/), linked into data/ by
  # scripts/sh/link-cellranger.sh
  tar_target(
    cellranger_matrix_dirs,
    file.path("data", paste0("cellranger_", sample_id_vec), "sample_filtered_feature_bc_matrix"),
    format = "file"
  ),
  tar_target(
    sample_name_long,
    with(tube_sample_map, setNames(sample_name_long, Tube))
  ),
  tar_target(
    sample_name_short,
    with(tube_sample_map, setNames(sample_name_short, Tube))
  ),

  # Load each CellRanger sample into a named list of Seurat objects.
  tar_target(
    mpn_per_sample,
    {
      srt_lst <- Map(read_10x_sample, cellranger_matrix_dirs, sample_id_vec,
                     MoreArgs = list(sample_name_short = sample_name_short,
                                     sample_name_long  = sample_name_long,
                                     transgene_features = unlist(!!config$features$transgene_features)))
      names(srt_lst) <- sample_id_vec
      srt_lst
    }
  ),

  # Merge all per-sample objects and run QC scoring on the combined object.
  tar_target(
    mpn_raw_all,
    merge_seurat_list(mpn_per_sample)
  ),

  # sex genes that are consistently upregulated between Het and WT mice due to sex difference
  tar_target(sex_degs,  c("Xist","Tsix", "Eif2s3y", "Xlr", "Xlr3b", "Ddx3y", "Uty", "Kdm5d")),
  # mt/rp/hb patterns come from qc_params; excluded_genes adds reporter genes.
  tar_target(mt_genes,       grep(!!config$qc$mt_pattern, rownames(mpn_raw_all), value = TRUE)),
  tar_target(rp_genes,       grep(!!config$qc$rp_pattern, rownames(mpn_raw_all), value = TRUE)),
  tar_target(hb_genes,       grep(!!config$qc$hb_pattern, rownames(mpn_raw_all), value = TRUE)),
  tar_target(hsp_genes,       grep(!!config$qc$hsp_pattern, rownames(mpn_raw_all), value = TRUE)),
  tar_target(excluded_genes, c(mt_genes, rp_genes, hsp_genes, "Malat1")),   # transgenes live in their own assay

  # Per-compartment QC flagging — subset mpn_raw_all, compute QC metrics and
  # pass/fail flags, and run pre-filter dimensionality reduction.
  # Produces: mpn_raw_stroma, mpn_raw_hsc, mpn_raw_mnc
  #           mpn_raw_processed_stroma, mpn_raw_processed_hsc, mpn_raw_processed_mnc
  tarchetypes::tar_map(
    values = tibble::tibble(
      name             = c("stroma",  "hsc",  "mnc"),
      compartment_name = c("Stroma",  "HSC",  "MNC")
    ),
    names = name,

    tar_target(demux_res,
               demux_per_rx(subset(mpn_raw_all, compartment ==compartment_name))
               ),
    tar_target(mpn_raw,
              flag_qc_metadata(subset(mpn_raw_all, compartment ==compartment_name),
                               min_cells_per_gene = !!config$qc$min_cells_per_gene,
                               min_counts_per_cell_per_gene = !!config$qc$min_counts_per_cell_per_gene,
                               mito_percent_max   = !!config$qc$mito_percent_max,
                               hb_percent_max     = !!config$qc$hb_percent_max,
                               mt_genes           = mt_genes,
                               rp_genes           = rp_genes,
                               hb_genes           = hb_genes,
                               stress_genes       = unlist(!!config$qc$stress_genes),
                               compartment_qc     = !!config$qc$compartment_qc,
                               demux_out = demux_res)
          ),

    tar_target(
      mpn_raw_processed,
      process_all(mpn_raw,
                  ref_assay     = !!config$reduction$ref_assay,
                  dims_i        = seq(!!config$reduction$n_pca_dims),
                  k             = !!config$reduction$k_neighbors,
                  res           = unlist(!!config$clustering$resolutions),
                  n_features    = !!config$reduction$n_features,
                  k.nn          = !!config$reduction$n_umap_neighbors,
                  exclude_genes = excluded_genes)
    )
  ),

  # ── Processing — stroma / hsc / mnc ─────────────────────────────
  # Each compartment: apply cell + gene QC filters → full processing
  # (SCTransform, PCA, SNN, clustering, UMAP) → slim DietSeurat →
  # find cluster markers.
  # Produces: mpn_clean_<x>, mpn_processed_<x>, mpn_mpn_processed_data_<x>,
  #           subset_marker_features_<x>, cluster_markers_<x>
  tarchetypes::tar_map(
    values = tibble::tibble(
      # Not `subset`: tar_map() replaces every symbol named like a values column, so subset() calls broke
      compartment = c("stroma",      "hsc",      "mnc"),
      raw_subset  = rlang::syms(c("mpn_raw_stroma", "mpn_raw_hsc", "mpn_raw_mnc"))
    ),
    names = compartment,

    tar_target(
      pass_gene,
      {
        # Genes detected (>= min_counts_per_cell_per_gene counts) in >= min_cells_per_gene QC-passing cells
        cells_pass  <- colnames(raw_subset)[raw_subset@meta.data$pass_cell]
        counts      <- Seurat::GetAssayData(raw_subset, assay = "RNA", layer = "counts")[, cells_pass]
        detected    <- counts >= !!config$qc$min_counts_per_cell_per_gene   # gene × cell, TRUE/FALSE
        n_cells_det <- Matrix::rowSums(detected)                          # per gene: cells detected in
        names(n_cells_det)[n_cells_det >= !!config$qc$min_cells_per_gene]
      }
    ),

    tar_target(
      mpn_clean,
      # QC-passing cells; RNA genes in pass_gene. Features are subset across all assays (an assay left
      # with none is dropped), so the Transgene and HTO features are listed to keep those assays whole.
      subset(raw_subset,
             cells    = colnames(raw_subset)[raw_subset$pass_cell],
             features = c(pass_gene,
                          rownames(raw_subset[["Transgene"]]),
                          rownames(raw_subset[["HTO"]])))
    ),

    tar_target(
      mpn_processed,
      process_all(mpn_clean,
                  ref_assay     = !!config$reduction$ref_assay,
                  dims_i        = seq(!!config$reduction$n_pca_dims),
                  k             = !!config$reduction$k_neighbors,
                  res           = unlist(!!config$clustering$resolutions),
                  n_features    = !!config$reduction$n_features,
                  k.nn          = !!config$reduction$n_umap_neighbors,
                  exclude_genes = excluded_genes)
    ),

    # Slim object: keep only the normalised data layer + reductions + graphs.
    # Downstream targets use mpn_processed_data to avoid retaining the raw counts.
    tar_target(
      mpn_processed_data,
      DietSeurat(mpn_processed,
                 layers    = "data",
                 assays    = !!config$reduction$ref_assay,
                 dimreducs = Seurat::Reductions(mpn_processed),
                 graphs    = SeuratObject::Graphs(mpn_processed))
    ),

    tar_target(
      subset_marker_features,
      rownames(mpn_processed_data)
    ),
    # Cluster markers using DElegate
    tar_target(delegate_res_cluster_rest_markers,
      {
        options(future.globals.maxSize = 8 * 1024^3)
        srt2       <- mpn_processed
        rm(mpn_processed); gc()
        res_vec <- !!config$clustering$resolutions
        res_vec <- paste0("sNN_pca.SCT_30_res.", res_vec)
        names(res_vec) <- res_vec

        lapply(res_vec, function(marker_res){
              cat(marker_res,"\n")
          Seurat::Idents(srt2) <- marker_res
          find_delegate_markers_rest_clusters(srt2, excluded_genes, group_column = marker_res)
        })
      }
    ),
    # PCA embedding + metadata data.frame for PCA distribution notebook.
    tar_target(
      pca_df,
      {
        srt        <- mpn_processed_data
        marker_res <- !!config$markers$resolution
        clust_col  <- grep(paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", marker_res), "$"),
                           colnames(srt@meta.data), value = TRUE)[1]
        Seurat::Idents(srt) <- clust_col
        extract_pca_df(srt,
                       pca_name = paste0("pca.", !!config$reduction$ref_assay),
                       n_pcs    = !!config$reduction$n_pca_dims)
      }
    )
  ),

  # ── Cell-cycle-corrected HSC processing ───────────────────────────────────
  # "full"  strategy: regress out S.Score + G2M.Score (removes all CC signal).
  # "diff"  strategy: regress out CC.Difference = S.Score - G2M.Score
  #         (removes G1-vs-cycling bias; preserves proliferation signal).
  # Both targets use mpn_clean_hsc so they are independent of mpn_processed_hsc.
  tar_target(
    mpn_processed_hsc_fullcc,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      process_all_cc(mpn_clean_hsc,
                     ref_assay   = !!config$reduction$ref_assay,
                     dims_i      = seq(!!config$reduction$n_pca_dims),
                     k           = !!config$reduction$k_neighbors,
                     res         = unlist(!!config$clustering$resolutions),
                     n_features  = !!config$reduction$n_features,
                     k.nn        = !!config$reduction$n_umap_neighbors,
                     exclude_genes = excluded_genes,
                     cc_strategy = "full")
    }
  ),

  tar_target(
    mpn_processed_hsc_diffcc,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      process_all_cc(mpn_clean_hsc,
                     ref_assay   = !!config$reduction$ref_assay,
                     dims_i      = seq(!!config$reduction$n_pca_dims),
                     k           = !!config$reduction$k_neighbors,
                     res         = unlist(!!config$clustering$resolutions),
                     n_features  = !!config$reduction$n_features,
                     k.nn        = !!config$reduction$n_umap_neighbors,
                     exclude_genes = excluded_genes,
                     cc_strategy = "diff")
    }
  ),

  tar_target(
    mpn_processed_hsc_wt_fullcc,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      srt_wt <- subset(mpn_clean_hsc, sample.name == "HSC_WT_CD41+")
      process_all_cc(srt_wt,
                     ref_assay   = !!config$reduction$ref_assay,
                     dims_i      = seq(!!config$reduction$n_pca_dims),
                     k           = !!config$reduction$k_neighbors,
                     res         = unlist(!!config$clustering$resolutions),
                     n_features  = !!config$reduction$n_features,
                     k.nn        = !!config$reduction$n_umap_neighbors,
                     exclude_genes = excluded_genes,
                     cc_strategy = "full")
    }
  ),
  tar_target(
    mpn_processed_hsc_mchneg_fullcc,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      srt_mchneg <- subset(mpn_clean_hsc, sample.name == "HSC_Het_CD41+mCh-")
      process_all_cc(srt_mchneg,
                     ref_assay   = !!config$reduction$ref_assay,
                     dims_i      = seq(!!config$reduction$n_pca_dims),
                     k           = !!config$reduction$k_neighbors,
                     res         = unlist(!!config$clustering$resolutions),
                     n_features  = !!config$reduction$n_features,
                     k.nn        = !!config$reduction$n_umap_neighbors,
                     exclude_genes = excluded_genes,
                     cc_strategy = "full")
    }
  ),
  tar_target(
    mpn_processed_hsc_mchpos_fullcc,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      srt_mchpos <- subset(mpn_clean_hsc, sample.name == "HSC_Het_CD41+mCh+")
      process_all_cc(srt_mchpos,
                     ref_assay   = !!config$reduction$ref_assay,
                     dims_i      = seq(!!config$reduction$n_pca_dims),
                     k           = !!config$reduction$k_neighbors,
                     res         = unlist(!!config$clustering$resolutions),
                     n_features  = !!config$reduction$n_features,
                     k.nn        = !!config$reduction$n_umap_neighbors,
                     exclude_genes = excluded_genes,
                     cc_strategy = "full")
    }
  ),

  #calculate Adjusted Rand Index for the cell cycle phases for the different strategies
  tar_target(mpn_hsc_ari_df,
             ari_df <- dplyr::bind_rows(
               compute_ari(mpn_processed_hsc@meta.data,        "Phase",
                           str_subset(colnames(mpn_processed_hsc@meta.data),"sNN_pca."))        |> dplyr::mutate(strategy = "none"),
               compute_ari(mpn_processed_hsc_fullcc@meta.data, "Phase",
                           str_subset(colnames(mpn_processed_hsc_fullcc@meta.data),"sNN_pca.")) |> dplyr::mutate(strategy = "full"),
               compute_ari(mpn_processed_hsc_diffcc@meta.data, "Phase",
                           str_subset(colnames(mpn_processed_hsc_diffcc@meta.data),"sNN_pca.")) |> dplyr::mutate(strategy = "diff")
             )
    ),
  #calculate gene correlation with S and G2M score
  tar_target(cc_cor,
             {
             full_data <- GetAssayData(mpn_processed_hsc_fullcc, assay = "SCT", layer = "data")
             cc_scores  <- mpn_processed_hsc_fullcc@meta.data[, c("S.Score", "G2M.Score")]

             # Compute gene × {S.Score, G2M.Score} correlations efficiently via crossprod
             common_cells <- intersect(colnames(full_data), rownames(cc_scores))
             expr_t        <- t(full_data[, common_cells])          # cells × genes (dense slice)
             cc_cor        <- cor(as.matrix(expr_t), as.matrix(cc_scores[common_cells, ]))
             }
             ),
  tar_target(cc_cor_top_genes,
             apply(cc_cor,1, function(x){min(x) >= 0.4})
             ),

  # ── hdWGCNA HSC — metacell aggregation → TOM network → modules ────────────
  # §2–§6 of notebooks/hdwgcna-hsc.qmd, cached to avoid re-running TOM.
  # GRN inference (MotifScan / ConstructTFNetwork) stays in the notebook.

  tar_target(
    hdwgcna_hsc_hvg,
    {
      excluded <- c(hsp_genes,
                    names(cc_cor_top_genes)[cc_cor_top_genes])
      setdiff(VariableFeatures(mpn_processed_hsc_fullcc), excluded)
    }
  ),

  tar_target(
    hdwgcna_hsc_setup,
    hdWGCNA::SetupForWGCNA(
      mpn_processed_hsc_fullcc,
      wgcna_name  = "hsc_modules",
      gene_select = "variable",
      features    = hdwgcna_hsc_hvg
    )
  ),

  tar_target(
    hdwgcna_hsc_metacells,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      srt <- hdWGCNA::MetacellsByGroups(
        seurat_obj       = hdwgcna_hsc_setup,
        group.by         = c("sample.name"),
        k                = 25,
        max_shared       = 15,
        target_metacells = 1000,
        ident.group      = "sample.name",
        wgcna_name       = "hsc_modules",
        reduction        = "pca.SCT"
      )
      hdWGCNA::NormalizeMetacells(srt, wgcna_name = "hsc_modules")
    }
  ),

  # SetDatExpr + TestSoftPowers + ConstructNetwork bundled:
  # ConstructNetwork writes the TOM file; soft power must be derived first.
  tar_target(
    hdwgcna_hsc_network,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      tom_dir <- file.path(!!config$paths$objects_dir, "hdwgcna")
      dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
      srt <- hdWGCNA::SetDatExpr(
        hdwgcna_hsc_metacells,
        group_name = unique(hdwgcna_hsc_metacells$sample.name),
        group.by   = "sample.name",
        assay      = "SCT",
        layer      = "data",
        wgcna_name = "hsc_modules"
      )
      srt      <- hdWGCNA::TestSoftPowers(srt, wgcna_name = "hsc_modules")
      sp_tbl   <- srt@misc[["hsc_modules"]]$wgcna_powerTable
      soft_pow <- sp_tbl$Power[min(which(sp_tbl$`SFT.R.sq` >= 0.80))]
      cat("Selected soft power:", soft_pow, "\n")
      hdWGCNA::ConstructNetwork(
        srt,
        soft_power    = soft_pow,
        wgcna_name    = "hsc_modules",
        tom_outdir    = tom_dir,
        overwrite_tom = TRUE,
        deepSplit = 4,
        minModuleSize = 20
      )
    }
  ),

  tar_target(
    hdwgcna_hsc_soft_power,
    {
      sp_tbl <- hdwgcna_hsc_network@misc[["hsc_modules"]]$wgcna_powerTable
      sp_tbl$Power[min(which(sp_tbl$`SFT.R.sq` >= 0.80))]
    }
  ),

  tar_target(
    hdwgcna_hsc_modules,
    {
      options(future.globals.maxSize = 8 * 1024^3)
      tom_dir <- file.path(!!config$paths$objects_dir, "hdwgcna")
      srt <- Seurat::ScaleData(hdwgcna_hsc_network, features = hdwgcna_hsc_hvg)
      srt <- hdWGCNA::ModuleEigengenes(srt, wgcna_name = "hsc_modules")
      hdWGCNA::ModuleConnectivity(srt,
                                   wgcna_name  = "hsc_modules",
                                   path_to_TOM = tom_dir)
    }
  ),

  # ── Pairwise comparison group definitions ─────────────────────────────────────
  tar_target(hsc_pairwise_comparisons,
             list(
               mch_pos_vs_wt  = c("HSC_Het_CD41+mCh+", "HSC_WT_CD41+"),
               mch_neg_vs_wt  = c("HSC_Het_CD41+mCh-", "HSC_WT_CD41+"),
               mch_pos_vs_neg = c("HSC_Het_CD41+mCh+", "HSC_Het_CD41+mCh-")
             )
  ),
  tar_target(mnc_pairwise_comparisons,
             list(
               mch_pos_vs_wt  = c("MNC_mCh+",  "MNC_WT"),
               mch_neg_vs_wt  = c("MNC_mCh-",  "MNC_WT"),
               mch_pos_vs_neg = c("MNC_mCh+",  "MNC_mCh-")
             )
  ),
  tar_target(stroma_pairwise_comparisons,
             list(
               het_vs_wt = c("Stroma_Het", "Stroma_WT")
             )
  ),

  # ── Cluster markers: each-vs-rest per cluster at each resolution ─────────────
  # Produces: cluster_de_each_vs_rest_{compartment}_{res_label}
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("hsc",                    "hsc",                    "mnc",                    "stroma"),
      res_label   = c("0_9",                    "1_5",                    "0_9",                    "0_9"),
      cluster_res = c("sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.1.5", "sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.0.9"),
      srt_target  = rlang::syms(c("mpn_processed_hsc_fullcc", "mpn_processed_hsc_fullcc",
                                   "mpn_processed_mnc",        "mpn_processed_stroma"))
    ),
    names = c(compartment, res_label),

    tar_target(
      cluster_de_each_vs_rest,
      {
        options(future.globals.maxSize = 8 * 1024^3)
        counts_mat <- srt_target[["SCT"]]@counts
        counts_mat <- counts_mat[!rownames(counts_mat) %in% excluded_genes, ]
        DElegate::findDE(
          counts_mat,
          meta_data        = srt_target@meta.data,
          group_column     = cluster_res,
          compare_is_ref   = FALSE
        )
      }
    )
  ),

  # ── Cluster markers: all-vs-all per cluster at each resolution ────────────────
  # Produces: cluster_de_all_vs_all_{compartment}_{res_label}
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("hsc",                    "hsc",                    "mnc",                    "stroma"),
      res_label   = c("0_9",                    "1_5",                    "0_9",                    "0_9"),
      cluster_res = c("sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.1.5", "sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.0.9"),
      srt_target  = rlang::syms(c("mpn_processed_hsc_fullcc", "mpn_processed_hsc_fullcc",
                                   "mpn_processed_mnc",        "mpn_processed_stroma"))
    ),
    names = c(compartment, res_label),

    tar_target(
      cluster_de_all_vs_all,
      {
        options(future.globals.maxSize = 8 * 1024^3)
        srt2 <- srt_target
        rm(srt_target); gc()
        Seurat::Idents(srt2) <- cluster_res
        find_delegate_markers_all_clusters2(srt2, excluded_genes, group_column = cluster_res)
      }
    )
  ),

  # ── Condition DE per cluster: pairwise at each resolution ────────────────────
  # Produces: condition_de_per_cluster_{compartment}_{res_label}
  # Returns list(comparison_name = df, ...) — one entry per comparison in
  # pairwise_comparisons; the outer cluster_levels wrapper is unwrapped via [[1L]].
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("hsc",                    "hsc",                    "mnc",                    "stroma"),
      res_label   = c("0_9",                    "1_5",                    "0_9",                    "0_9"),
      cluster_res = c("sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.1.5", "sNN_pca.SCT_30_res.0.9", "sNN_pca.SCT_30_res.0.9"),
      srt_target  = rlang::syms(c("mpn_processed_hsc_fullcc", "mpn_processed_hsc_fullcc",
                                   "mpn_processed_mnc",        "mpn_processed_stroma")),
      comparisons = rlang::syms(c("hsc_pairwise_comparisons",  "hsc_pairwise_comparisons",
                                   "mnc_pairwise_comparisons",  "stroma_pairwise_comparisons"))
    ),
    names = c(compartment, res_label),

    tar_target(
      condition_de_per_cluster,
      {
        options(future.globals.maxSize = 8 * 1024^3)
        clust_vec <- setNames(
          as.character(srt_target@meta.data[[cluster_res]]),
          colnames(srt_target)
        )
        find_delegate_pairwise(
          srt            = srt_target,
          excluded_genes = excluded_genes,
          comparisons    = comparisons,
          cluster_levels = list(res = clust_vec)
        )[[1L]]
      }
    )
  ),

  # ── Moran's I per sample.name condition — HSC ────────────────────────────────
  # Computes graph-based spatial autocorrelation for variable features within
  # each condition subset using the SNN graph as spatial weights.
  # moran_sparse() stays fully sparse — no n² dense matrix allocation.
  tarchetypes::tar_map(
    values = tibble::tibble(
      condition_name  = c("neg",       "pos",       "wt"),
      condition_label = c("HSC_Het_CD41+mCh-", "HSC_Het_CD41+mCh+", "HSC_WT_CD41+")
    ),
    names = condition_name,

    tar_target(
      moransi_hsc,
      {
        srt_sub <- subset(mpn_processed_hsc_fullcc, sample.name == condition_label)
        pos     <- Embeddings(srt_sub, "umap.nn")
        genes   <- VariableFeatures(mpn_processed_hsc_fullcc)
        data    <- as.matrix(GetAssayData(srt_sub, assay = "SCT", layer = "data")[genes, ])
        res     <- Seurat::RunMoransI(data = data, pos = pos)
        tibble::tibble(
          gene      = rownames(res),
          morans_i  = res$observed,
          p_value   = res$p.value,
          p_adj     = p.adjust(res$p.value, method = "BH"),
          condition = condition_label
        ) |> dplyr::arrange(dplyr::desc(morans_i))
      }
    ),

    tar_target(
      moransi_hsc_path,
      save_csv(moransi_hsc,
               !!config$paths$output_table_dir,
               paste0("moransi_hsc_", condition_name, ".csv")),
      format = "file"
    )
  ),

  # ── Moran's I per sample.name condition — HSC (PCA, first 10 PCs) ───────────
  # Same as moransi_hsc but uses PCA coordinates (dims 1–10) as spatial weights
  # instead of UMAP, giving a higher-dimensional neighbourhood signal.
  tarchetypes::tar_map(
    values = tibble::tibble(
      condition_name  = c("neg",                "pos",                "wt"),
      condition_label = c("HSC_Het_CD41+mCh-",  "HSC_Het_CD41+mCh+",  "HSC_WT_CD41+")
    ),
    names = condition_name,

    tar_target(
      moransi_pca_hsc,
      {
        srt_sub <- subset(mpn_processed_hsc_fullcc, sample.name == condition_label)
        pos     <- Embeddings(srt_sub, "pca.SCT")[, 1:10]
        genes   <- VariableFeatures(mpn_processed_hsc_fullcc)
        data    <- as.matrix(GetAssayData(srt_sub, assay = "SCT", layer = "data")[genes, ])
        res     <- Seurat::RunMoransI(data = data, pos = pos)
        tibble::tibble(
          gene      = rownames(res),
          morans_i  = res$observed,
          p_value   = res$p.value,
          p_adj     = p.adjust(res$p.value, method = "BH"),
          condition = condition_label
        ) |> dplyr::arrange(dplyr::desc(morans_i))
      }
    ),

    tar_target(
      moransi_pca_hsc_path,
      save_csv(moransi_pca_hsc,
               !!config$paths$output_table_dir,
               paste0("moransi_pca_hsc_", condition_name, ".csv")),
      format = "file"
    )
  ),

  # ── SCENT potency scores ──────────────────────────────────────────
  # Compute signaling entropy potency scores per compartment.
  # Returns data.frame with columns: cell.name, potency, scaled_potency.
  tar_target(scent_ppi_net,rename_scent_net()),
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",                    "hsc",                    "mnc"),
      srt_target = rlang::syms(c("mpn_processed_data_stroma", "mpn_processed_data_hsc", "mpn_processed_data_mnc"))
    ),
    names = compartment,

    tar_target(
      scent_scores,
      get_scent_potency(srt_target, scent_ppi_net)
    )
  ),

  # ── SCENT metadata + scores ────────────────────────────────────────
  # Lightweight dataframe with metadata and potency scores for visualization.
  # Columns: cell.name, sample.name, nFeature_RNA, potency, scaled_potency
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment  = c("stroma",                    "hsc",                    "mnc"),
      srt_target   = rlang::syms(c("mpn_processed_data_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_data_mnc")),
      scent_target = rlang::syms(c("scent_scores_stroma", "scent_scores_hsc", "scent_scores_mnc"))
    ),
    names = compartment,

    tar_target(
      scent_meta,
      {
        tibble::tibble(
          cell.name     = colnames(srt_target),
          sample.name   = srt_target$sample.name,
          nFeature_RNA  = srt_target$nFeature_RNA
        ) %>%
          dplyr::left_join(scent_target, by = "cell.name")
      }
    )
  ),
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(cytotrace,
               CytoTRACE2::cytotrace2(srt_target@assays$SCT@counts,
                                      species = "mouse",
                                      ncores = 1)
               )
    ),
  tar_target(cytotrace_hsc_wt,
             CytoTRACE2::cytotrace2(mpn_processed_hsc_wt_fullcc@assays$SCT@counts,
                                    species = "mouse",
                                    ncores = 1)
             ),
  # ── MNC QC curve clusters ────────────────────────────────────────
  # Cluster MNC cells by nCount_RNA, nFeature_RNA, percent.rp to capture
  # the two complexity curves visible in the feature scatter.
  # Returns a named character vector: cluster label per cell.
  tar_target(
    mnc_qc_curves_clusters,
    {
      cols <- c("nCount_RNA", "nFeature_RNA", "percent.rp")
      mat  <- scale(as.matrix(mpn_raw_processed_mnc@meta.data[, cols, drop = FALSE]))
      rownames(mat) <- colnames(mpn_raw_processed_mnc)
      colnames(mat) <- paste0("qc", seq_along(cols))
      srt <- mpn_raw_processed_mnc
      srt[["qc.pca"]] <- Seurat::CreateDimReducObject(
        embeddings = mat, key = "QC_", assay = DefaultAssay(srt))
      srt <- Seurat::FindNeighbors(srt, reduction = "qc.pca",
                                   dims = seq_len(ncol(mat)), k.param = 100,
                                   graph.name = c("qc_kNN", "qc_sNN"))
      srt <- Seurat::FindClusters(srt, graph.name = "qc_sNN",
                                  resolution = 0.1, group.singletons = TRUE)
      stats::setNames(as.character(srt$seurat_clusters), colnames(srt))
    }
  ),

  # ── QC output ───────────────────────────────────────────────────
  # metadata_flagged_all: combined metadata table across all compartments,
  # with QC metrics and pass/fail flags. Lightweight target (data.frame only).
  tar_target(
    metadata_flagged_all,
    rbind(mpn_raw_stroma@meta.data,
          mpn_raw_hsc@meta.data,
          mpn_raw_mnc@meta.data)
  ),

  tar_target(
    qc_summary,
    {
      before   <- table(mpn_raw_all$sample.name)
      pass_tbl <- table(metadata_flagged_all$sample.name[metadata_flagged_all$pass_cell])
      data.frame(
        sample        = names(before),
        cells_before  = as.integer(before),
        cells_after   = as.integer(pass_tbl[names(before)]),
        cells_removed = as.integer(before) - as.integer(pass_tbl[names(before)])
      )
    }
  ),

  tar_target(
    raw_cell_summary,
    data.frame(
      sample  = names(table(mpn_raw_all$sample.name)),
      n_cells = as.integer(table(mpn_raw_all$sample.name)),
      n_genes = nrow(mpn_raw_all)
    )
  ),

  tar_target(qc_stats_summary, summarise_qc_metrics(metadata_flagged_all)),

  tar_target(doublet_summary,  summarise_doublets(metadata_flagged_all)),

  tar_target(qc_summary_path,
             save_csv(qc_summary,       !!config$paths$output_table_dir, "qc_summary.csv"),
             format = "file"),

  tar_target(raw_cell_summary_path,
             save_csv(raw_cell_summary, !!config$paths$output_table_dir, "raw_cell_summary.csv"),
             format = "file"),

  tar_target(qc_stats_summary_path,
             save_csv(qc_stats_summary, !!config$paths$output_table_dir, "qc_stats_summary.csv"),
             format = "file"),

  tar_target(doublet_summary_path,
             save_csv(doublet_summary,  !!config$paths$output_table_dir, "doublet_summary.csv"),
             format = "file"),

  # ── Sub-clustering ───────────────────────────────────────────────
  # Driven by `subclusters:` in metadata/project_config.yaml.  Each entry defines:
  #   name:        short identifier (used in target names + output paths)
  #   compartment: stroma | hsc | mnc  → sources mpn_processed_<compartment>
  #   resolution:  cluster resolution to subset on (e.g. 0.1)
  #   clusters:    list of cluster IDs to keep (e.g. ["5", "9"])
  #
  # Produces per entry: mpn_subclust_<name>  (subsetted + reprocessed Seurat)
  # To add a run: append an entry under subclusters: in metadata/project_config.yaml.

  tarchetypes::tar_map(
    values = {
      entries <- config$subclusters
      tibble::tibble(
        name        = vapply(entries, `[[`, "", "name"),
        compartment = tolower(vapply(entries, `[[`, "", "compartment")),   # config says HSC; switch() below uses hsc
        resolution  = vapply(entries, function(x) as.numeric(x$resolution), 0),
        clusters    = lapply(entries, function(x) as.character(unlist(x$clusters)))
      )
    },
    names  = name,

    tar_target(
      mpn_subclust,
      {
        src <- switch(compartment,
                      stroma = mpn_processed_stroma,
                      hsc    = mpn_processed_hsc,
                      mnc    = mpn_processed_mnc,
                      stop("Unknown compartment: ", compartment))
        sub <- subset_clusters(src, clusters = clusters, resolution = resolution)
        rm(src); gc()
        process_all(sub,
                    ref_assay     = !!config$reduction$ref_assay,
                    dims_i        = seq(!!config$reduction$n_pca_dims),
                    k             = !!config$reduction$k_neighbors,
                    res           = unlist(!!config$clustering$resolutions),
                    n_features    = !!config$reduction$n_features,
                    k.nn          = !!config$reduction$n_umap_neighbors,
                    exclude_genes = excluded_genes)
      }
    )
  ),

  # ── Reference datasets ──────────────────────────────────────────
  # Each reference is split into two targets:
  #   <name>           — processed Seurat object (format = "qs")
  #   <name>_azmth_path — Azimuth folder built from it (format = "file")
  # Projection targets receive <name>_azmth_path as ref_dir.
  #
  # Tabula Muris Senis (cellxgene.census, SCT only):
  #   Marrow + Spleen: biologically relevant; Lung + Kidney: negative controls.
  # tar_target(tabula_muris_marrow,
  #   seurat_ref_process(fetch_ref_census("0bd1a1de-3aee-40e0-b2ec-86c7a30c7149")),
  #   format = "qs"),
  # tar_target(tabula_muris_spleen,
  #   seurat_ref_process(fetch_ref_census("e2b469d4-b5c3-4a35-9d19-ee71ce61cae0")),
  #   format = "qs"),
  # tar_target(tabula_muris_lung,
  #   seurat_ref_process(fetch_ref_census("0fb7916e-7a68-4a4c-a441-3ab3989f29a7")),
  #   format = "qs"),
  # tar_target(tabula_muris_kidney,
  #   seurat_ref_process(fetch_ref_census("b8c618e5-4b3d-4566-8a3f-7e40047f5c54")),
  #   format = "qs"),
  # Mouse Bone Marrow Cell Atlas (syn53222728): pre-built refDR/refUMAP/refAssay.
  tar_target(ferchen_26_mm_bm_path, "data/ferchen_2026_nat_immunol/azimuth/", format = "file"),
  tar_target(ferchen_26_mm_bm_umap, extract_azimuth_umap_df(ferchen_26_mm_bm_path)),

  # MCA references — azimuth/ subdirectory built by data/mca_mouse_cell_atlas/10x_to_azimuth.R
  tar_target(mca_mm_pbmc_path,       "data/mca_mouse_cell_atlas/Peripheral-Blood/azimuth/",     format = "file"),
  tar_target(mca_mm_bm_path,         "data/mca_mouse_cell_atlas/Bone-Marrow/azimuth/",          format = "file"),
  tar_target(mca_mm_bm_ckit_path,    "data/mca_mouse_cell_atlas/Bone-Marrow_c-kit/azimuth/",    format = "file"),
  tar_target(mca_mm_bm_mesench_path, "data/mca_mouse_cell_atlas/Embryonic-Mesenchyme/azimuth/", format = "file"),
  tar_target(mca_mm_pbmc_umap,       extract_azimuth_umap_df(mca_mm_pbmc_path)),
  tar_target(mca_mm_bm_umap,         extract_azimuth_umap_df(mca_mm_bm_path)),
  tar_target(mca_mm_bm_ckit_umap,    extract_azimuth_umap_df(mca_mm_bm_ckit_path)),
  tar_target(mca_mm_bm_mesench_umap, extract_azimuth_umap_df(mca_mm_bm_mesench_path)),

  # zeng_25_hm_bm — Zeng 2025 BoneMarrowMap human reference (cross-species projection).
  tar_target(zeng_25_hs_bm_path,
             "data/zeng_2025_blood_cancer_discov/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds",
             format = "file"
             ),
  tar_target(zeng_25_hs_bm_umap,
             {
               srt <- readRDS(zeng_25_hs_bm_path)
               srt %>%
                 Seurat::Embeddings("umap") %>%
                 as.data.frame() %>%
                 rownames_to_column("cell_id") %>%
                 left_join(srt@meta.data %>% rownames_to_column("cell_id"))
             }
             ),
  #Hematlas: bone-marrow (bm), aorta-gonad-mesonephros (agm), and fetal liver (fl)
  tar_target(hematlas_mm_bm_path, "data/hematlas/mouse_bm/azimuth/", format = "file"),
  tar_target(hematlas_mm_agm_path, "data/hematlas/mouse_agm/azimuth/", format = "file"),
  tar_target(hematlas_mm_fl_path, "data/hematlas/mouse_fl/azimuth/", format = "file"),
  tar_target(hematlas_mm_bm, extract_azimuth_umap_df(hematlas_mm_bm_path)),
  tar_target(hematlas_mm_agm, extract_azimuth_umap_df(hematlas_mm_agm_path)),
  tar_target(hematlas_mm_fl, extract_azimuth_umap_df(hematlas_mm_fl_path)),

    # PanSci references: stromal and immune WT subsets with Sub_cell_type_organ labels.
  # tar_target(pansci_stromal_rds, "external_data/pansci_stromal_wt.rds", format = "file"),
  # tar_target(pansci_stromal,
  #   load_ref_pansci(pansci_stromal_rds, label = "PanSci Stromal WT"),
  #   format = "qs"),
  # tar_target(pansci_immune_rds, "external_data/pansci_immune_wt.rds", format = "file"),
  # tar_target(pansci_immune,
  #   load_ref_pansci(pansci_immune_rds, label = "PanSci Immune WT"),
  #   format = "qs"),
  # ── Annotation: reclustering at multiple resolutions ────────────
  # Runs FindClusters on the cached SNN graph at each resolution.
  # Output: CSV with cell_id + one column per resolution (res_0_1, res_0_2, ...).
  tar_target(purity_cluster_res,
             c(0.5, 1, 1.5, 2, 2.5, 3,5)
             ),
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,

    tar_target(
      annotation_recluster,
      recluster_resolutions_df(srt_target,
                               resolutions = purity_cluster_res)
    ),

    tar_target(
      annotation_recluster_path,
      save_csv(annotation_recluster,
               !!config$paths$output_table_dir,
               paste0("annotation_recluster_", compartment, ".csv")),
      format = "file"
    ),

    tar_target(
      condition_purity,
      compute_condition_purity(annotation_recluster, srt_target, compartment, "sample.name")
    ),

    tar_target(
      condition_purity_path,
      save_csv(condition_purity,
               !!config$paths$output_table_dir,
               paste0("condition_purity_", compartment, ".csv")),
      format = "file"
    )
  ),

  # ── Purity group labels per cell — HSC ──────────────────────────────────────
  # At the chosen resolution, clusters above purity_cutoff are labelled by their
  # dominant condition; the rest are labelled "Impure".
  # Output includes UMAP1/UMAP2 so the notebook needs no Seurat object.
  tar_target(
    purity_group_hsc,
    assign_purity_groups(
      condition_purity  = condition_purity_hsc,
      recluster_df      = annotation_recluster_hsc,
      srt               = mpn_processed_hsc_fullcc,
      purity_res        = 5,
      purity_cutoff     = 0.6
    )
  ),

  tar_target(
    purity_group_hsc_path,
    save_csv(purity_group_hsc,
             !!config$paths$output_table_dir,
             "purity_group_hsc.csv"),
    format = "file"
  ),

  # ── Annotation: RunAzimuth label transfer per compartment × reference ──────
  # One tar_map per reference — compartments as rows.
  # Args exposed: annotation.levels, query.modality, assay.

  # tms_bm — Tabula Muris Senis bone marrow (biologically relevant)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_tms_bm,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "tabula_muris_marrow"))
  #     ),
  #   tar_target(proj_path_tms_bm,
  #     save_csv(proj_tms_bm, !!config$paths$output_table_dir,
  #              paste0("proj_tms_bm_", compartment, ".csv")),
  #     format = "file")
  #   ),
  #
  # # tms_sp — Tabula Muris Senis spleen (negative control)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_tms_sp,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "tabula_muris_spleen"))
  #     ),
  #   tar_target(proj_path_tms_sp,
  #     save_csv(proj_tms_sp, !!config$paths$output_table_dir,
  #              paste0("proj_tms_sp_", compartment, ".csv")), format = "file")
  # ),
  #
  # # tms_lu — Tabula Muris Senis lung (negative control)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_tms_lu,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "tabula_muris_lung"))
  #     ),
  #   tar_target(proj_path_tms_lu,
  #     save_csv(proj_tms_lu, !!config$paths$output_table_dir,
  #              paste0("proj_tms_lu_", compartment, ".csv")), format = "file")
  # ),
  #
  # # tms_ki — Tabula Muris Senis kidney (negative control)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_tms_ki,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "tabula_muris_kidney"))
  #     ),
  #   tar_target(proj_path_tms_ki,
  #     save_csv(proj_tms_ki, !!config$paths$output_table_dir,
  #              paste0("proj_tms_ki_", compartment, ".csv")), format = "file")
  # ),
  #
  # bm_atlas — Mouse Bone Marrow Cell Atlas (pre-built Azimuth reference)
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_ferchen_26_mm_bm,
               run_azimuth_df(srt_target,
                              ref_dir = ferchen_26_mm_bm_path)
    ),
    tar_target(proj_ferchen_26_mm_bm_path,
               save_csv(proj_ferchen_26_mm_bm, !!config$paths$output_table_dir,
                        paste0("proj_ferchen_26_mm_bm", compartment, ".csv")), format = "file")
  ),

  # mca_pb — MCA Peripheral Blood
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_mca_mm_pbmc,
      run_azimuth_df(srt_target,
                     ref_dir = mca_mm_pbmc_path)
      ),
    tar_target(proj_path_mca_mm_pbmc,
      save_csv(proj_mca_mm_pbmc, !!config$paths$output_table_dir,
               paste0("proj_mca_mm_pbmc_", compartment, ".csv")), format = "file")
  ),

  # mca_bm — MCA Bone Marrow
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_mca_mm_bm,
               run_azimuth_df(srt_target, ref_dir = mca_mm_bm_path)),
    tar_target(proj_path_mca_mm_bm,
               save_csv(proj_mca_mm_bm, !!config$paths$output_table_dir,
                        paste0("proj_mca_mm_bm_", compartment, ".csv")), format = "file")
  ),

  # mca_ckit — MCA BM c-Kit enriched (progenitor-enriched)
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_mca_mm_bm_ckit,
               run_azimuth_df(srt_target, ref_dir = mca_mm_bm_ckit_path)),
    tar_target(proj_path_mca_mm_bm_ckit,
               save_csv(proj_mca_mm_bm_ckit, !!config$paths$output_table_dir,
                        paste0("proj_mca_mm_bm_ckit_", compartment, ".csv")), format = "file")
  ),

  # mca_mes — MCA BM Mesenchyme (stromal/mesenchymal)
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_mca_mm_bm_mesench,
               run_azimuth_df(srt_target, ref_dir = mca_mm_bm_mesench_path)),
    tar_target(proj_path_mca_mm_bm_mesench,
               save_csv(proj_mca_mm_bm_mesench, !!config$paths$output_table_dir,
                        paste0("proj_mca_mm_bm_mesench_", compartment, ".csv")), format = "file")
  ),

  # zeng_25_hm_bm — Zeng 2025 BoneMarrowMap human reference (cross-species projection).
  # Mouse gene names are first renamed to human orthologs (orthogene::convert_orthologs),
  # then LogNormalised on the RNA assay, then projected via FindTransferAnchors.
  # Transfers CellType and CellType_Broad labels from the reference metadata.

  # Named vector: names = mouse gene symbols, values = human ortholog symbols.
  # Built once from the full gene universe; NAs are filtered inside rename_rna_feats().
  tar_target(
    ortholog_map_mm_to_hm,
    {
      orth_df <- orthogene::convert_orthologs(
        gene_df        = rownames(mpn_raw_all),
        input_species  = "mouse",
        output_species = "human",
        method         = "gprofiler",
        agg_fun = "stats"
      )
      # orth_df rownames = input mouse gene, Gene.Symbol = human ortholog
      setNames(rownames(orth_df), orth_df$input_gene)
    }
  ),

  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,

    tar_target(
      proj_zeng_25_hs_bm,
      project_mm_onto_hs_seurat(ref = readRDS(zeng_25_hs_bm_path),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(
                                  CellType  = "CellType",
                                  CellType_Broad = "CellType_Broad")
                                )
      ),
    tar_target(
      proj_path_zeng_25_hs_bm,
      save_csv(proj_zeng_25_hs_bm, !!config$paths$output_table_dir,
               paste0("proj_zeng_25_hs_bm_", compartment, ".csv")),
      format = "file"
    )
  ),

  #suo
  ##Megakaryocyte/erythroid cells
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,

    tar_target(
      proj_suo_22_hs_mem,
      project_mm_onto_hs_seurat(ref = readRDS("data/suo_2022/suo_mem_prog_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(
                                  celltype  = "celltype",
                                  organ = "organ")
      )
    ),
    tar_target(
      proj_path_suo_22_hs_mem,
      save_csv(proj_suo_22_hs_mem, !!config$paths$output_table_dir,
               paste0("proj_suo_22_hs_mem_", compartment, ".csv")),
      format = "file"
    )
  ),

  ##HSC/progenitor cells
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,

    tar_target(
      proj_suo_22_hs_hsc,
      project_mm_onto_hs_seurat(ref = readRDS("data/suo_2022/suo_hsc_prog_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(
                                  celltype  = "celltype",
                                  organ = "organ")
      )
    ),
    tar_target(
      proj_path_suo_22_hs_hsc,
      save_csv(proj_suo_22_hs_hsc, !!config$paths$output_table_dir,
               paste0("proj_suo_22_hs_hsc_", compartment, ".csv")),
      format = "file"
    )
  ),
  ##Myeloid
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,

    tar_target(
      proj_suo_22_hs_myel,
      project_mm_onto_hs_seurat(ref = readRDS("data/suo_2022/suo_myeloid_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(
                                  celltype  = "celltype",
                                  organ = "organ")
      )
    ),
    tar_target(
      proj_path_suo_22_hs_myel,
      save_csv(proj_suo_22_hs_myel, !!config$paths$output_table_dir,
               paste0("proj_suo_22_hs_myel_", compartment, ".csv")),
      format = "file"
    )
  ),
  #jardine
  ## fig1b — full fetal BM atlas
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(
      proj_jardine_21_hs_fig1b,
      project_mm_onto_hs_seurat(ref = readRDS("data/jardine_2021/jardine_fig1b_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(celltype = "celltype"))
    ),
    tar_target(
      proj_path_jardine_21_hs_fig1b,
      save_csv(proj_jardine_21_hs_fig1b, !!config$paths$output_table_dir,
               paste0("proj_jardine_21_hs_fig1b_", compartment, ".csv")),
      format = "file"
    )
  ),

  ## prog — progenitors only
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(
      proj_jardine_21_hs_prog,
      project_mm_onto_hs_seurat(ref = readRDS("data/jardine_2021/jardine_prog_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(celltype = "celltype"))
    ),
    tar_target(
      proj_path_jardine_21_hs_prog,
      save_csv(proj_jardine_21_hs_prog, !!config$paths$output_table_dir,
               paste0("proj_jardine_21_hs_prog_", compartment, ".csv")),
      format = "file"
    )
  ),

  ## mnc — mononuclear cells (CITE-seq)
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(
      proj_jardine_21_hs_mnc,
      project_mm_onto_hs_seurat(ref = readRDS("data/jardine_2021/jardine_mnc_ref.rds"),
                                query = srt_target,
                                ortho_table = mm_hs_ortho_long,
                                assay = "RNA",
                                tranfer_list = list(celltype = "celltype"))
    ),
    tar_target(
      proj_path_jardine_21_hs_mnc,
      save_csv(proj_jardine_21_hs_mnc, !!config$paths$output_table_dir,
               paste0("proj_jardine_21_hs_mnc_", compartment, ".csv")),
      format = "file"
    )
  ),

  # hematlas_fl — Hematlas Fetal Liver
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_hematlas_mm_fl,
               run_azimuth_df(srt_target,ref_dir = hematlas_mm_fl_path)
    ),
    tar_target(proj_path_hematlas_mm_fl,
               save_csv(proj_hematlas_mm_fl, !!config$paths$output_table_dir,
                        paste0("proj_hematlas_mm_fl_", compartment, ".csv")), format = "file")
  ),

  # hematlas_agm — Hematlas aorta-gonad-mesonephros
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_hematlas_mm_agm,
               run_azimuth_df(srt_target, ref_dir = hematlas_mm_agm_path)),
    tar_target(proj_path_hematlas_mm_agm,
               save_csv(proj_hematlas_mm_agm, !!config$paths$output_table_dir,
                        paste0("proj_hematlas_mm_agm_", compartment, ".csv")), format = "file")
  ),

  # hematlas_bm — Hematlas Bone marrow
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(proj_hematlas_mm_bm,
               run_azimuth_df(srt_target, ref_dir = hematlas_mm_bm_path)),
    tar_target(proj_path_hematlas_mm_bm,
               save_csv(proj_hematlas_mm_bm, !!config$paths$output_table_dir,
                        paste0("proj_hematlas_mm_bm_", compartment, ".csv")), format = "file")
  ),

  # # pansci_stromal — PanSci stromal WT reference (Sub_cell_type_organ labels)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_pansci_stromal,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "pansci_stromal"))
  #     ),
  #   tar_target(proj_path_pansci_stromal,
  #     save_csv(proj_pansci_stromal, !!config$paths$output_table_dir,
  #              paste0("proj_pansci_stromal_", compartment, ".csv")), format = "file")
  # ),

  # # pansci_immune — PanSci immune WT reference (Sub_cell_type_organ labels)
  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc", "mpn_processed_mnc"))
  #   ),
  #   names = compartment,
  #   tar_target(proj_pansci_immune,
  #     run_azimuth_df(srt_target,
  #                    ref_dir = file.path(!!config$paths$azimuth_dir, "pansci_immune"))
  #     ),
  #   tar_target(proj_path_pansci_immune,
  #     save_csv(proj_pansci_immune, !!config$paths$output_table_dir,
  #              paste0("proj_pansci_immune_", compartment, ".csv")), format = "file")
  # ),

  # ── Annotation notebook flat tables ────────────────────────────────────────
  # These replace Seurat object loading in notebooks. Each target is a plain
  # data.frame (CSV-serialisable) containing only what the notebook needs.

  # Per-compartment: umap.nn coords + cluster columns + sample/HTO metadata.
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_mnc"))
    ),
    names = compartment,
    tar_target(anno_meta,   extract_anno_meta(srt_target)),
    tar_target(anno_meta_path,
      save_csv(anno_meta, !!config$paths$output_table_dir, paste0("anno_meta_", compartment, ".csv")),
      format = "file")
  ),

  # Claude-generated marker list
  tar_target(hem_markers_claude_df,
             read.csv("data/marker_list_claude/mouse_hematopoiesis_markers.csv")),
  # Per-compartment: SCT expression for all canonical marker genes.
  tar_target(canonical_markers_lst, list(!!!config$annotation$cell_types)),
  tar_target(canonical_markers_vec, canonical_markers_lst %>% unlist() %>% unique()),
  tar_target(anno_expr_stroma,
             extract_expr_df(mpn_processed_stroma,
                             genes = canonical_markers_vec)
             ),
  tar_target(anno_expr_hsc,
             extract_expr_df(mpn_processed_hsc,
                             genes = canonical_markers_vec)
             ),
  tar_target(anno_expr_mnc,
             extract_expr_df(mpn_processed_mnc,
                             genes = canonical_markers_vec)
             ),

  # Pre-filter (raw) metadata tables — umap.nn coords + PC1 + all QC columns.
  # Replaces loading mpn_raw_processed_* Seurat objects in the milestone report.
  tar_target(anno_raw_meta_stroma,
    extract_raw_meta(mpn_raw_processed_stroma, metadata_flagged_all)),
  tar_target(anno_raw_meta_hsc,
    extract_raw_meta(mpn_raw_processed_hsc,    metadata_flagged_all)),
  tar_target(anno_raw_meta_mnc,
    extract_raw_meta(mpn_raw_processed_mnc,    metadata_flagged_all)),

  # Per-reference: UMAP coords + cell_type label from the reference Seurat.
  # Tiny table (~100k rows × 3 cols) — replaces loading the full Seurat in notebooks.
  # tar_target(ref_umap_df_tms_bm,   extract_ref_umap_df(tabula_muris_marrow)),
  # tar_target(ref_umap_df_tms_sp,   extract_ref_umap_df(tabula_muris_spleen)),
  # tar_target(ref_umap_df_tms_lu,   extract_ref_umap_df(tabula_muris_lung)),
  # tar_target(ref_umap_df_tms_ki,   extract_ref_umap_df(tabula_muris_kidney)),
  # tar_target(ref_umap_df_mca_pb,   extract_ref_umap_df(mca_pbmc)),
  # tar_target(ref_umap_df_mca_bm,   extract_ref_umap_df(mca_bm)),
  # tar_target(ref_umap_df_mca_ckit, extract_ref_umap_df(mca_bm_ckit)),
  # tar_target(ref_umap_df_mca_mes,  extract_ref_umap_df(mca_bm_mesench)),

  # ── Dim-reduction parameter grid (SCT) ──────────────────────────────────────
  # Pre-computed data for the dim-reduction-and-clustering-SCT notebook.
  # $pre: elbow + PCA scatter + pipeline UMAP + Harmony UMAP (data frames).
  # $runs: lightweight data frames (UMAP + cluster cols) per dims×k combination.
  # SCTransform is called once per subset in $runs — NOT once per grid cell.

  # tarchetypes::tar_map(
  #   values = tibble::tibble(
  #     compartment = c("stroma",               "hsc",               "mnc"),
  #     srt_proc    = rlang::syms(c("mpn_processed_stroma", "mpn_processed_hsc",   "mpn_processed_mnc")),
  #     srt_clean   = rlang::syms(c("mpn_clean_stroma",     "mpn_clean_hsc",       "mpn_clean_mnc"))
  #   ),
  #   names = compartment,
  #
  #   tar_target(
  #     dim_grid_sct_pre,
  #     compute_dim_grid_sct_pre(srt_proc,
  #                               n_dims   = !!config$reduction$n_pca_dims,
  #                               k_umap   = !!config$reduction$n_umap_neighbors,
  #                               pca_name = paste0("pca.", !!config$reduction$ref_assay)),
  #     format = "qs"
  #   ),
  #
  #   tar_target(
  #     dim_grid_sct_runs,
  #     compute_dim_grid_sct_runs(srt_clean,
  #                                ref_assay     = !!config$reduction$ref_assay,
  #                                dims_grid     = unlist(!!config$dim_grid$dims),
  #                                k_grid        = unlist(!!config$dim_grid$k),
  #                                res_scan      = unlist(!!config$dim_grid$res_scan),
  #                                n_features    = !!config$reduction$n_features,
  #                                k_umap        = !!config$reduction$n_umap_neighbors,
  #                                exclude_genes = excluded_genes),
  #     format = "qs"
  #   )
  # ),

  # ── Milo differential abundance ─────────────────────────────────────────────
  # Each target stores list(milo, da_results, nhood_sizes, cluster_col).
  # The milo object is needed by plotNhoodGraphDA in the notebook.
  # cluster_col is derived once here so the notebook does not grep metadata.

  tar_target(
    hsc_milo_cluster,
    "sNN_pca.SCT_30_res.1.5"
    ),
  tar_target(milo_hsc_obj,
             .build_milo(mpn_processed_hsc_fullcc,
                         subset = "hsc",
                         graph_name = "sNN_pca.SCT_30",
                         k = !!config$reduction$k_neighbors,
                         d = !!config$reduction$n_pca_dims,
                         prop = 0.1,
                         condition = "sample.name",
                         replicates = "dm2_class",
                         dim_reduc = "pca.SCT",
                         umap = "umap.nn")
             ),
  tar_target(hsc_milo_tests,
             list(mChNeg_vs_WT = list(test_variable = "condition",
                                      test_control = "WT",
                                      test_condition = "Het_mChNeg"),
                  mChPos_vs_WT = list(test_variable = "condition",
                                      test_control = "WT",
                                      test_condition = "Het_mChPos"),
                  mChPos_vs_mChNeg = list(test_variable = "condition",
                                          test_control = "Het_mChNeg",
                                          test_condition = "Het_mChPos"))
             ),
  tar_target(milo_hsc_da,
             lapply(hsc_milo_tests, function(test_x){
               test_milo_pair(test_vars = test_x,
                              milo_hsc_obj,
                              reduced_dim = "pca.SCT",
                              cluster_col = hsc_milo_cluster,
                              condition = "sample.name",
                              norm_method = "TMM")
               } )
             ),

  # Per-cell mean logFC: average neighbourhood logFC across all nhoods each cell
  # belongs to. Cells in no neighbourhood receive NA. Returns a named list of
  # three data frames (one per comparison) with columns: cell_id, mean_lfc, n_nhoods.
  tar_target(
    milo_hsc_cell_scores,
    {
      nhood_mat    <- miloR::nhoods(milo_hsc_obj)          # cells × nhoods (sparse)
      n_membership <- Matrix::rowSums(nhood_mat)           # nhoods per cell

      lapply(milo_hsc_da, function(res) {
        da      <- res$da_results
        lfc_vec <- setNames(da$logFC, as.character(da$Nhood))
        lfc_aligned <- as.numeric(lfc_vec[as.character(seq_len(ncol(nhood_mat)))])
        lfc_aligned[is.na(lfc_aligned)] <- 0

        cell_sum  <- as.numeric(nhood_mat %*% lfc_aligned)
        cell_mean <- ifelse(n_membership > 0, cell_sum / n_membership, NA_real_)

        data.frame(
          cell_id  = rownames(nhood_mat),
          mean_lfc = cell_mean,
          n_nhoods = as.integer(n_membership)
        )
      })
    }
  ),

  # ── Stroma Milo (HSC-style: shared Milo object + per-contrast DA) ───────────
  tar_target(
    stroma_milo_cluster,
    "sNN_pca.SCT_30_res.0.9"
    ),
  tar_target(
    milo_stroma_obj,
    .build_milo(mpn_processed_stroma,
                subset     = "stroma",
                graph_name = "sNN_pca.SCT_30",
                k          = !!config$reduction$k_neighbors,
                d          = !!config$reduction$n_pca_dims,
                prop       = 0.1,
                condition  = "sample.name",
                replicates = "dm2_class",
                dim_reduc  = "pca.SCT",
                umap       = "umap.nn")
  ),
  tar_target(
    stroma_milo_tests,
    list(Het_vs_WT = list(test_variable  = "condition",
                          test_control   = "WT",
                          test_condition = "Het"))
  ),
  tar_target(
    milo_stroma_da,
    lapply(stroma_milo_tests, function(test_x) {
      test_milo_pair(test_vars   = test_x,
                     milo        = milo_stroma_obj,
                     reduced_dim = "pca.SCT",
                     cluster_col = stroma_milo_cluster,
                     condition   = "condition",
                     norm_method = "TMM")
    })
  ),

  # ── MNC Milo (HSC-style: shared Milo object + per-contrast DA) ──────────────
  tar_target(
    mnc_milo_cluster,
    "sNN_pca.SCT_30_res.0.9"
    ),
  tar_target(
    milo_mnc_obj,
    .build_milo(mpn_processed_mnc,
                subset     = "mnc",
                graph_name = "sNN_pca.SCT_30",
                k          = !!config$reduction$k_neighbors,
                d          = !!config$reduction$n_pca_dims,
                prop       = 0.1,
                condition  = "sample.name",
                replicates = "dm2_class",
                dim_reduc  = "pca.SCT",
                umap       = "umap.nn")
  ),
  tar_target(
    mnc_milo_tests,
    list(mChNeg_vs_WT     = list(test_variable  = "condition",
                                 test_control   = "WT",
                                 test_condition = "mChNeg"),
         mChPos_vs_WT     = list(test_variable  = "condition",
                                 test_control   = "WT",
                                 test_condition = "mChPos"),
         mChPos_vs_mChNeg = list(test_variable  = "condition",
                                 test_control   = "mChNeg",
                                 test_condition = "mChPos"))
  ),
  tar_target(
    milo_mnc_da,
    lapply(mnc_milo_tests, function(test_x) {
      test_milo_pair(test_vars   = test_x,
                     milo        = milo_mnc_obj,
                     reduced_dim = "pca.SCT",
                     cluster_col = mnc_milo_cluster,
                     condition   = "condition",
                     norm_method = "TMM")
    })
  ),

  # anno_meta with condition column pre-computed — used by cluster proportions
  # section of the compositional-analysis notebook.
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma", "hsc", "mnc"),
      meta_target = rlang::syms(c("anno_meta_stroma", "anno_meta_hsc", "anno_meta_mnc"))
    ),
    names = compartment,
    tar_target(
      anno_meta_cond,
      {
        df           <- add_condition_df(meta_target, compartment, "sample.name")
        cluster_col  <- grep(
          paste0("sNN_pca\\..*_res\\.",
                 gsub("\\.", "\\\\.", !!config$markers$resolution), "$"),
          colnames(df), value = TRUE)[1]
        names(df)[names(df) == cluster_col] <- "cluster"
        df
      }
    )
  ),
  #miloDE
  tar_target(milo_within_ngb_pos_de,
             miloDE::de_test_neighbourhoods(milo_hsc_obj[VariableFeatures(mpn_processed_hsc_fullcc),] ,
                                    sample_id = "dm2_class",
                                    design = ~0+condition,
                                    correct_pvalues = FALSE,
                                    contrasts = c("conditionHet_mChPos - conditionWT"),
                                    covariates = c("condition"))
             ),
  tar_target(milo_within_ngb_neg_de,
             miloDE::de_test_neighbourhoods(milo_hsc_obj[VariableFeatures(mpn_processed_hsc_fullcc),] ,
                                            sample_id = "dm2_class",
                                            design = ~0+condition,
                                            correct_pvalues = FALSE,
                                            contrasts = c("conditionHet_mChNeg - conditionWT"),
                                            covariates = c("condition"))
             ),
  tar_target(milo_within_ngb_pos_neg_de,
             miloDE::de_test_neighbourhoods(milo_hsc_obj[VariableFeatures(mpn_processed_hsc_fullcc),] ,
                                            sample_id = "dm2_class",
                                            design = ~0+condition,
                                            correct_pvalues = FALSE,
                                            contrasts = c("conditionHet_mChPos - conditionHet_mChNeg"),
                                            covariates = c("condition"))
             ),

  # ── DE notebook: allvsall genesets, ORA, module scores, feature expr ────────
  # These replace computation that was previously done inline in
  # notebooks/differential-expression.qmd.

  # Per-subset: scatter_df — significant genes from pairwise DE comparisons.
  # Significance: padj <= padj_thresh & log_fc >= lfc_thresh & rate1 >= rate_cutoff (up)
  #            OR padj <= padj_thresh & log_fc <= -lfc_thresh & rate2 >= rate_cutoff (down)
  # scatter_df_hsc: named list with one entry per resolution (res_0.9, res_1.5).
  # mch_neg_vs_wt + mch_pos_vs_wt are stacked so make_scatter_df can compute
  # the two-comparison scatter axes.
  tar_target(
    scatter_df_hsc,
    list(
      res_0.9 = make_scatter_df(
        dplyr::bind_rows(condition_de_per_cluster_hsc_0_9$mch_neg_vs_wt,
                         condition_de_per_cluster_hsc_0_9$mch_pos_vs_wt) |>
          dplyr::filter(!feature %in% sex_degs),
        group_neg   = "HSC_Het_CD41+mCh-",
        group_ref   = "HSC_WT_CD41+",
        group_pos   = "HSC_Het_CD41+mCh+",
        padj_thresh = 0.001,
        lfc_thresh  = 0.5,
        rate_cutoff = 0.2
      ),
      res_1.5 = make_scatter_df(
        dplyr::bind_rows(condition_de_per_cluster_hsc_1_5$mch_neg_vs_wt,
                         condition_de_per_cluster_hsc_1_5$mch_pos_vs_wt) |>
          dplyr::filter(!feature %in% sex_degs),
        group_neg   = "HSC_Het_CD41+mCh-",
        group_ref   = "HSC_WT_CD41+",
        group_pos   = "HSC_Het_CD41+mCh+",
        padj_thresh = 0.001,
        lfc_thresh  = 0.5,
        rate_cutoff = 0.2
      )
    )
  ),

  # RETIRED (not ported as active targets): scatter_df_mnc / scatter_df_stroma read
  # delegate_allvsall_markers_mnc / _stroma, which no longer exist in this pipeline, so they and
  # everything built on them (genesets_*, ora_results_*, biased_feature_expr_* for MNC and Stroma)
  # would fail on a fresh run. scatter_df_mnc_09 / scatter_df_stroma_09 below are the current versions;
  # re-point these targets to them to bring them back. Only notebooks/archive/degs_gsea.qmd used them.
  # tar_target(
  #   scatter_df_mnc,
  #   make_scatter_df(
  #     delegate_allvsall_markers_mnc,
  #     group_neg   = "MNC_mCh-",
  #     group_ref   = "MNC_WT",
  #     group_pos   = "MNC_mCh+",
  #     padj_thresh = 0.001,
  #     lfc_thresh  = 1,
  #     rate_cutoff = 0.2
  #   )
  # ),

  # tar_target(
  #   scatter_df_stroma,
  #   make_scatter_df(
  #     delegate_allvsall_markers_stroma,
  #     group_neg   = "Stroma_Het",
  #     group_ref   = "Stroma_WT",
  #     padj_thresh = 0.001,
  #     lfc_thresh  = 1,
  #     rate_cutoff = 0.2
  #   )
  # ),

  # Per-subset scatter at res 0.9 — consumed by the DE notebooks for MNC/Stroma.
  tar_target(
    scatter_df_mnc_09,
    make_scatter_df(
      dplyr::bind_rows(
        condition_de_per_cluster_mnc_0_9$mch_neg_vs_wt,
        condition_de_per_cluster_mnc_0_9$mch_pos_vs_wt
      ) |> dplyr::filter(!feature %in% sex_degs),
      group_neg   = "MNC_mCh-",
      group_ref   = "MNC_WT",
      group_pos   = "MNC_mCh+",
      padj_thresh = 0.001,
      lfc_thresh  = 1,
      rate_cutoff = 0.2
    )
  ),

  tar_target(
    scatter_df_stroma_09,
    make_scatter_df(
      condition_de_per_cluster_stroma_0_9$het_vs_wt |>
        dplyr::filter(!feature %in% sex_degs),
      group_neg   = "Stroma_Het",
      group_ref   = "Stroma_WT",
      padj_thresh = 0.001,
      lfc_thresh  = 1,
      rate_cutoff = 0.2
    )
  ),

  # Per-subset: genesets — quadrant/band classification for two-comparison subsets,
  # up_Het/down_Het for stroma.  Returns list(genesets_df, genesets, universe).
  tar_target(
    genesets_hsc,
    make_genesets(scatter_df_hsc$res_0.9)
  ),

  # tar_target(
  #   genesets_mnc,
  #   make_genesets(scatter_df_mnc)
  # ),

  # tar_target(
  #   genesets_stroma,
  #   make_genesets(scatter_df_stroma)
  # ),

  # Per-subset: ORA results (fgsea::fora against MSigDB).
  tar_target(
    ora_results_hsc,
    compute_ora_results(genesets_hsc$genesets, genesets_hsc$universe)
  ),

  # tar_target(
  #   ora_results_mnc,
  #   compute_ora_results(genesets_mnc$genesets, genesets_mnc$universe)
  # ),

  # tar_target(
  #   ora_results_stroma,
  #   compute_ora_results(genesets_stroma$genesets, genesets_stroma$universe)
  # ),
  # Per-subset: expression of DE-notable genes for UMAP feature plots.
  # Gene lists are fixed from the each_vs_rest section of the DE notebook.
  tar_target(
    de_feature_genes_hsc,
    c("Xist", "Uty", "Kdm5d", "Ddx3y", "Xlr", "Igkc", "Iglc1",
      "Dynlt1b", "Mmrn1", "Itga2b", "Dnm3", "Hells", "Meis1", "Meis2")
  ),

  tar_target(
    de_feature_genes_mnc,
    c("Xist", "Uty", "Kdm5d", "Ddx3y", "Xlr")
  ),

  tar_target(
    de_feature_genes_stroma,
    c("Xist", "Uty", "Kdm5d", "Ddx3y", "Xlr", "Snca", "Nrgn", "Car1")
  ),

  tar_target(
    de_feature_expr_hsc,
    extract_feature_expr_long(mpn_processed_hsc, de_feature_genes_hsc)
  ),

  tar_target(
    de_feature_expr_mnc,
    extract_feature_expr_long(mpn_processed_mnc, de_feature_genes_mnc)
  ),

  tar_target(
    de_feature_expr_stroma,
    extract_feature_expr_long(mpn_processed_stroma, de_feature_genes_stroma)
  ),

  # Per-subset: expression of top-20 biased genes per geneset for biased-expr UMAPs.
  # Gene list is data-driven (top-20 off-diagonal by distance from origin).
  tar_target(
    biased_feature_expr_hsc,
    {
      biased_genes <- genesets_hsc$genesets_df %>%
        dplyr::filter(grepl("biased", band)) %>%
        dplyr::group_by(geneset) %>%
        dplyr::slice_max(sqrt(log_fc_neg^2 + log_fc_pos^2), n = 20, with_ties = FALSE) %>%
        dplyr::ungroup()
      if (nrow(biased_genes) == 0) return(data.frame())
      genes <- unique(biased_genes$feature)
      df    <- extract_feature_expr_long(mpn_processed_hsc, genes)
      merge(df, dplyr::select(biased_genes, feature, geneset),
            by.x = "gene", by.y = "feature", all.x = TRUE)
    }
  ),

  # tar_target(
  #   biased_feature_expr_mnc,
  #   {
  #     biased_genes <- genesets_mnc$genesets_df %>%
  #       dplyr::filter(grepl("biased", band)) %>%
  #       dplyr::group_by(geneset) %>%
  #       dplyr::slice_max(sqrt(log_fc_neg^2 + log_fc_pos^2), n = 20, with_ties = FALSE) %>%
  #       dplyr::ungroup()
  #     if (nrow(biased_genes) == 0) return(data.frame())
  #     genes <- unique(biased_genes$feature)
  #     df    <- extract_feature_expr_long(mpn_processed_mnc, genes)
  #     merge(df, dplyr::select(biased_genes, feature, geneset),
  #           by.x = "gene", by.y = "feature", all.x = TRUE)
  #   }
  # ),

  # tar_target(
  #   biased_feature_expr_stroma,
  #   {
  #     biased_genes <- genesets_stroma$genesets_df %>%
  #       dplyr::group_by(geneset) %>%
  #       dplyr::slice_max(abs(log_fc), n = 20, with_ties = FALSE) %>%
  #       dplyr::ungroup()
  #     if (nrow(biased_genes) == 0) return(data.frame())
  #     genes <- unique(biased_genes$feature)
  #     df    <- extract_feature_expr_long(mpn_processed_stroma, genes)
  #     merge(df, dplyr::select(biased_genes, feature, geneset),
  #           by.x = "gene", by.y = "feature", all.x = TRUE)
  #   }
  # ),
  # ── SuperCell — HSC ──────────────────────────────────────────────────────────
  # Cells are split by sample.name so no supercell spans multiple samples.
  # gamma = 10 → ~10 cells per supercell on average.
  tar_target(
    supercells_hsc,
    compute_supercells(mpn_processed_data_hsc,
                       gamma      = 50,
                       k.knn      = 5,
                       n.pc       = !!config$reduction$n_pca_dims,
                       pca_name   = paste0("pca.", !!config$reduction$ref_assay),
                       split_by   = "sample.name",
                       fast.merge = TRUE)
  ),

  tar_target(
    supercell_meta_hsc,
    {
      marker_res <- !!config$markers$resolution
      clust_col  <- grep(paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", marker_res), "$"),
                         colnames(mpn_processed_data_hsc@meta.data), value = TRUE)[1]
      extract_supercell_meta(supercells_hsc, mpn_processed_data_hsc, clust_col)
    }
  ),

  tar_target(
    supercell_ge_hsc,
    extract_supercell_ge(supercells_hsc, mpn_processed_data_hsc,
                         supercell_meta_hsc,
                         assay = !!config$reduction$ref_assay)
  ),

  tar_target(mart_hs, read.delim("data/ensembl/mart_hs_18_04_26_ensembl115.txt",sep = ",")),
  tar_target(mart_ms_hs_ortho, read.delim("data/ensembl/mart_ms_hs_ortho_18_04_26_ensembl115.txt",sep = ",")),
  tar_target(mm_hs_ortho_ensembl,
             inner_join(mart_ms_hs_ortho,mart_hs) %>%
               dplyr::select(Mouse.gene.name,
                             Gene.name,
                             Gene.Synonym,
                             HGNC.symbol) %>%
               pivot_longer(cols = c(Gene.name, Gene.Synonym, HGNC.symbol),
                            values_to = "human_ortho") %>%
               dplyr::select(Mouse.gene.name,
                             human_ortho) %>%
               distinct_all()),
  tar_target(mm_hs_ortho_mgi, read.delim("data/mgi/HOM_MouseHumanSequence.rpt") %>%
               dplyr::select(DB.Class.Key, Symbol, Mouse.MGI.ID ,   HGNC.ID) %>%
               tidyr::unite("symbol_id",Mouse.MGI.ID, HGNC.ID,sep = "")  %>%
               tidyr::separate(symbol_id,
                               into = c("species",NULL),sep = ":") %>%
               group_by(DB.Class.Key) %>%
               ungroup() %>%
               pivot_wider(names_from = "species",
                           values_from = "Symbol") %>%
               unnest(HGNC)%>%
               unnest(MGI) %>%
               dplyr::select(Mouse.gene.name = MGI ,
                             human_ortho = HGNC) ),
  tar_target(mm_hs_ortho_long, full_join(mm_hs_ortho_mgi, mm_hs_ortho_ensembl)),
  # As CSV for src/python/analysis/perturbgen_run.py
  tar_target(mm_hs_ortho_long_path,
             save_csv(mm_hs_ortho_long, !!config$paths$output_table_dir, "mm_hs_ortho_long.csv"),
             format = "file"),

  # ── External geneset resources ──────────────────────────────────────────────
  # One row per resource — consumed by the geneset-exploration notebook index.
  tar_target(
    geneset_study_table,
    tibble::tribble(
      ~study,                             ~year, ~species, ~tissue,
      "Nestorowa et al.",                 2017,  "Mouse",  "BM HSC",
      "Kucinski et al.",                  2024,  "Mouse",  "BM HSPC",
      "Isobe et al.",                     2023,  "Mouse",  "BM HSPC",
      "HemAtlas — HSPC subclusters",     2024,  "Mouse",  "BM HSPC",
      "HemAtlas — Mouse lineage ref.",   2024,  "Mouse",  "BM multi-lineage",
      "Ferchen et al. — marker genes",   2026,  "Mouse",  "BM",
      "Ferchen et al. — TF activities",  2026,  "Mouse",  "BM",
      "Gavish et al. — marker genes",           2023,  "Mouse",  "BM",
      "Gavish et al. — TF activities",          2023,  "Mouse",  "BM",
      "Psaila & Wang et al. (hs→mm)",           2020,  "Human",  "HSC/MPN",
      "Pan et al. — NMF programs (hs→mm)",      2025,  "Human",  "HSPC",
      "Gavish et al. — SCENIC regulators (hs→mm)", 2023, "Human", "BM",
      "Gavish et al. — meta programs (hs→mm)",  2023,  "Human",  "BM",
      "Gavish et al. — cancer MP & NMF (hs→mm)",2023,  "Human",  "BM"
    )
  ),

  tar_target(geneset_study_table_path,
             save_csv(geneset_study_table, !!config$paths$output_table_dir,
                      "geneset_study_table.csv"),
             format = "file"),

  # File targets — re-invalidate geneset targets if any CSV is edited.
  tar_target(nestorowa_17_bl_mm_csv,   "data/nestorowa_2017_blood/lineage_genes_tidy.csv",                                      format = "file"),
  tar_target(kucinski_24_csc_mm_csv,   "data/kucinski_2024_cell_stem_cell/hspc_marker_genes.csv",                               format = "file"),
  tar_target(isobe_23_cg_mm_csv,       "data/isobe_2023_cell_genom/trajectory_lineage_markers.csv",                             format = "file"),
  tar_target(hematlas_hspc_mm_csv,     "data/hematlas/hspc_subcluster_analysis_package/hspc_subcluster_markers.csv",            format = "file"),
  tar_target(hematlas_lineage_mm_csv,  "data/hematlas/mouse_lineage_reference_package/mouse_cross_stage_reference_markers.csv", format = "file"),
  tar_target(ferchen_26_mm_marker_csv, "data/ferchen_2026_nat_immunol/markers_marker_genes_tidy.csv",                           format = "file"),
  tar_target(ferchen_26_mm_tf_csv,     "data/ferchen_2026_nat_immunol/markers_TF_activities_tidy.csv",                          format = "file"),
  tar_target(gavish_23_mm_marker_csv,       "data/gavish_2023_nature/markers_marker_genes_tidy.csv",                           format = "file"),
  tar_target(gavish_23_mm_tf_csv,           "data/gavish_2023_nature/markers_TF_activities_tidy.csv",                          format = "file"),
  tar_target(psaila_20_hs_csv,              "data/psaila_wang_2020_molcell/Psaila2020_SupplTable5_gene_signatures_tidy.csv",   format = "file"),
  tar_target(pan_25_hs_csv,                 "data/pan_2025_cell_rep/hspc_nmf_programs_tidy.csv",                               format = "file"),
  tar_target(gavish_23_hs_scenic_csv,       "data/gavish_2023_nature/mp_scenic_regulators_tidy.csv",                          format = "file"),
  tar_target(gavish_23_hs_meta_program_csv, "data/gavish_2023_nature/meta_programs_tidy.csv",                                  format = "file"),
  tar_target(gavish_23_hs_cancer_mp_csv,    "data/gavish_2023_nature/cancer_mp_and_nmf_programs_tidy.csv",                    format = "file"),
  tar_target(yoshida_2020_ocr_csv,          "data/yoshida_2020_cell/ImmGenATAC18_AllOCRsInfo.csv",                             format = "file"),
  tar_target(hpa_secretome_tsv,            "data/human_protein_atlas/secretome.tsv",                                          format = "file"),

  # Human Protein Atlas secretome — human gene symbols with secretion classification.
  # ~2520 genes covering all secretome classes; two annotation columns retained:
  #   secretome_location: "Secreted to blood", "Secreted to extracellular matrix",
  #                        "Secreted in other tissues", "Intracellular and membrane", etc.
  #   secretome_function: "Enzyme inhibitor", "Enzyme", "No annotated function", etc.
  tar_target(
    hpa_secretome,
    readr::read_tsv(hpa_secretome_tsv, show_col_types = FALSE) %>%
      dplyr::select(
        gene_hs            = Gene,
        secretome_location = `Secretome location`,
        secretome_function = `Secretome function`
      ) %>%
      dplyr::filter(!is.na(gene_hs))
  ),

  # Mouse orthologs of HPA secretome genes mapped via mm_hs_ortho_long.
  # Use secretome_location to subset to a class of interest downstream:
  #   secreted  <- hpa_secretome_mm %>% filter(!secretome_location %in% "Intracellular and membrane")
  tar_target(
    hpa_secretome_mm,
    dplyr::inner_join(
      hpa_secretome,
      mm_hs_ortho_long,
      by = c("gene_hs" = "human_ortho")
    ) %>%
      dplyr::rename(gene_mm = Mouse.gene.name) %>%
      dplyr::distinct(gene_mm, gene_hs, secretome_location, secretome_function) %>%
      dplyr::filter(!is.na(gene_mm))
  ),

  # Naba et al. 2012 (Mol Cell Proteomics) — mouse matrisome masterlist.
  # ~1110 genes across 6 categories in 2 divisions:
  #   Core matrisome:        ECM Glycoproteins, Collagens, Proteoglycans
  #   Matrisome-associated:  Secreted Factors, ECM Regulators, ECM-affiliated Proteins
  tar_target(naba_matrisome_mm_csv, "data/naba_2012_mcp/Mm_Matrisome_Masterlist.csv", format = "file"),
  tar_target(
    naba_matrisome_mm,
    readr::read_csv(naba_matrisome_mm_csv, skip = 1, show_col_types = FALSE) %>%
      dplyr::select(
        matrisome_division = `Matrisome Division`,
        matrisome_category = `Matrisome Category`,
        gene_mm            = `Gene Symbol`
      ) %>%
      dplyr::filter(!is.na(gene_mm))
  ),

  # ImmGen ATAC-seq OCR atlas (Yoshida et al. 2020 Cell) — annotation columns only
  tar_target(
    yoshida_2020_ocr_df,
    {
      df <- readr::read_csv(yoshida_2020_ocr_csv, col_select = 1:8,
                            show_col_types = FALSE)
      names(df)[1] <- "peakID"
      df <- dplyr::filter(df, `_-log10_bestPvalue` > 2)
      df
    }
  ),

  tar_target(
    yoshida_2020_ocr_gr,
    {
      gr <- GenomicRanges::GRanges(
        seqnames = yoshida_2020_ocr_df$chrom,
        ranges   = IRanges::IRanges(
          start = yoshida_2020_ocr_df$Summit - 250L,
          end   = yoshida_2020_ocr_df$Summit + 250L
        )
      )
      names(gr) <- paste0("atac:", yoshida_2020_ocr_df$peakID)
      gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
      GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gr)) <- "mm10"
      gr
    }
  ),

  tar_target(
    mm10_protein_coding_genes_gr,
    build_protein_coding_genes_gr(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79, species_genome = "mm10")
  ),

  tar_target(
    mm10_protein_coding_promoters_gr,
    build_protein_coding_promoters_gr(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                      mm10_protein_coding_genes_gr,
                                      species_genome = "mm10")
  ),

  tar_target(
    mm10_non_blood_promoters_gr,
    IRanges::subsetByOverlaps(
      mm10_protein_coding_promoters_gr,
      yoshida_2020_ocr_gr,
      invert = TRUE
    )
  ),

  tar_target(
    yoshida_2020_ocr_linked_orig,
    parse_yoshida_orig_links(yoshida_2020_ocr_gr, yoshida_2020_ocr_df)
  ),

  tar_target(
    yoshida_2020_ocr_linked_100kb,
    link_ranges_to_genes(yoshida_2020_ocr_gr, mm10_protein_coding_genes_gr, link_kb = 100)
  ),

  tar_target(
    yoshida_2020_ocr_linked_50kb,
    link_ranges_to_genes(yoshida_2020_ocr_gr, mm10_protein_coding_genes_gr, link_kb = 50)
  ),

  # Named gene lists — one element per cell type / lineage.
  tar_target(nestorowa_17_bl_mm_genesets,
             readr::read_csv(nestorowa_17_bl_mm_csv, show_col_types = FALSE) |>
               with(split(gene, lineage))),
  tar_target(kucinski_24_csc_mm_genesets,
             readr::read_csv(kucinski_24_csc_mm_csv, show_col_types = FALSE) |>
               with(split(marker_genes, cell_identity))),
  tar_target(isobe_23_cg_mm_genesets,
             readr::read_csv(isobe_23_cg_mm_csv, show_col_types = FALSE) |>
               with(split(gene, trajectory))),
  tar_target(hematlas_hspc_mm_genesets,
             readr::read_csv(hematlas_hspc_mm_csv, show_col_types = FALSE) |>
               with(split(gene_symbol, subcluster))),
  tar_target(hematlas_lineage_mm_genesets,
             readr::read_csv(hematlas_lineage_mm_csv, show_col_types = FALSE) |>
               with(split(gene_symbol, cell_type))),
  tar_target(ferchen_26_mm_marker_genesets,
             readr::read_csv(ferchen_26_mm_marker_csv, show_col_types = FALSE) |>
               with(split(marker_gene, level_4))),
  tar_target(ferchen_26_mm_tf_genesets,
             readr::read_csv(ferchen_26_mm_tf_csv, show_col_types = FALSE) |>
               with(split(TF, level_4))),
  tar_target(gavish_23_mm_marker_genesets,
             readr::read_csv(gavish_23_mm_marker_csv, show_col_types = FALSE) |>
               with(split(marker_gene, level_4))),
  tar_target(gavish_23_mm_tf_genesets,
             readr::read_csv(gavish_23_mm_tf_csv, show_col_types = FALSE) |>
               with(split(TF, level_4))),

  # Human → mouse ortholog-mapped genesets.
  # Each CSV is joined to mm_hs_ortho_long before splitting into gene lists.
  tar_target(psaila_20_hs_genesets,
             readr::read_csv(psaila_20_hs_csv, show_col_types = FALSE) |>
               dplyr::inner_join(mm_hs_ortho_long, by = c("gene" = "human_ortho")) |>
               with(split(Mouse.gene.name, signature))),

  tar_target(pan_25_hs_genesets,
             readr::read_csv(pan_25_hs_csv, show_col_types = FALSE) |>
               dplyr::inner_join(mm_hs_ortho_long, by = c("gene" = "human_ortho")) |>
               with(split(Mouse.gene.name, pattern))),

  tar_target(gavish_23_hs_scenic_genesets,
             readr::read_csv(gavish_23_hs_scenic_csv, show_col_types = FALSE) |>
               dplyr::mutate(program = paste0(meta_program, "_", regulon_type)) |>
               dplyr::inner_join(mm_hs_ortho_long, by = c("TF" = "human_ortho")) |>
               with(split(Mouse.gene.name, program))),

  tar_target(gavish_23_hs_meta_program_genesets,
             readr::read_csv(gavish_23_hs_meta_program_csv, show_col_types = FALSE) |>
               dplyr::inner_join(mm_hs_ortho_long, by = c("gene" = "human_ortho")) |>
               with(split(Mouse.gene.name, meta_program))),

  tar_target(gavish_23_hs_cancer_mp_genesets,
             readr::read_csv(gavish_23_hs_cancer_mp_csv, show_col_types = FALSE) |>
               dplyr::mutate(sheet   = stringr::str_remove(sheet, " MP.?"),
                             program = paste(sheet, program)) |>
               dplyr::inner_join(mm_hs_ortho_long, by = c("gene" = "human_ortho")) |>
               with(split(Mouse.gene.name, program))),

  # UCell-based module score targets — rank-based, no background bias.
  # Each target is self-contained: UMAP coords + sample.name included so
  # geneset-exploration.qmd does not need to load the Seurat object.
  tar_target(
    ucell_scores_nestorowa_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = nestorowa_17_bl_mm_genesets)
  ),

  tar_target(
    ucell_scores_isobe_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = isobe_23_cg_mm_genesets)
  ),

  tar_target(
    ucell_scores_psaila_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = psaila_20_hs_genesets)
  ),

  tar_target(
    ucell_scores_pan_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = pan_25_hs_genesets)
  ),

  tar_target(
    ucell_scores_gavish_meta_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = gavish_23_hs_meta_program_genesets)
  ),

  tar_target(
    ucell_scores_gavish_cancer_hsc,
    extract_ucell_scores(mpn_processed_hsc, genesets = gavish_23_hs_cancer_mp_genesets)
  ),

  # ── RNA velocity: h5ad export ────────────────────────────────────────────────
  # Exports the processed HSC Seurat object to AnnData h5ad for use in the
  # Python velocity-hsc.ipynb notebook (scVelo + CellRank).
  # Regenerated automatically if mpn_processed_data_hsc changes.
  tar_target(
    hsc_h5ad_path,
    {
      srt <- mpn_processed_hsc_fullcc
      # Add CytoTRACE2 scores to cell metadata for CytoTRACEKernel / PseudotimeKernel
      idx <- match(colnames(srt), rownames(cytotrace_hsc))
      srt$cytotrace2_score    <- cytotrace_hsc$CytoTRACE2_Score[idx]
      srt$cytotrace2_potency  <- cytotrace_hsc$CytoTRACE2_Potency[idx]
      srt$cytotrace2_relative <- cytotrace_hsc$CytoTRACE2_Relative[idx]
      export_seurat_h5ad(
        srt  = srt,
        path = file.path(!!config$paths$objects_dir, "hsc.h5ad")
      )
    },
    format = "file"
  ),

  # ── MAGIC imputation — HSC ────────────────────────────────────────────────────
  # Runs MAGIC on the SCT normalised data of the full cell-cycle HSC object.
  # Returns a list: $raw and $magic are cells × HVG matrices for before/after.
  # Rmagic runs the Python package magic-impute from the project's .venv (RETICULATE_PYTHON, top of file).
  tar_target(
    magic_hsc,
    {
      genes <- Seurat::VariableFeatures(mpn_processed_hsc_fullcc)
      compute_magic(mpn_processed_hsc_fullcc, genes = genes, t = 3, knn = 5)
    }
  ),

  # ── SCENIC: per-condition AUC + aggregation (scripts/slurm/analysis/scenic.sh runs 4 subsets) ────
  # Suffixes: full | mch_neg | mch_pos | wt_ctrl
  # To switch cluster resolution: edit scenic_agg_column.
  tar_target(scenic_agg_column, "sNN_pca.SCT_30_res.1.5"),
  tar_target(scenic_meta_hsc,   scenic_build_meta(anno_meta_hsc, scenic_agg_column)),

  tarchetypes::tar_map(
    values = tibble::tibble(
      suffix   = c("full", "mch_neg", "mch_pos", "wt_ctrl"),
      auc_file = paste0("output/objects/scenic/auc_mtx_clean_",       suffix, ".csv"),
      reg_file = paste0("output/objects/scenic/regulons_tidy_clean_", suffix, ".csv")
    ),
    names = suffix,
    tar_target(scenic_auc_path_hsc,    auc_file, format = "file"),
    tar_target(scenic_reg_path_hsc,    reg_file, format = "file"),
    tar_target(scenic_auc_mat_hsc,     scenic_load_auc(scenic_auc_path_hsc)),
    tar_target(scenic_reg_tidy_hsc,    readr::read_csv(scenic_reg_path_hsc, show_col_types = FALSE)),
    tar_target(scenic_agg_hsc,         scenic_aggregate_auc(scenic_auc_mat_hsc, scenic_meta_hsc, scenic_agg_column)),
    tar_target(scenic_prop_mat_hsc,    scenic_cond_prop_mat(scenic_meta_hsc, scenic_agg_column, rownames(scenic_agg_hsc))),
    tar_target(scenic_hclust_rows_hsc, hclust(dist(t(scale(scenic_agg_hsc))), method = "ward.D2")),
    tar_target(scenic_hclust_cols_hsc, hclust(dist(scale(scenic_agg_hsc)),    method = "ward.D2")),
    tar_target(scenic_cutree_rows_hsc, cutree(scenic_hclust_rows_hsc, k = 8))
  ),

  # Transposed full-cohort AUC as a Seurat Assay5 (regulons × cells).
  # Used in scenic-hsc.qmd section 5 to run PCA/UMAP on regulon activity.
  tar_target(
    scenic_assay_hsc,
    {
      common_cells <- intersect(rownames(scenic_auc_mat_hsc_full), colnames(mpn_processed_hsc_fullcc))
      auc_t <- t(scenic_auc_mat_hsc_full[common_cells, ])   # regulons × cells
      SeuratObject::CreateAssay5Object(data = auc_t)
    }
  ),

  # ── Alternative dim-reductions: PHATE + FA2 ──────────────────────────────────
  # Returns embedding matrices (cells × 2) only — the full Seurat object is not
  # cached.  Notebooks add them back via CreateDimReducObject().
  tarchetypes::tar_map(
    values = tibble::tibble(
      compartment = c("stroma",               "hsc",               "mnc"),
      srt_target = rlang::syms(c("mpn_processed_data_stroma", "mpn_processed_hsc_fullcc", "mpn_processed_data_mnc"))
    ),
    names = compartment,

    tar_target(
      phate_coords,
      compute_phate_coords(srt_target,
                           reduction = "pca.SCT",
                           dims      = 1:!!config$reduction$n_pca_dims,
                           knn       = 5, decay = 40, t = "auto")
    ),

    tar_target(
      fa2_coords,
      compute_fa2_coords(srt_target, graph = NULL, niter = 500)
    )
  ),

  # ── Cluster-level expression summaries — HSC res 1.5 ─────────────────────────
  # Average SCT expression: mean of log-normalised SCT data layer per cluster.
  tar_target(
    cluster_avg_expr_hsc_res1_5,
    {
      srt <- mpn_processed_hsc_fullcc
      mat <- Seurat::AverageExpression(
        srt,
        assays        = "SCT",
        group.by      = "sNN_pca.SCT_30_res.1.5",
        layer         = "data",
        return.seurat = FALSE
      )[["SCT"]]
      df <- tibble::rownames_to_column(as.data.frame(as.matrix(mat)), "gene")
      rm(srt, mat); gc()
      df
    }
  ),

  tar_target(
    cluster_avg_expr_hsc_res1_5_path,
    save_csv(cluster_avg_expr_hsc_res1_5,
             !!config$paths$output_table_dir,
             "cluster_avg_expr_hsc_res1_5.csv"),
    format = "file"
  ),

  # Pseudobulk raw counts: sum of RNA counts per cluster.
  tar_target(
    cluster_pseudobulk_hsc_res1_5,
    {
      srt <- mpn_processed_hsc_fullcc
      mat <- Seurat::AggregateExpression(
        srt,
        assays        = "RNA",
        group.by      = "sNN_pca.SCT_30_res.1.5",
        layer         = "counts",
        return.seurat = FALSE
      )[["RNA"]]
      df <- tibble::rownames_to_column(as.data.frame(as.matrix(mat)), "gene")
      rm(srt, mat); gc()
      df
    }
  ),

  tar_target(
    cluster_pseudobulk_hsc_res1_5_path,
    save_csv(cluster_pseudobulk_hsc_res1_5,
             !!config$paths$output_table_dir,
             "cluster_pseudobulk_counts_hsc_res1_5.csv"),
    format = "file"
  ),

  # Pseudobulk TMM-normalized log-CPM (edgeR).
  tar_target(
    cluster_pseudobulk_norm_hsc_res1_5,
    {
      counts_mat <- as.matrix(cluster_pseudobulk_hsc_res1_5[, -1])
      rownames(counts_mat) <- cluster_pseudobulk_hsc_res1_5$gene
      dge      <- edgeR::DGEList(counts = counts_mat)
      dge      <- edgeR::calcNormFactors(dge, method = "TMM")
      norm_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
      tibble::rownames_to_column(as.data.frame(norm_mat), "gene")
    }
  ),

  tar_target(
    cluster_pseudobulk_norm_hsc_res1_5_path,
    save_csv(cluster_pseudobulk_norm_hsc_res1_5,
             !!config$paths$output_table_dir,
             "cluster_pseudobulk_norm_hsc_res1_5.csv"),
    format = "file"
  ),

  # ── Cross-compartment merged objects (general-purpose) ───────────────────────
  tar_target(mpn_all,         SeuratObject:::merge.Seurat(mpn_processed_hsc_fullcc, list(mpn_processed_mnc, mpn_processed_stroma), merge.data = FALSE)),
  tar_target(mpn_wt_all,      SeuratObject:::subset.Seurat(mpn_all, subset = grepl("_WT_", sample.name, fixed = TRUE))),
  tar_target(mpn_mch_neg_all, SeuratObject:::subset.Seurat(mpn_all, subset = grepl("mCh-", sample.name, fixed = TRUE))),
  tar_target(mpn_mch_pos_all, SeuratObject:::subset.Seurat(mpn_all, subset = grepl("mCh+", sample.name, fixed = TRUE) | (compartment == "Stroma" & grepl("Het", sample.name, fixed = TRUE)))),

  # ── Shared colour palettes ────────────────────────────────────────────────────
  tar_target(
    gradientn_pals,
    list(
      blue    = c("grey88", "#14B3FF", "#2E8BE0", "#2A5DB0", "#10316E"),
      red     = c("grey88", "#FF6B6B", "#F03B3B", "#C81E1E", "#7A0C0C"),
      heat    = c("grey88", "#FDB31A", "#FB6A4A", "#E42A2A", "#A31D1D", "#67000D"),
      green   = c("grey88", "#4FD16B", "#22A53C", "#157A28", "#0A4715"),
      orange  = c("grey88", "#FFA63D", "#FB7A1E", "#E25410", "#8C2D04"),
      purple  = c("grey88", "#B26BFF", "#8A3FE0", "#6A24B0", "#3D106E"),
      teal    = c("grey88", "#2EE0C0", "#14B8A6", "#0E8074", "#0A4A44"),
      magenta = c("grey88", "#FF5CC8", "#E0289E", "#A8156E", "#5E0A3E")
    )
  ),

  # ── Reports ─────────────────────────────────────────────────────────────────
  # One folder per report under notebooks/, rendered to output/reports/latest/ (quarto_report(),
  # above the list). Each lists the files it reads outside targets, so it re-renders when they change.

  quarto_report(report_demultiplex_diagnostics, "demultiplex-diagnostics", utilities_r),

  # TEST (temporary): delete these two targets, notebooks/test_report_target/ and
  # output/reports/latest/notebooks/test_report_target/ afterwards
  tar_target(test_report_data, {
    set.seed(1)
    data.frame(group = rep(c("A", "B", "C"), each = 50),
               value = c(rnorm(50, 0), rnorm(50, 1), rnorm(50, 2.5)))
  }),
  quarto_report(report_test_report_target, "test_report_target")   # last element of the list: no comma

  # QC and demultiplexing
  # quarto_report(report_first_analysis_milestone_report, "first-analysis-milestone-report", utilities_r),
  # quarto_report(report_qc_exploration, "qc-exploration", utilities_r, config_yaml)
  # quarto_report(report_transgene_detection, "transgene-detection", utilities_r, config_yaml),
  # quarto_report(report_demultiplex_before_qc, "demultiplex_before_qc", utilities_r),
  # quarto_report(report_demultiplex_after_qc, "demultiplex_after_qc", utilities_r),
  # quarto_report(report_gene_biotypes, "gene-biotypes", utilities_r),

  # Dimensionality reduction, clustering, cell cycle
  # quarto_report(report_dim_reduction_and_clustering_RNA, "dim-reduction-and-clustering-RNA", utilities_r, config_yaml),
  # quarto_report(report_dim_reduction_and_clustering_SCT, "dim-reduction-and-clustering-SCT"),
  # quarto_report(report_pca_cell_distribution, "pca-cell-distribution", utilities_r),
  # quarto_report(report_alt_dim_reduction, "alt-dim-reduction", utilities_r),
  # quarto_report(report_alt_dim_reduction_hsc, "alt-dim-reduction-hsc", utilities_r),
  # quarto_report(report_cell_cycle, "cell-cycle", utilities_r),
  # quarto_report(report_cell_cycle_hsc, "cell-cycle-hsc", utilities_r),
  # quarto_report(report_hsc_per_condition_embedding, "hsc-per-condition-embedding", utilities_r),
  # quarto_report(report_subcluster_exploration, "subcluster-exploration", utilities_r, config_yaml),
  # quarto_report(report_supercell_exploration, "supercell-exploration"),
  # quarto_report(report_graph_analysis, "graph-analysis"),
  # quarto_report(report_magic_denoising_hsc, "magic-denoising-hsc"),
  # quarto_report(report_mpig6b_expression, "mpig6b-expression", utilities_r, config_yaml),

  # Annotation and gene sets
  # quarto_report(report_cluster_annotation, "cluster-annotation", utilities_r, config_yaml),
  # quarto_report(report_reference_annotation, "reference-annotation", utilities_r),
  # quarto_report(report_reference_annotation_hsc, "reference-annotation-hsc", utilities_r),
  # quarto_report(report_reference_annotation_mnc, "reference-annotation-mnc", utilities_r),
  # quarto_report(report_reference_annotation_stroma, "reference-annotation-stroma", utilities_r),
  # quarto_report(report_geneset_exploration, "geneset-exploration"),

  # Differential abundance (Milo)
  # quarto_report(report_milo_hsc, "milo-hsc", utilities_r),
  # quarto_report(report_milo_mnc, "milo-mnc", utilities_r),
  # quarto_report(report_milo_stroma, "milo-stroma", utilities_r),
  # quarto_report(report_miloR_to_genes, "miloR_to_genes", utilities_r),

  # Gene regulation and cell communication
  # quarto_report(report_hdwgcna_hsc, "hdwgcna-hsc"),
  # quarto_report(report_hdwgcna_hsc_modules, "hdwgcna-hsc-modules"),
  # quarto_report(report_hdwgcna_hsc_grn_vanilla, "hdwgcna-hsc-grn-vanilla"),
  # quarto_report(report_hdwgcna_hsc_grn_promoters, "hdwgcna-hsc-grn-promoters"),
  # quarto_report(report_hdwgcna_hsc_grn_non_blood_promoters, "hdwgcna-hsc-grn-non-blood-promoters"),
  # quarto_report(report_hdwgcna_hsc_grn_distal, "hdwgcna-hsc-grn-distal"),
  # quarto_report(report_hdwgcna_hsc_grn_atac, "hdwgcna-hsc-grn-atac"),
  # quarto_report(report_hdwgcna_hsc_grn_both, "hdwgcna-hsc-grn-both"),
  # quarto_report(report_scenic_hsc, "scenic-hsc", utilities_r),
  # quarto_report(report_scenic_hsc_conditions, "scenic-hsc-conditions", utilities_r),
  # quarto_report(report_upstream_regulators_hsc, "upstream-regulators-hsc"),
  # quarto_report(report_cellchat_wt, "cellchat-wt"),
  # quarto_report(report_cellchat_mch_neg, "cellchat-mch-neg"),
  # quarto_report(report_cellchat_mch_pos, "cellchat-mch-pos"),

  # Potency and trajectories (paga-*, dpt-*, cellrank-*: outputs of the cluster jobs in scripts/slurm/)
  # quarto_report(report_differentiation_potential, "differentiation_potential", utilities_r),
  # quarto_report(report_trajectory_slingshot, "trajectory-slingshot"),
  # quarto_report(report_trajectory_slingshot_fullcc_res09, "trajectory-slingshot-fullcc-res09", utilities_r),
  # quarto_report(report_trajectory_slingshot_fullcc_res09_distsl, "trajectory-slingshot-fullcc-res09-distsl", utilities_r),
  # quarto_report(report_trajectory_slingshot_fullcc_res09_sub, "trajectory-slingshot-fullcc-res09_sub", utilities_r),
  # quarto_report(report_trajectory_slingshot_fullcc_res09_sub_subset, "trajectory-slingshot-fullcc-res09_sub_subset", utilities_r),
  # quarto_report(report_trajectory_slingshot_fullcc_res15, "trajectory-slingshot-fullcc-res15", utilities_r),
  # quarto_report(report_trajectory_slingshot_fullcc_res15_distsl, "trajectory-slingshot-fullcc-res15-distsl", utilities_r),
  # quarto_report(report_trajectory_slingshot_hsc_c6c8_to_c4_res09, "trajectory-slingshot-hsc-c6c8-to-c4-res09", utilities_r),
  # quarto_report(report_trajectory_slingshot_mchneg_res09, "trajectory-slingshot-mchneg-res09", utilities_r),
  # quarto_report(report_trajectory_slingshot_mchneg_res15, "trajectory-slingshot-mchneg-res15", utilities_r),
  # quarto_report(report_trajectory_slingshot_mchpos_res09, "trajectory-slingshot-mchpos-res09", utilities_r),
  # quarto_report(report_trajectory_slingshot_mchpos_res15, "trajectory-slingshot-mchpos-res15", utilities_r),
  # quarto_report(report_transition_matrices, "transition-matrices"),
  # quarto_report(report_paga_hsc, "paga-hsc"),
  # quarto_report(report_paga_traj_analysis_0_9, "paga-traj-analysis-0.9", "output/objects/paga/paga_edges_sNN_pca.SCT_30_res.0.9.csv", "output/objects/paga_conditions/paga_edge_comparison_sNN_pca.SCT_30_res.0.9.csv"),
  # quarto_report(report_paga_traj_analysis_1_0, "paga-traj-analysis-1.0", "output/objects/paga/paga_edges_sNN_pca.SCT_30_res.1.csv", "output/objects/paga_conditions/paga_edge_comparison_sNN_pca.SCT_30_res.1.csv"),
  # quarto_report(report_paga_traj_analysis_1_5, "paga-traj-analysis-1.5", "output/objects/paga/paga_edges_sNN_pca.SCT_30_res.1.5.csv", "output/objects/paga_conditions/paga_edge_comparison_sNN_pca.SCT_30_res.1.5.csv"),
  # quarto_report(report_dpt_traj_analysis, "dpt-traj-analysis", utilities_r, "output/objects/dpt/dpt_cell_table.csv"),
  # quarto_report(report_cellrank_hsc, "cellrank-hsc", "output/objects/cellrank/cytotrace"),
  # quarto_report(report_cellrank_hsc_cytotrace, "cellrank-hsc-cytotrace", "output/objects/cellrank/cytotrace_res1_5"),
  # quarto_report(report_cellrank_hsc_dpt, "cellrank-hsc-dpt", "output/objects/cellrank/dpt_res1_5")
)
