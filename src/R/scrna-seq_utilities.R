# scrna-seq_utilities.R
# scRNA-seq helper functions for the MPN 10x project.
# Adapted from Shoeb2025_CG-AMKL (https://github.com/cancerbits/Shoeb2025_CG-AMKL)

# ---------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------

#' Load a single 10x CellRanger sample → named Seurat object
#'
#' @param path       The sample's filtered matrix folder (Cell Ranger 10 multi:
#'                   per_sample_outs/<id>/sample_filtered_feature_bc_matrix).
#' @param sample_id  Single sample ID (the tube, e.g. "Rx_1").
#' @param transgene_features  Features of the transduced vector (metadata/project_config.yaml).
#'                   They are moved out of RNA into a Transgene assay, so that they cannot reach
#'                   variable features, PCA or clustering, and nCount_RNA/nFeature_RNA measure the
#'                   transcriptome only. Keep them as counts: three features carry no size factors
#'                   of their own, so normalise against nCount_RNA if you need a rate.
#' @return Seurat object with cell.name, sample.name, seq_batch metadata.
read_10x_sample <- function(path, sample_id, sample_name_short, sample_name_long,
                            transgene_features = character(0)) {
  if (!dir.exists(path))
    stop("CellRanger directory not found: ", path)
  mat <- Seurat::Read10X(path)
  # Read full feature table (col1=ensembl_id, col2=gene_name, col3=feature_type)
  feat <- read.table(file.path(path, "features.tsv.gz"),
                     sep = "\t", header = FALSE,
                     col.names = c("ensembl_id", "gene_name", "feature_type"))
  rownames(feat) <- make.unique(feat$gene_name)
  # Gene Expression assay, with the vector's features held back for their own assay
  gex <- if (is.list(mat)) mat[["Gene Expression"]] else mat
  tg_rows <- intersect(transgene_features, rownames(gex))
  missing_tg <- setdiff(transgene_features, tg_rows)
  if (length(missing_tg))
    warning(sample_id, ": transgene feature(s) not in the reference: ",
            paste(missing_tg, collapse = ", "))
  srt <- Seurat::CreateSeuratObject(counts = gex[setdiff(rownames(gex), tg_rows), ],
                                    project = sample_id,
                                    min.cells = 0, min.features = 0)
  if (length(tg_rows))
    srt[["Transgene"]] <- Seurat::CreateAssayObject(counts = gex[tg_rows, , drop = FALSE])
  srt[["cell.name"]]   <- colnames(srt)
  srt[["tube"]]        <- sample_id
  srt[["sample.name"]] <- sample_name_short[srt$tube] %>% unname()
  srt[["sample.name.long"]] <- sample_name_long[srt$tube] %>% unname()

  # Compartment derived from sample name: Stroma | HSC | MNC
  srt[["compartment"]] <- dplyr::case_when(
    grepl("Stroma",  srt[["sample.name"]], ignore.case = TRUE) ~ "Stroma",
    grepl("HSC",     srt[["sample.name"]], ignore.case = TRUE) ~ "HSC",
    grepl("MNC",     srt[["sample.name"]], ignore.case = TRUE) ~ "MNC",
    TRUE ~ NA_character_
  )

  srt[["RNA"]]@meta.data[["ensembl_id"]] <- feat[rownames(srt[["RNA"]]), "ensembl_id"]
  # HTO assay (Antibody Capture) - rename count matrix rows to antibody names (e.g. TotalSeq_B0301) instead of sample-specific names (e.g. WT_GFP_NC)
  # Remove sample information to report on antibody level (Antonio's idea)
  sample_ab_id <- rownames(mat[["Antibody Capture"]])
  rownames(mat[["Antibody Capture"]]) <- feat[sample_ab_id, "ensembl_id"]
  srt[["HTO"]] <- Seurat::CreateAssayObject(counts = mat[["Antibody Capture"]])
  ab_meta_features <- feat[sample_ab_id, ] %>% magrittr::set_rownames(stringr::str_replace(.$ensembl_id, "_", "-"))
  srt[["HTO"]]@meta.features <- ab_meta_features

  srt
}

#' Merge a named list of per-sample Seurat objects and join layers (Seurat v5)
#'
#' @param srt_lst  Named list of Seurat objects (names = sample IDs).
#' @return Merged Seurat object with joined layers.
merge_seurat_list <- function(srt_lst) {
  sample_ids <- names(srt_lst)
  merged <- merge(srt_lst[[1]],
                  y            = srt_lst[-1],
                  add.cell.ids = sample_ids,  # e.g. Rx_1-ACGTACGTACGT-1
                  project      = "MPN_10x")
  rm(srt_lst); gc()
  SeuratObject::JoinLayers(merged)
}

# ---------------------------------------------------------------
# Sub-clustering
# ---------------------------------------------------------------

#' Subset a processed Seurat object to specific clusters at a given resolution.
#'
#' Identifies the cluster column matching `resolution` in `@meta.data`, then
#' keeps only cells whose cluster identity is in `clusters`. The original
#' cluster column is preserved as `parent_cluster` so provenance is traceable.
#'
#' @param srt        Seurat object (post-`process_all`); must contain the
#'                   cluster column for `resolution`.
#' @param clusters   Character vector of cluster IDs to keep (e.g. c("5", "9")).
#' @param resolution Numeric resolution used to select the cluster column
#'                   (matched via `sNN_pca.*_res.<resolution>` grep).
#'
#' @return Seurat object containing only the selected cells, with a
#'         `parent_cluster` metadata column added.
subset_clusters <- function(srt, clusters, resolution) {
  clust_col <- grep(
    paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", resolution), "$"),
    colnames(srt@meta.data), value = TRUE
  )[1]
  if (is.na(clust_col))
    stop("No cluster column found for resolution ", resolution,
         ". Available columns: ",
         paste(grep("^sNN_", colnames(srt@meta.data), value = TRUE), collapse = ", "))

  keep <- srt@meta.data[[clust_col]] %in% as.character(clusters)
  message(sprintf("subset_clusters: keeping %d / %d cells (clusters %s at res %s)",
                  sum(keep), ncol(srt), paste(clusters, collapse = ","), resolution))

  srt_sub <- srt[, keep]
  srt_sub$parent_cluster <- srt_sub@meta.data[[clust_col]]
  srt_sub
}

# ---------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------

#' deMULTIplex2 per tube (orig.ident).
#'
#' plot.umap = "umi" is needed for the UMAP coordinates in the result, but makes demultiplexTags()
#' always write a PDF (to plot.path, default getwd()); it goes to a temporary folder deleted on exit.
#' The plots are redrawn from the stored result with dm2_summary_plot() / dm2_diagnostic_plots().
#' @return Named list (one per tube) of demultiplexTags() results.
demux_per_rx <- function(cg_srt){
  plot_dir <- tempfile("demultiplex2-")
  dir.create(plot_dir)
  on.exit(unlink(plot_dir, recursive = TRUE), add = TRUE)

  lst <- Seurat::SplitObject(cg_srt, split.by = "orig.ident")
  lst_dm2 <- lapply(lst, function(rx_srt){
    hto_counts <- Seurat::GetAssayData(rx_srt, assay = "HTO", layer = "counts")
    dm2_res    <- deMULTIplex2::demultiplexTags(
      tag_mtx          = BiocGenerics::t(hto_counts),
      plot.umap        = "umi",
      plot.diagnostics = FALSE,
      plot.path        = plot_dir
    )
    dm2_res
  })

  return(lst_dm2)
}

# Data frame deMULTIplex2 plots from: UMAP coordinates + assignment table, barcode_count binned
.dm2_umap_df <- function(dm2_res) {
  umap_df <- cbind(as.data.frame(dm2_res$umap), dm2_res$assign_table)
  umap_df$barcode_count <- as.character(umap_df$barcode_count)
  umap_df$barcode_count[as.numeric(umap_df$barcode_count) >= 3] <- ">=3"
  umap_df$barcode_count <- factor(umap_df$barcode_count, c("0", "1", "2", ">=3"))
  umap_df
}

#' deMULTIplex2's summary page (UMAP of the tag counts by assignment), from a stored result.
#' Same as the first page of demultiplexTags()'s PDF. Returns a grob: draw with grid::grid.draw().
dm2_summary_plot <- function(dm2_res, point.size = 1, label.size = 3, min.tag.show = 50) {
  deMULTIplex2:::plotSummary(.dm2_umap_df(dm2_res), point.size = point.size,
                             label.size = label.size, min.tag.show = min.tag.show)
}

#' deMULTIplex2's per-tag diagnostic pages (plot.diagnostics = TRUE), from a stored result.
#' prob.cut must match the one used for the assignment (demultiplexTags() default 0.5).
#' Returns a named list of grobs, one per tag: draw with grid::grid.draw().
dm2_diagnostic_plots <- function(dm2_res, prob.cut = 0.5, point.size = 1) {
  umap_df  <- .dm2_umap_df(dm2_res)
  mappings <- list(c("log(bc.umi)", "cos.umi", "cos.umi"),
                   c("log(tt.umi)", "log(bc.umi)", "prob.pos"),
                   c("log(tt.umi)", "log(tt.umi-bc.umi)", "prob.pos"),
                   c("log(tt.umi)", "res", "prob.pos"),
                   c("UMAP_1", "UMAP_2", "prob.pos"),
                   c("UMAP_1", "UMAP_2", "barcode_count"))
  tags <- names(dm2_res$df_list)
  stats::setNames(lapply(tags, function(bc) {
    fit <- dm2_res$df_list[[bc]]
    df  <- data.frame(bc.umi   = fit$bc.umi,
                      res      = as.numeric(dm2_res$res_mtx[, bc]),
                      cos.umi  = fit$cos.umi,
                      prob.pos = dm2_res$prob_mtx[, bc],
                      tt.umi   = fit$tt.umi)
    df <- cbind(df, umap_df)
    df$barcode_count[!(dm2_res$prob_mtx[, bc] > prob.cut)] <- NA   # cells not called for this tag
    deMULTIplex2:::plot.all.diagnostics(df, mappings, bc = bc, prob.cut = prob.cut,
                                        point.size = point.size, ncol = 3)
  }), tags)
}

#' Compute QC metrics, HTO demultiplexing, per-metric boolean flags, and
#' doublet scores — attaching all results to the Seurat object's metadata.
#'
#' Steps:
#'   1. QC metrics (percent.mt/rp/hb, count_feat_ratio, percent.stress)
#'   2. Cell-cycle scoring
#'   3. HTO demultiplexing (dm2_class)
#'   4. Per-metric boolean flags (nFeature_min and ribo_percent_min are applied
#'      per-compartment using the compartment column in cg_srt@meta.data):
#'        pass_mt       — percent.mt   <= mito_percent_max
#'        pass_nFeature — nFeature_RNA in [nFeature_min, nFeature_max]
#'        pass_rp       — percent.rp   in [ribo_percent_min, ribo_percent_max]
#'        pass_hb       — percent.hb   <= hb_percent_max
#'        pass_qc       — all four metric flags combined
#'        pass_hto      — dm2_class not in remove_hto_classes
#'        pass_cell     — pass_qc & pass_hto
#'   5. scDblFinder on pass_qc cells; scores written to all cells (NA otherwise)
#'
#' @param compartment_qc Named list keyed by compartment label (e.g. "Stroma",
#'   "HSC", "MNC"), each with fields nFeature_min, nFeature_max, ribo_percent_min, and ribo_percent_max.
#' @return Seurat object with QC metadata columns added.
flag_qc_metadata <- function(cg_srt,
                             min_cells_per_gene,
                             min_counts_per_cell_per_gene,
                             mito_percent_max,
                             hb_percent_max,
                             mt_genes,
                             rp_genes,
                             hb_genes,
                             stress_genes,
                             compartment_qc,
                             demux_out) {
  # ── 1. QC metrics ───────────────────────────────────────────────
  cg_srt[["percent.mt"]]     <- Seurat::PercentageFeatureSet(cg_srt, features = mt_genes)
  cg_srt[["percent.rp"]]     <- Seurat::PercentageFeatureSet(cg_srt, features = rp_genes)
  cg_srt[["percent.hb"]]     <- Seurat::PercentageFeatureSet(cg_srt, features = hb_genes)
  cg_srt[["count_feat_ratio"]] <- cg_srt@meta.data$nCount_RNA / cg_srt@meta.data$nFeature_RNA
  cg_srt[["percent.stress"]] <- Seurat::PercentageFeatureSet(cg_srt, features = stress_genes)

  # ── 2. Cell-cycle scoring ────────────────────────────────────────
  cg_srt <- Seurat::NormalizeData(cg_srt, verbose = FALSE)
  cg_srt <- Seurat::CellCycleScoring(
    cg_srt,
    g2m.features = stringr::str_to_title(cc.genes$g2m.genes),
    s.features   = stringr::str_to_title(cc.genes$s.genes)
  )
  cg_srt[["RNA"]]$data <- NULL

  # ── 3. HTO demultiplexing ────────────────────────────────────────
  max_prop <- lapply(demux_out,`[[`, "prob_mtx") %>% do.call(rbind,.) %>% apply(1,max)
  max_prop <- max_prop[colnames(cg_srt)]
  names(max_prop) <- colnames(cg_srt)
  max_prop <- ifelse(is.na(max_prop),0,max_prop)
  cg_srt@meta.data$dm2_max_prop <-  max_prop

  max_resid <- lapply(demux_out,`[[`, "res_mtx") %>% do.call(rbind,.) %>% apply(1,max)
  max_resid <- max_resid[colnames(cg_srt)]
  names(max_resid) <- colnames(cg_srt)
  max_resid <- ifelse(is.na(max_resid),0,max_resid)
  cg_srt@meta.data$dm2_max_resid <-  max_resid[colnames(cg_srt)]

  assign_vec <- lapply(demux_out,`[[`, "final_assign") %>% unlist
  names(assign_vec) <- str_remove(names(assign_vec),"^Rx_..")
  assign_vec <- assign_vec[colnames(cg_srt)]
  names(assign_vec) <- colnames(cg_srt)
  assign_vec <- ifelse(is.na(assign_vec), "negative", assign_vec)
  cg_srt@meta.data$dm2_class <-  assign_vec

  # ── 4. Per-metric boolean flags ──────────────────────────────────
  md <- cg_srt@meta.data

  # nFeature_min/max, ribo_percent_min, and ribo_percent_max are looked up per compartment
  comp         <- md$compartment
  nFeature_min <- vapply(comp, function(c) compartment_qc[[c]]$nFeature_min,     numeric(1))
  nFeature_max <- vapply(comp, function(c) compartment_qc[[c]]$nFeature_max,     numeric(1))
  ribo_pct_min <- vapply(comp, function(c) compartment_qc[[c]]$ribo_percent_min, numeric(1))
  ribo_pct_max <- vapply(comp, function(c) compartment_qc[[c]]$ribo_percent_max, numeric(1))

  cg_srt@meta.data$pass_mt       <- md$percent.mt   <= mito_percent_max
  cg_srt@meta.data$pass_nFeature <- md$nFeature_RNA >= nFeature_min &
                                     md$nFeature_RNA <= nFeature_max
  cg_srt@meta.data$pass_rp       <- md$percent.rp   >= ribo_pct_min &
                                     md$percent.rp   <= ribo_pct_max
  cg_srt@meta.data$pass_hb       <- md$percent.hb   <= hb_percent_max
  cg_srt@meta.data$pass_qc       <- cg_srt@meta.data$pass_mt       &
                                     cg_srt@meta.data$pass_nFeature &
                                     cg_srt@meta.data$pass_rp       &
                                     cg_srt@meta.data$pass_hb
  cg_srt@meta.data$pass_hto      <- md$dm2_class != "multiplet"

  # ── 5. Doublet scoring on pass_qc cells × pass_gene features ────
  cell_mask  <- cg_srt@meta.data$pass_qc & cg_srt@meta.data$pass_hto
  counts_all <- Seurat::GetAssayData(cg_srt, assay = "RNA", layer = "counts")
  pass_gene  <- rownames(counts_all)[rowSums(counts_all[, cell_mask] > min_counts_per_cell_per_gene) >= min_cells_per_gene]
  counts_sub <- counts_all[pass_gene, cell_mask, drop = FALSE]
  rm(counts_all); gc()

  set.seed(1)
  sce        <- scDblFinder::scDblFinder(counts_sub,
                                         samples = cg_srt$orig.ident[cell_mask])
  score_vec  <- rep(NA_real_,   ncol(cg_srt))
  class_vec  <- rep("unscored", ncol(cg_srt))
  score_vec[cell_mask] <- sce$scDblFinder.score
  class_vec[cell_mask] <- as.character(sce$scDblFinder.class)
  cg_srt@meta.data$scDblFinder.score <- score_vec
  cg_srt@meta.data$scDblFinder.class <- class_vec
  # pass_dbl: singlet or unscored (cells that failed pass_qc were not scored)
  cg_srt@meta.data$pass_dbl  <- class_vec %in% c("singlet", "unscored")
  cg_srt@meta.data$pass_cell <- cg_srt@meta.data$pass_qc & cg_srt@meta.data$pass_hto #& cg_srt@meta.data$pass_dbl

  return(cg_srt)
}

#' Query Ensembl BioMart for gene biotypes and return an annotated data.frame.
#'
#' @param cg_srt Seurat object; gene symbols taken from rownames.
#' @return data.frame with columns: mgi_symbol, ensembl_gene_id, gene_biotype,
#'   chromosome_name. Genes absent from BioMart have gene_biotype = "not_found".
fetch_gene_biotypes <- function(cg_srt) {
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  genes <- rownames(cg_srt)
  mart  <- biomaRt::useEnsembl(biomart = "genes",
                               dataset = "mmusculus_gene_ensembl",
                               host = "https://may2025.archive.ensembl.org")
  gene_attributes <- c("mgi_symbol", "external_gene_name", "external_synonym",
                       "ensembl_gene_id", "gene_biotype", "chromosome_name")
  gene_filters <- c("mgi_symbol", "external_gene_name", "external_synonym")
  l <- lapply(gene_filters, function(filter_x){
    bm <- biomaRt::getBM(
      attributes = gene_attributes,
      filters    = filter_x,
      values     = genes,
      mart       = mart
    )
    bm$query_symbol <- bm[[filter_x]]
    bm
  })
  names(l) <- gene_filters
  dplyr::bind_rows(l, .id = "gene_filter")
  }

# ---------------------------------------------------------------
# Normalisation & integration
# ---------------------------------------------------------------

#' SCTransform normalisation + PCA
#'
#' SCTransform picks the variable features itself (ranked by residual variance); the excluded genes
#' are then removed from them, so the PCA runs on exactly `n_features` genes.
#'
#' @param cg_srt          Seurat object.
#' @param ref_assay       Must be "SCT".
#' @param n_features      Number of variable features to end up with (after exclusion).
#' @param exclude_genes   Exact gene names to exclude from variable features.
#' @param vars_to_regress Passed to SCTransform's vars.to.regress.
preprocess_seurat <- function(cg_srt, ref_assay, n_features, exclude_genes = character(0),
                              vars_to_regress = NULL) {
  # SCTransform uses future internally; set sequential plan and disable the
  # globals size check (the check fires even in sequential mode because SCT
  # captures the full Seurat object in its closure).
  future::plan("sequential")
  options(future.globals.maxSize = Inf)

  # Ask for extra features so that n_features remain after the excluded genes are removed
  cg_srt <- Seurat::SCTransform(cg_srt,
                                assay               = "RNA",
                                do.correct.umi      = FALSE,
                                variable.features.n = n_features + length(exclude_genes),
                                vars.to.regress     = vars_to_regress)

  # head(): if fewer excluded genes were among the variable features, trim to exactly n_features
  Seurat::VariableFeatures(cg_srt, assay = "SCT") <-
    head(setdiff(Seurat::VariableFeatures(cg_srt, assay = "SCT"), exclude_genes), n_features)

  cg_srt <- Seurat::RunPCA(cg_srt,
                           assay          = "SCT",
                           reduction.name = "pca.SCT")
  cg_srt
}

# ---------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------

#' Build kNN/SNN graph and find clusters at multiple resolutions.
#' Also stores the NN object (required for umap.nn).
#'
#' @param dims_i  PCA dimensions to use (e.g. seq(30)).
#' @param ref_assay  Assay name — reduction used is paste0("pca.", ref_assay).
#' @param res     Scalar or vector of clustering resolutions.
#' @param k       Number of nearest neighbours.
find_neighbors_clusters <- function(cg_srt, dims_i, ref_assay, res, k) {
  reduction    <- paste0("pca.", ref_assay)
  graph_prefix <- paste0(c("kNN_", "sNN_"), reduction, "_", k)

  cg_srt <- Seurat::FindNeighbors(cg_srt,
                                  dims        = dims_i,
                                  reduction   = reduction,
                                  k.param     = k,
                                  nn.method   = "rann", #Exact search returns the same neighbours
                                  compute.SNN = TRUE,
                                  graph.name  = graph_prefix)

  cg_srt <- Seurat::FindClusters(cg_srt,
                                 resolution       = res,
                                 group.singletons = FALSE,
                                 graph.name       = paste0("sNN_", reduction, "_", k))

  # store NN object — required for umap.nn
  cg_srt <- Seurat::FindNeighbors(cg_srt,
                                  dims            = dims_i,
                                  reduction       = reduction,
                                  nn.method   = "rann", #Exact search returns the same neighbours
                                  graph.name      = paste0("neighbors_", reduction, "_", k),
                                  return.neighbor = TRUE)

  return(cg_srt)
}

# ---------------------------------------------------------------
# Full single-sample processing wrapper
# ---------------------------------------------------------------

#' SCTransform → PCA → SNN clustering → UMAP (nn-based)
#'
#' @param cg_merge_clean  QC-filtered merged Seurat object.
#' @param ref_assay       Assay; must be "SCT".
#' @param dims_i          PCA dimensions to use for neighbour graph and UMAP.
#' @param k               Number of nearest neighbours for SNN graph.
#' @param res             Clustering resolutions (scalar or vector).
#' @param n_features      Number of variable features for SCTransform.
#' @param k.nn            k for umap.nn.
#' @param exclude_patterns Regex patterns to exclude from variable features.
#' @param exclude_feats   Exact gene names to exclude from variable features.
process_all <- function(cg_merge_clean,
                        ref_assay,
                        dims_i,
                        k,
                        res,
                        n_features,
                        k.nn,
                        exclude_genes   = character(0),
                        vars_to_regress = NULL) {

  srt <- preprocess_seurat(cg_merge_clean,
                           ref_assay       = ref_assay,
                           n_features      = n_features,
                           exclude_genes   = exclude_genes,
                           vars_to_regress = vars_to_regress)

  srt <- find_neighbors_clusters(srt,
                                 dims_i    = dims_i,
                                 ref_assay = ref_assay,
                                 k         = k,
                                 res       = res)

  srt <- Seurat::RunUMAP(srt,
                  n.neighbors    = k.nn,
                  nn.name        = Seurat::Neighbors(srt),
                  reduction.name = "umap.nn",
                  return.model   = TRUE)
  srt
}

#' Adjusted Rand Index between a reference column and one or more comparison
#' columns in a metadata data.frame.
#'
#' For each column in \code{compare_cols}, computes
#' ARI(ref_col, compare_col) using only rows where both values are non-NA.
#'
#' @param meta         data.frame of cell metadata (e.g. \code{srt@meta.data}).
#' @param ref_col      Name of the ground-truth column (e.g. \code{"Phase"}).
#' @param compare_cols Character vector of column names to compare against
#'                     \code{ref_col} (e.g. cluster columns at several
#'                     resolutions).
#' @return data.frame with one row per element of \code{compare_cols}:
#'         \code{compare_col}, \code{n_levels}, \code{ari}.
compute_ari <- function(meta, ref_col, compare_cols) {
  missing <- setdiff(c(ref_col, compare_cols), colnames(meta))
  if (length(missing))
    stop("Column(s) not found in meta: ", paste(missing, collapse = ", "))

  ref_vec <- meta[[ref_col]]

  do.call(rbind, lapply(compare_cols, function(cc) {
    cmp_vec <- as.character(meta[[cc]])
    keep    <- !is.na(cmp_vec) & !is.na(ref_vec)
    data.frame(
      compare_col = cc,
      n_levels    = length(unique(cmp_vec[keep])),
      ari         = mclust::adjustedRandIndex(ref_vec[keep], cmp_vec[keep]),
      stringsAsFactors = FALSE
    )
  }))
}

#' Cell-cycle-corrected processing wrapper around process_all().
#'
#' Applies one of three cell-cycle regression strategies before SCTransform:
#' \describe{
#'   \item{"none"}{No regression (identical to plain process_all).}
#'   \item{"full"}{Regress out S.Score and G2M.Score.}
#'   \item{"diff"}{Add CC.Difference = S.Score - G2M.Score; regress that out.
#'         Preserves proliferation signal while removing G1 vs cycling bias.}
#' }
#' S.Score and G2M.Score must already be present in \code{srt@meta.data}
#' (e.g. from \code{Seurat::CellCycleScoring()}).
#'
#' @inheritParams process_all
#' @param cc_strategy One of "none", "full", or "diff".
#' @return Processed Seurat object (identical structure to process_all output).
process_all_cc <- function(srt,
                           ref_assay,
                           dims_i,
                           k,
                           res,
                           n_features,
                           k.nn,
                           exclude_genes = character(0),
                           cc_strategy   = c("full", "diff", "none")) {
  cc_strategy <- match.arg(cc_strategy)

  vars_reg <- switch(cc_strategy,
    none = NULL,
    full = c("S.Score", "G2M.Score"),
    diff = {
      srt[["CC.Difference"]] <- srt$S.Score - srt$G2M.Score
      "CC.Difference"
    }
  )

  process_all(srt,
              ref_assay       = ref_assay,
              dims_i          = dims_i,
              k               = k,
              res             = res,
              n_features      = n_features,
              k.nn            = k.nn,
              exclude_genes   = exclude_genes,
              vars_to_regress = vars_reg)
}

# ---------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------

# Default discrete colour palette for UMAP plots (murine BM cell types)
# Cell-type colour palette — loaded from metadata/project_config.yaml at source time
# (from the project root, where the pipeline and the reports run).
# Use UMAP_COLOURS["CellType"] in plotting functions.
UMAP_COLOURS <- unlist(yaml::read_yaml("metadata/project_config.yaml")$plotting$colours)

#' Generate n maximally-distinct colours for discrete UMAP plots.
#' Uses Polychrome::palette36 (via paletteer), skipping the near-black and
#' near-white entries (positions 1–2), then interleaved to maximise perceptual
#' distance between consecutive indices.
#' Falls back to RColorBrewer Set1/Set2/Set3 if paletteer is not installed.
.discrete_colours <- function(n) {
  if (requireNamespace("paletteer", quietly = TRUE)) {
    base <- as.character(
      paletteer::paletteer_d("Polychrome::palette36")
    )[-c(1, 2)]   # drop near-black (#5A5156) and near-white (#E4E1E3)
    # Interleave: take every other index from two halves so that visually
    # similar neighbours in the original sequence end up far apart.
    half  <- ceiling(length(base) / 2)
    idx   <- as.integer(rbind(seq_len(half),
                               seq(half + 1, length(base), length.out = half)))
    idx   <- idx[idx <= length(base)]
    base  <- base[idx]
    if (n <= length(base)) return(base[seq_len(n)])
    return(rep_len(base, n))
  }
  # fallback: combine RColorBrewer qualitative palettes
  sets <- c(RColorBrewer::brewer.pal(9, "Set1"),
            RColorBrewer::brewer.pal(8, "Set2"),
            RColorBrewer::brewer.pal(12, "Set3"))
  rep_len(unique(sets), n)
}

#' Add custom UMAP-style axes (arrow annotations replacing axis lines).
#' Returns a list of ggplot2 annotate() calls to be added to a plot.
.umap_axes <- function(d, xcol, ycol, label_size = 3, arrow_len = 2) {
  xmin <- min(d[[xcol]], na.rm = TRUE)
  ymin <- min(d[[ycol]], na.rm = TRUE)
  arr  <- ggplot2::arrow(angle = 30, ends = "last", type = "closed",
                         length = ggplot2::unit(0.1, "inches"))
  list(
    ggplot2::annotate("segment",
                      x = xmin - 1, xend = xmin + arrow_len,
                      y = ymin - 1, yend = ymin - 1,
                      colour = "black", linewidth = 0.6, arrow = arr),
    ggplot2::annotate("segment",
                      x = xmin - 1, xend = xmin - 1,
                      y = ymin - 1, yend = ymin + arrow_len,
                      colour = "black", linewidth = 0.6, arrow = arr),
    ggplot2::annotate("text", size = label_size, label = "UMAP 1",
                      x = xmin + 0.2, y = ymin - 1.6),
    ggplot2::annotate("text", size = label_size, label = "UMAP 2",
                      x = xmin - 1.6, y = ymin + 0.2, angle = 90)
  )
}

#' Base UMAP theme — no axes, no grid, clean white background.
.umap_theme <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio        = 1,
      legend.title        = ggplot2::element_text(size = 14),
      legend.text         = ggplot2::element_text(size = 13),
      plot.title          = ggplot2::element_text(size = 14, hjust = 0.5)
    )
}

#' Extract plot data from a Seurat DimPlot or FeaturePlot and rebuild as a
#' lean ggplot2 object that carries only the minimal coordinate data.
#'
#' Two-layer rendering (black outline + colour on top) with custom arrow axes,
#' inspired by Shoeb2025_CG-AMKL style. Serializes to < 2 MB (no Seurat env).
#'
#' @param p            ggplot from Seurat::DimPlot or ::FeaturePlot.
#' @param colour_col   Column in p$data used for colour.
#' @param cols         Named colour vector (NULL = ggplot2 default / viridis for numeric).
#' @param legend_title Legend title string.
#' @param pt_size      Size of coloured points (outline layer is pt_size * 2).
lean_dimplot <- function(p, colour_col, cols = NULL,
                         legend_title = colour_col, pt_size = 0.7,
                         label = FALSE, label_size = 3,
                         viridis_option = "mako",
                         cell_order = "shuffle",
                         axes = NULL,
                         outline = TRUE) {
  d    <- p$data
  xcol <- names(d)[1]
  ycol <- names(d)[2]
  # Seurat sanitizes column names (e.g. "Ighv5-4" → "Ighv5.4", "mt-Co3" → "mt.Co3").
  # If colour_col is not found verbatim, prefer the first numeric non-embedding column
  # (FeaturePlot case) before falling back to column [3] (DimPlot discrete case).
  if (!colour_col %in% names(d)) {
    non_emb <- names(d)[-(1:2)]
    numeric_cols <- non_emb[sapply(non_emb, function(cn) is.numeric(d[[cn]]))]
    colour_col <- if (length(numeric_cols) > 0) numeric_cols[1] else non_emb[1]
  }
  # FeaturePlot stores values as numeric; DimPlot stores groups as factor/char.
  # Factors are always treated as discrete — even when levels are numeric strings
  # (e.g. cluster ids "0","1","2"). Only plain numeric vectors use the continuous scale.
  col_vals <- d[[colour_col]]
  is_num <- is.numeric(col_vals) && !is.factor(col_vals)

  if (cell_order == "shuffle") {
    set.seed(42)
    d <- d[sample(nrow(d)), ]
  } else if (cell_order == "expression" && is_num) {
    # High-expressing cells last so they render on top
    d <- d[order(d[[colour_col]], na.last = FALSE), ]
  }

  # Centroid data for labels (only used when label = TRUE and colour is discrete)
  if (label && !is_num) {
    d_lab <- do.call(rbind, lapply(split(d, d[[colour_col]]), function(g) {
      data.frame(
        x     = median(g[[xcol]]),
        y     = median(g[[ycol]]),
        label = as.character(g[[colour_col]][1]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    d_lab <- NULL
  }

  # Pre-compute discrete colours outside eval() — .discrete_colours() is not
  # visible inside the eval env whose parent is getNamespace("ggplot2").
  if (!is_num && is.null(cols)) {
    cols <- .discrete_colours(length(unique(d[[colour_col]])))
  }

  args <- list(d = d, d_lab = d_lab, xcol = xcol, ycol = ycol,
               colour_col = colour_col, legend_title = legend_title,
               pt_size = pt_size, cols = cols, is_num = is_num,
               label = label, label_size = label_size,
               viridis_option = viridis_option, axes = axes,
               outline = outline)

  eval(quote({
    p2 <- ggplot2::ggplot(d,
            ggplot2::aes(x = !!rlang::sym(xcol), y = !!rlang::sym(ycol))) +
      # black outline layer
      (if (outline) ggplot2::geom_point(colour = "black", size = pt_size * 2,
                                        stroke = 0, show.legend = FALSE)) +
      # coloured layer on top
      ggplot2::geom_point(ggplot2::aes(colour = !!rlang::sym(colour_col)),
                          size = pt_size, stroke = 0) +
      ggplot2::labs(colour = legend_title, x = NULL, y = NULL) +
      .umap_theme()

    if (is_num) {
      p2 <- p2 +
        ggplot2::scale_colour_viridis_c(option = viridis_option, na.value = "grey85") +
        ggplot2::guides(colour = ggplot2::guide_colourbar(
          barwidth = 0.8, barheight = 6))
    } else if (!is.null(cols)) {
      p2 <- p2 +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(size = 6)))
    }

    if (label && !is.null(d_lab)) {
      x_range <- range(d[[xcol]], na.rm = TRUE)
      y_range <- range(d[[ycol]], na.rm = TRUE)
      d_lab_vis <- d_lab[
        d_lab$x >= x_range[1] & d_lab$x <= x_range[2] &
        d_lab$y >= y_range[1] & d_lab$y <= y_range[2], ]
      p2 <- p2 + ggplot2::geom_label(
        data         = d_lab_vis,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes  = FALSE,
        size         = label_size,
        label.size   = 0,
        fill         = scales::alpha("white", 0.7)
      )
    }

    # append arrow axes — use pre-computed axes if supplied (e.g. from split panels)
    axes_to_add <- if (!is.null(axes)) axes else .umap_axes(d, xcol, ycol)
    for (ann in axes_to_add) p2 <- p2 + ann
    p2
  }), envir = list2env(args, parent = getNamespace("ggplot2")))
}

#' DimPlot wrapper — returns a lean ggplot (< 2 MB serialized)
plot_umap_groups <- function(srt_obj,
                             group_i,
                             cols         = NULL,
                             reduction    = "umap.nn",
                             legend_title = group_i,
                             pt.size      = 0.7,
                             label        = FALSE,
                             label_size   = 3) {
  p <- Seurat::DimPlot(srt_obj, group.by = group_i,
                       reduction = reduction, pt.size = pt.size)
  lean_dimplot(p, colour_col = group_i, cols = cols,
               legend_title = legend_title, pt_size = pt.size,
               label = label, label_size = label_size)
}

#' FeaturePlot wrapper — returns a lean ggplot (< 2 MB serialized)
plot_umap_feature <- function(srt_obj,
                              feature,
                              reduction      = "umap.nn",
                              pt.size        = 0.7,
                              viridis_option = "mako",
                              cell_order     = c("shuffle", "expression")) {
  cell_order <- match.arg(cell_order)
  p <- Seurat::FeaturePlot(srt_obj, features = feature,
                           reduction = reduction, pt.size = pt.size)
  lean_dimplot(p, colour_col = feature, cols = NULL,
               legend_title = feature, pt_size = pt.size,
               viridis_option = viridis_option, cell_order = cell_order)
}

.umap_highlight <- function(srt_obj, cells.highlight, reduction, pt.size,
                             legend_title) {
  emb  <- as.data.frame(Seurat::Embeddings(srt_obj, reduction))
  xcol <- colnames(emb)[1]; ycol <- colnames(emb)[2]
  emb$highlighted <- rownames(emb) %in% cells.highlight
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = emb[!emb$highlighted, ],
                        ggplot2::aes(.data[[xcol]], .data[[ycol]]),
                        colour = "grey80", size = pt.size, stroke = 0) +
    ggplot2::geom_point(data = emb[emb$highlighted, ],
                        ggplot2::aes(.data[[xcol]], .data[[ycol]]),
                        colour = "#e31a1c", size = pt.size * 1.5, stroke = 0) +
    ggplot2::labs(title = legend_title %||% "", x = NULL, y = NULL) +
    .umap_theme()
  for (ann in .umap_axes(emb, xcol, ycol)) p <- p + ann
  p
}

.umap_split <- function(srt_obj, group_i, cols, reduction, legend_title,
                         pt.size, bg.size, split.by, label, label_size) {
  p_full <- Seurat::DimPlot(srt_obj, group.by = group_i, reduction = reduction,
                             pt.size = pt.size, label = FALSE)
  d_full     <- p_full$data
  colour_col <- if (group_i %in% colnames(d_full)) group_i else "ident"
  d_full[[split.by]] <- srt_obj@meta.data[rownames(d_full), split.by]

  all_levels <- levels(d_full[[colour_col]]) %||%
                sort(unique(as.character(d_full[[colour_col]])))
  if (is.null(cols)) cols <- stats::setNames(.discrete_colours(length(all_levels)), all_levels)

  xcol     <- colnames(d_full)[1]
  ycol     <- colnames(d_full)[2]
  levels_x <- unique(d_full[[split.by]])
  # Pre-compute axes from full embedding so all panels share the same origin
  shared_axes <- .umap_axes(d_full, xcol, ycol)
  bg_black <- ggplot2::geom_point(
    data        = d_full,
    mapping     = ggplot2::aes(.data[[xcol]], .data[[ycol]]),
    colour      = "black", size = bg.size * 2, stroke = 0,
    inherit.aes = FALSE
  )
  bg_grey <- ggplot2::geom_point(
    data        = d_full,
    mapping     = ggplot2::aes(.data[[xcol]], .data[[ycol]]),
    colour      = "grey85", size = bg.size, stroke = 0,
    inherit.aes = FALSE
  )

  plots <- lapply(seq_along(levels_x), function(i) {
    lvl      <- levels_x[i]
    p_sub    <- p_full
    p_sub$data <- d_full[d_full[[split.by]] == lvl, ]
    pan <- lean_dimplot(p_sub, colour_col = colour_col, cols = cols,
                        legend_title = legend_title %||% group_i,
                        pt_size = pt.size, label = label,
                        label_size = label_size, axes = shared_axes) +
      ggplot2::labs(title = lvl)
    pan$layers <- c(list(bg_black, bg_grey), pan$layers)
    if (i < length(levels_x)) pan <- pan + ggplot2::theme(legend.position = "none")
    pan
  })
  patchwork::wrap_plots(plots)
}

#' DimPlot wrapper with split.by / cells.highlight support.
plot_umap_groups_ext <- function(srt_obj,
                                 group_i         = NULL,
                                 cols            = NULL,
                                 reduction       = "umap.nn",
                                 legend_title    = group_i,
                                 pt.size         = 0.7,
                                 bg.size         = pt.size,
                                 split.by        = NULL,
                                 cells.highlight = NULL,
                                 label           = FALSE,
                                 label_size      = 3) {
  if (!is.null(cells.highlight))
    return(.umap_highlight(srt_obj, cells.highlight, reduction, pt.size, legend_title))
  if (!is.null(split.by))
    return(.umap_split(srt_obj, group_i, cols, reduction, legend_title,
                       pt.size, bg.size, split.by, label, label_size))
  plot_umap_groups(srt_obj, group_i = group_i, cols = cols,
                   reduction = reduction, legend_title = legend_title %||% group_i,
                   pt.size = pt.size, label = label, label_size = label_size)
}

#' FeaturePlot wrapper with split.by / multiple features support.
#' Multiple features are returned as a patchwork grid.
#' split.by facets each feature panel by a metadata column.
plot_umap_feature_ext <- function(srt_obj,
                                  features,
                                  reduction      = "umap.nn",
                                  pt.size        = 0.7,
                                  split.by       = NULL,
                                  viridis_option = "mako",
                                  cell_order     = c("shuffle", "expression"),
                                  outline        = TRUE,
                                  ncol           = NULL) {
  cell_order <- match.arg(cell_order)
  if (!is.null(split.by)) {
    # Extract embedding data from the full object per feature — avoids subsetting
    # which can silently drop reductions in Seurat v5.
    plots <- lapply(features, function(feat) {
      p_full <- Seurat::FeaturePlot(srt_obj, features = feat,
                                    reduction = reduction, pt.size = pt.size)
      d_full <- p_full$data
      d_full[[split.by]] <- srt_obj@meta.data[rownames(d_full), split.by]
      levels_x <- unique(d_full[[split.by]])
      row_plots <- lapply(seq_along(levels_x), function(i) {
        lvl   <- levels_x[i]
        d_sub <- d_full[d_full[[split.by]] == lvl, ]
        p_sub <- p_full
        p_sub$data <- d_sub
        pan <- lean_dimplot(p_sub, colour_col = feat, cols = NULL,
                            legend_title = feat, pt_size = pt.size,
                            viridis_option = viridis_option,
                            cell_order = cell_order,
                            axes = list(), outline = outline) +
          ggplot2::labs(title = paste0(feat, " — ", lvl))
        # suppress legend on all panels except the last in the row
        if (i < length(levels_x)) pan <- pan + ggplot2::theme(legend.position = "none")
        pan
      })
      patchwork::wrap_plots(row_plots, nrow = 1)
    })
    return(patchwork::wrap_plots(plots, ncol = 1))
  }
  plots <- lapply(features, function(feat) {
    p <- Seurat::FeaturePlot(srt_obj, features = feat,
                             reduction = reduction, pt.size = pt.size)
    lean_dimplot(p, colour_col = feat, cols = NULL,
                 legend_title = feat, pt_size = pt.size,
                 viridis_option = viridis_option, cell_order = cell_order,
                 axes = list(), outline = outline)
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

#' Extract a tidy data frame of PCA embeddings + metadata for PCA distribution plots.
#'
#' Returns a data.frame with columns PC1…PCn, sample.name, Phase, and the active
#' cluster identity. Intended as a lightweight targets output (no Seurat object).
#'
#' @param srt       Seurat object with PCA and cluster metadata.
#' @param pca_name  Name of the PCA reduction (e.g. \code{"pca.SCT"}).
#' @param n_pcs     Number of PCs to retain.
#' @return data.frame
extract_pca_df <- function(srt, pca_name, n_pcs) {
  pca_df <- as.data.frame(Seurat::Embeddings(srt, pca_name)[, seq(n_pcs)])
  colnames(pca_df)    <- paste0("PC", seq(n_pcs))
  pca_df$sample.name  <- srt@meta.data$sample.name
  pca_df$Phase        <- srt@meta.data$Phase
  pca_df$cluster      <- as.character(Seurat::Idents(srt))
  pca_df$UMAP1        <- Seurat::Embeddings(srt, "umap.nn")[, 1]
  pca_df$UMAP2        <- Seurat::Embeddings(srt, "umap.nn")[, 2]
  pca_df
}

#' Parallel-coordinates plot over top PCs, coloured by a grouping variable.
#'
#' Lines are individual cells (downsampled to \code{n_max}); bold lines are
#' per-group means overlaid on top.
#'
#' @param df         data.frame from \code{extract_pca_df}.
#' @param colour_var Column name to colour by.
#' @param title      Plot title.
#' @param palette    Optional named colour vector.
#' @param n_max      Max cells to draw as background lines (default 5000).
plot_pc_parcoor <- function(df, colour_var, title, palette = NULL, n_max = 5000L) {
  pc_cols  <- grep("^PC", colnames(df), value = TRUE)
  pc_range <- apply(df[, pc_cols, drop = FALSE], 2, range)

  avg_long <- df |>
    tidyr::pivot_longer(dplyr::all_of(pc_cols), names_to = "variable", values_to = "value") |>
    dplyr::group_by(.data[[colour_var]], variable) |>
    dplyr::summarise(value = mean(value), .groups = "drop") |>
    dplyr::mutate(variable = factor(variable, levels = pc_cols),
                  value    = (value - pc_range[1, variable]) /
                               (pc_range[2, variable] - pc_range[1, variable]))

  if (nrow(df) > n_max) df <- df[sample(nrow(df), n_max), ]

  p <- GGally::ggparcoord(
    data        = df,
    columns     = match(pc_cols, colnames(df)),
    groupColumn = colour_var,
    scale       = "uniminmax",
    alpha       = 0.05,
    showPoints  = FALSE
  ) +
    ggplot2::geom_point(data = avg_long,
               ggplot2::aes(x = variable, y = value, colour = .data[[colour_var]]),
               size = 2.5, shape = 19, alpha = 1, inherit.aes = FALSE) +
    ggplot2::geom_line(data = avg_long,
              ggplot2::aes(x = variable, y = value,
                  colour = .data[[colour_var]], group = .data[[colour_var]]),
              linewidth = 0.8, alpha = 0.9, inherit.aes = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = "scaled embedding", colour = colour_var) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
          legend.text     = ggplot2::element_text(size = 7),
          axis.text.x     = ggplot2::element_text(angle = 0)) +
    ggplot2::guides(color = ggplot2::guide_legend(
        override.aes = list(linewidth = 3, alpha = 1), nrow = 1))

  if (!is.null(palette)) p <- p + ggplot2::scale_colour_manual(values = palette)
  p
}

#' Grid of UMAP panels coloured by PC score for each of the top PCs.
#'
#' @param df     data.frame from \code{extract_pca_df} (must contain UMAP1, UMAP2).
#' @param n_pcs  Number of PCs to plot.
#' @param n_cols Number of columns in the patchwork grid (default 5).
#' @param n_max  Max cells per panel (default 50000).
plot_pc_umap_grid <- function(df, n_pcs, n_cols = 5L, n_max = 50000L) {
  umap_df <- df[, c("UMAP1", "UMAP2")]
  plots <- lapply(seq(n_pcs), function(i) {
    scores <- df[[paste0("PC", i)]]
    d      <- cbind(umap_df, score = scores)
    if (nrow(d) > n_max) d <- d[sample(nrow(d), n_max), ]
    d <- d[order(abs(d$score)), ]
    lim <- max(abs(d$score)) * c(-1, 1)
    ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = score)) +
      ggplot2::geom_point(size = 0.3, alpha = 0.6) +
      ggplot2::scale_colour_gradientn(
        colours = c("#2166ac", "grey90", "#d73027"),
        limits  = lim, name = "score") +
      ggplot2::labs(title = paste0("PC", i)) +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(plot.title       = ggplot2::element_text(size = 8, hjust = 0.5),
            legend.key.width = ggplot2::unit(0.3, "cm"),
            legend.text      = ggplot2::element_text(size = 6))
  })
  patchwork::wrap_plots(plots, ncol = n_cols)
}

#' Distribution of selected metadata levels along embedding components.
#'
#' For each requested component, shows how cells belonging to
#' \code{levels_interest} in \code{group_col} are distributed along that
#' component, split by \code{split_col}. Useful for checking whether a
#' specific cluster sits differently across conditions along PCA/UMAP axes.
#'
#' @param reduction      A Seurat \code{DimReduc} object \emph{or} a cell ×
#'                       component matrix / data.frame (rownames = cell IDs).
#' @param metadata       data.frame of cell metadata (rownames = cell IDs).
#' @param group_col      Metadata column whose levels are filtered
#'                       (e.g. a cluster column).
#' @param levels_interest Character vector of levels to keep from
#'                       \code{group_col} (e.g. \code{c("3", "5")}).
#' @param split_col      Metadata column that defines the groups whose
#'                       distributions are compared (e.g. \code{"sample.name"}).
#' @param n_components   Number of leading components to plot (default 10).
#' @param plot_type      \code{"violin"} (default) or \code{"boxplot"}.
#' @param palette        Optional named colour vector for \code{split_col}
#'                       levels. If \code{NULL}, uses the default ggplot2
#'                       palette.
#' @param ncol           Number of facet columns (default 2).
#' @return A \code{ggplot} object.
plot_embedding_dist <- function(reduction,
                                metadata,
                                group_col,
                                levels_interest,
                                split_col,
                                n_components = 10L,
                                plot_type    = c("violin", "boxplot"),
                                palette      = NULL,
                                ncol         = 2L) {
  plot_type <- match.arg(plot_type)

  # ── Extract embedding matrix ───────────────────────────────────────────────
  if (inherits(reduction, "DimReduc")) {
    emb <- Seurat::Embeddings(reduction)
  } else {
    emb <- as.matrix(reduction)
  }

  n_components <- min(n_components, ncol(emb))
  emb          <- emb[, seq_len(n_components), drop = FALSE]

  # ── Align cells ───────────────────────────────────────────────────────────
  common <- intersect(rownames(emb), rownames(metadata))
  if (length(common) == 0L)
    stop("No shared cell IDs between reduction and metadata.")
  emb  <- emb[common, , drop = FALSE]
  meta <- metadata[common, , drop = FALSE]

  # ── Filter to levels of interest ──────────────────────────────────────────
  keep <- as.character(meta[[group_col]]) %in% as.character(levels_interest)
  if (!any(keep))
    stop("None of levels_interest found in column '", group_col, "'.")
  emb  <- emb[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  # ── Long format ───────────────────────────────────────────────────────────
  split_vec <- meta[[split_col]]
  df <- as.data.frame(emb) |>
    tibble::rownames_to_column("cell_id") |>
    tidyr::pivot_longer(-cell_id, names_to = "component", values_to = "value") |>
    dplyr::mutate(
      component    = factor(component, levels = colnames(emb)),
      !!split_col := split_vec[match(cell_id, rownames(meta))]
    )

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x    = .data[[split_col]],
    y    = value,
    fill = .data[[split_col]]
  ))

  if (plot_type == "violin") {
    p <- p +
      ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.8, colour = NA) +
      ggplot2::geom_boxplot(width = 0.08, outlier.size = 0,
                            fill = "white", colour = "grey30", alpha = 0.7)
  } else {
    p <- p +
      ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.8)
  }

  if (!is.null(palette))
    p <- p + ggplot2::scale_fill_manual(values = palette)

  p +
    ggplot2::facet_wrap(~component, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      x        = NULL,
      y        = "Embedding value",
      title    = paste0(group_col, " \u2208 {",
                        paste(levels_interest, collapse = ", "), "}"),
      subtitle = paste0("Split by: ", split_col,
                        "  \u00b7  ", n_components, " components")
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}

#' UMAP grid: sample.name (rows) × dm2_class (columns).
#'
#' Returns a list with two ggplots:
#'   \code{$plain}     — all cells grey; shows UMAP structure per combination.
#'   \code{$highlight} — gghighlight: each facet colours its own cells by
#'                       dm2_class; non-matching cells dimmed to grey.
#'
#' @param srt       Seurat object with \code{sample.name} and \code{dm2_class}
#'                  in \code{@meta.data}.
#' @param reduction Reduction name (default \code{"umap.nn"}).
#' @return Named list: \code{plain}, \code{highlight}.
plot_umap_grid_sample_dm2 <- function(srt, reduction = "umap.nn") {
  emb  <- as.data.frame(Seurat::Embeddings(srt, reduction = reduction))
  xcol <- colnames(emb)[1]
  ycol <- colnames(emb)[2]

  meta <- srt@meta.data[rownames(emb), c("sample.name", "dm2_class")]
  d    <- cbind(emb, meta)
  d$sample.name <- factor(d$sample.name)
  d$dm2_class   <- factor(d$dm2_class)

  base_theme <-
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio  = 1,
      strip.text    = ggplot2::element_text(size = 9),
      plot.title    = ggplot2::element_text(size = 11, hjust = 0.5),
      panel.spacing = ggplot2::unit(0.2, "lines")
    )

  p_plain <- ggplot2::ggplot(d,
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])) +
    ggplot2::geom_point(colour = "grey40", size = 0.3, stroke = 0, alpha = 0.4) +
    ggplot2::facet_grid(sample.name ~ dm2_class) +
    ggplot2::labs(title = "sample \u00d7 dm2_class (plain)") +
    base_theme

  # Manual highlight: grey background layer uses data stripped of the faceting
  # columns so ggplot2 does not subset it per panel — every panel sees all
  # points as grey. The coloured foreground layer uses the full data and IS
  # subsetted by facet_grid, so each panel only shows its own cells in colour.
  d_bg <- d[, c(xcol, ycol)]   # no sample.name / dm2_class → not facet-subsetted
  p_highlight <- ggplot2::ggplot(d,
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])) +
    ggplot2::geom_point(data = d_bg, colour = "grey85", size = 0.3, alpha = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = dm2_class), size = 0.3) +
    ggplot2::facet_grid(sample.name ~ dm2_class) +
    ggplot2::labs(title = "sample \u00d7 dm2_class (highlight)", colour = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 4, alpha = 1))) +
    base_theme

  list(plain = p_plain, highlight = p_highlight)
}

#' Dot plot from flat data.frame targets (no Seurat needed).
#' @param expr_df  anno_expr data.frame: cell_id + one col per gene (SCT).
#' @param meta_df  anno_meta data.frame: cell_id + res_* cluster columns.
#' @param genes    Character vector of genes to show.
#' @param res      Clustering resolution scalar (e.g. 0.1).
#' @param title    Plot title (NULL = "res <res>").
#' @return ggplot or NULL if genes / cluster column absent.
plot_dot_df <- function(expr_df, meta_df, genes, res, title = NULL) {
  clust_col <- paste0("res_", gsub("\\.", "_", res))
  if (!clust_col %in% colnames(meta_df)) return(invisible(NULL))
  genes_ok <- intersect(genes, colnames(expr_df))
  if (length(genes_ok) == 0) return(invisible(NULL))

  d <- merge(meta_df[, c("cell_id", clust_col)],
             expr_df[, c("cell_id", genes_ok), drop = FALSE],
             by = "cell_id")
  cluster <- d[[clust_col]]

  dot_df <- do.call(rbind, lapply(genes_ok, function(g) {
    x   <- d[[g]]
    avg <- tapply(x, cluster, mean)
    pct <- tapply(x, cluster, function(v) mean(v > 0) * 100)
    data.frame(gene = g, cluster_id = names(avg),
               avg_expr = as.numeric(avg), pct_expr = as.numeric(pct),
               stringsAsFactors = FALSE)
  }))
  dot_df$cluster_id <- factor(dot_df$cluster_id,
    levels = sort(unique(dot_df$cluster_id), method = "radix"))
  dot_df$gene <- factor(dot_df$gene, levels = rev(genes_ok))

  ggplot2::ggplot(dot_df,
                  ggplot2::aes(cluster_id, gene,
                               size = pct_expr, colour = avg_expr)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "grey90", high = "#08306b",
                                   name = "Avg expr") +
    ggplot2::scale_size_continuous(range = c(0.5, 5), name = "% expressed") +
    ggplot2::labs(x = "Cluster", y = NULL,
                  title = if (!is.null(title)) title else paste("res", res)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   axis.text.y = ggplot2::element_text(size = 8))
}

#' Save a ggplot to both PNG (for HTML report) and PDF (for manuscript editing).
#' Returns the PNG path, which is what targets caches as the file target value.
#'
#' @param fig    A ggplot object.
#' @param stem   Full path without extension, e.g. "output/run/umap/umap_sample".
#' @param width  Width in inches.
#' @param height Height in inches.
#' @param dpi    Resolution for PNG (default 200).
#' @return Invisible PNG path (character).
save_plot <- function(fig, stem, width, height, dpi = 200) {
  ggplot2::ggsave(paste0(stem, ".pdf"), fig, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), fig, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(paste0(stem, ".png"))
}

#' VlnPlot wrapper for QC metrics — NOT cached in targets; saved directly to PDF
#'
#' @param md  data.frame with QC metric columns (e.g. qc_metadata from tar_map).
plot_qc_violin <- function(md,
                           features = c("nFeature_RNA", "nCount_RNA",
                                        "percent.mt", "percent.rp",
                                        "percent.hb", "count_feat_ratio",
                                        "percent.stress"),
                           group_by = "sample.name",
                           ncol     = 2) {
  n_bottom <- ncol

  plots <- lapply(seq_along(features), function(i) {
    m       <- features[i]
    show_x  <- i > (length(features) - n_bottom)
    p <- ggplot2::ggplot(md, ggplot2::aes(x = .data[[group_by]], y = .data[[m]])) +
      ggplot2::geom_violin(fill = "steelblue", alpha = 0.5) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.7) +
      ggplot2::labs(title = m, x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x  = if (show_x) ggplot2::element_text(angle = 45, hjust = 1) else ggplot2::element_blank(),
        axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank()
      )
    if (m == "percent.rp") p <- p + ggplot2::scale_y_continuous(limits = c(0, 100))
    p
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

# ---------------------------------------------------------------
# Summary helpers for Quarto report
# ---------------------------------------------------------------

#' Summarise QC metrics per sample (pre-filter)
#'
#' Returns median and IQR for nFeature_RNA, nCount_RNA, percent.mt,
#' percent.rp, percent.hb per sample. Called on a metadata data.frame
#' (e.g. bound qc_metadata_* tables from tar_map).
#'
#' @param meta  data.frame with QC metrics columns and sample.name.
#' @return data.frame with one row per sample.
summarise_qc_metrics <- function(meta) {
  metrics <- c("nFeature_RNA", "nCount_RNA",
               "percent.mt", "percent.rp", "percent.hb")
  meta <- meta[, c("sample.name", metrics)]

  result <- do.call(rbind, lapply(
    split(meta, meta$sample.name),
    function(df) {
      row <- data.frame(sample = df$sample.name[1], stringsAsFactors = FALSE)
      for (m in metrics) {
        row[[paste0(m, "_median")]] <- round(median(df[[m]], na.rm = TRUE), 2)
        row[[paste0(m, "_IQR")]]    <- round(IQR(df[[m]],    na.rm = TRUE), 2)
      }
      row
    }
  ))
  rownames(result) <- NULL
  result
}

#' Summarise doublet counts and rates per sample.
#' Only scored cells (scDblFinder.class != "unscored") are included.
#'
#' @param meta  data.frame with pass_dbl and sample.name columns (qc_metadata).
#' @return data.frame: sample, n_singlet, n_doublet, doublet_rate_pct.
summarise_doublets <- function(meta) {
  meta    <- meta[meta$scDblFinder.class != "unscored", ]
  tbl     <- table(meta$sample.name, meta$pass_dbl)
  singlet <- if (TRUE  %in% colnames(tbl)) as.integer(tbl[, "TRUE"])  else rep(0L, nrow(tbl))
  doublet <- if (FALSE %in% colnames(tbl)) as.integer(tbl[, "FALSE"]) else rep(0L, nrow(tbl))
  data.frame(
    sample           = rownames(tbl),
    n_singlet        = singlet,
    n_doublet        = doublet,
    doublet_rate_pct = round(100 * doublet / (singlet + doublet), 2),
    stringsAsFactors = FALSE
  )
}

#' Summarise cluster sizes per resolution from a processed Seurat object
#'
#' @param srt  Seurat object (mpn_processed).
#' @return Named list of data.frames; names are resolution labels (e.g. "res.0.3").
#' Canonical marker bubble plot
#'
#' Reads the canonical_markers.csv produced by the pipeline, scales avg_exp
#' per gene, and returns a ggplot bubble plot (avg expression × % expressing).
#'
#' @param csv_path   Path to canonical_markers.csv.
#' @param gene_order Character vector of genes in desired y-axis order (top-down).
#' @param title      Plot title string.
#' @return ggplot object.
plot_canonical_bubble <- function(csv_path, gene_order, title = "") {
  bubble <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(bubble) == 0) {
    message("plot_canonical_bubble: empty CSV, returning blank plot")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = paste0(title, " (no data)")))
  }

  # keep only genes present in the data, preserving order
  gene_order <- gene_order[gene_order %in% bubble$gene]
  bubble$gene    <- factor(bubble$gene,    levels = rev(unique(gene_order)))
  bubble$cluster <- factor(bubble$cluster, levels = sort(unique(as.integer(bubble$cluster))))

  bubble <- bubble |>
    dplyr::group_by(gene) |>
    dplyr::mutate(avg_exp_scaled = as.numeric(scale(avg_exp))) |>
    dplyr::ungroup()

  ggplot2::ggplot(bubble,
                  ggplot2::aes(x = cluster, y = gene,
                               colour = avg_exp_scaled,
                               size   = pct_exp * 100)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(low      = "steelblue",
                                    mid      = "white",
                                    high     = "firebrick",
                                    midpoint = 0,
                                    name     = "Avg exp\n(scaled)") +
    ggplot2::scale_size_continuous(range = c(0.5, 6), name = "% expressing") +
    ggplot2::labs(title = title, x = "Cluster", y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x      = ggplot2::element_text(size = 9),
                   axis.text.y      = ggplot2::element_text(size = 8),
                   panel.grid.major = ggplot2::element_line(colour = "grey90"))
}

summarise_cluster_sizes <- function(srt) {
  clust_cols <- grep("^sNN_pca\\.",
                     colnames(srt@meta.data), value = TRUE)
  clust_cols <- clust_cols[grepl("_res\\.", clust_cols)]
  if (length(clust_cols) == 0) {
    warning("No cluster columns found matching expected pattern.")
    return(list())
  }
  result <- lapply(clust_cols, function(col) {
    tbl <- sort(table(srt@meta.data[[col]]), decreasing = TRUE)
    data.frame(
      cluster   = names(tbl),
      n_cells   = as.integer(tbl),
      pct_cells = round(100 * as.integer(tbl) / sum(tbl), 2),
      stringsAsFactors = FALSE
    )
  })
  names(result) <- sub("^sNN_pca\\.[^_]+_\\d+_", "", clust_cols)
  result
}

# ---------------------------------------------------------------
# Output helpers — plots and tables
# ---------------------------------------------------------------
# Each function creates the output directory, produces its output,
# and returns the file path for targets format = "file" tracking.

#' Write a data.frame to CSV and return the path.
save_csv <- function(df, dir, filename) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  write.csv(df, path, row.names = FALSE)
  path
}

#' UMAP coloured by sample.name.
plot_umap_sample <- function(srt, subset, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem  <- file.path(out_dir, "umap_sample")
  p_raw <- Seurat::DimPlot(srt, group.by = "sample.name",
                            reduction = "umap.nn", pt.size = 0.7)
  fig   <- lean_dimplot(p_raw, colour_col = "sample.name",
                        legend_title = "Sample") +
    ggplot2::labs(title = paste0("Samples — ", subset))
  save_plot(fig, stem, width = 7, height = 7)
}

#' UMAP coloured by cluster at the given resolution.
plot_umap_cluster <- function(srt, subset, marker_res, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem      <- file.path(out_dir, paste0("umap_cluster_res", marker_res))
  clust_col <- grep(paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", marker_res), "$"),
                    colnames(srt@meta.data), value = TRUE)[1]
  p_raw <- Seurat::DimPlot(srt, group.by = clust_col,
                            reduction = "umap.nn", pt.size = 0.7,
                            label = TRUE, label.size = 4)
  fig   <- lean_dimplot(p_raw, colour_col = clust_col,
                        legend_title = "Cluster") +
    ggplot2::labs(title = paste0("Clusters (res ", marker_res, ") — ", subset))
  save_plot(fig, stem, width = 7, height = 7)
}

#' UMAP coloured by cell cycle phase.
plot_umap_phase <- function(srt, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out_dir, "umap_cell_cycle_phase")
  fig  <- plot_umap_groups(srt,
                            group_i      = "Phase",
                            cols         = c("G1" = "#4daf4a", "S" = "#225ea8",
                                             "G2M" = "#e41a1c"),
                            reduction    = "umap.nn",
                            legend_title = "Cell cycle phase")
  save_plot(fig, stem, width = 7, height = 7)
}

#' UMAP coloured by scDblFinder doublet class.
plot_umap_doublet <- function(srt, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out_dir, "umap_doublet_class")
  fig  <- plot_umap_groups(srt,
                            group_i      = "scDblFinder.class",
                            cols         = c("singlet" = "#2166ac", "doublet" = "#d73027"),
                            reduction    = "umap.nn",
                            legend_title = "Doublet class")
  save_plot(fig, stem, width = 7, height = 7)
}

#' UMAP coloured by mitochondrial percentage (continuous feature).
plot_umap_mito <- function(srt, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out_dir, "umap_percent_mt")
  fig  <- plot_umap_feature(srt, feature = "percent.mt", reduction = "umap.nn")
  save_plot(fig, stem, width = 7, height = 7)
}

#' PCA elbow plot with a vertical line at n_dims.
#' Uses eval/list2env to avoid capturing the Seurat object in the ggplot environment.
plot_elbow <- function(srt, ref_assay, n_dims, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem  <- file.path(out_dir, "elbow_plot")
  edata <- Seurat::ElbowPlot(srt, ndims = 50,
                              reduction = paste0("pca.", ref_assay))$data
  args  <- list(edata = edata, n_dims = n_dims)
  fig   <- eval(quote({
    ggplot2::ggplot(edata, ggplot2::aes(x = dims, y = stdev)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = n_dims, linetype = "dashed",
                          colour = "firebrick") +
      ggplot2::labs(title    = "PCA elbow plot (SCT)",
                    subtitle = paste0("Dashed line = ", n_dims, " dims used"),
                    x = "PC", y = "Std Dev") +
      cowplot::theme_cowplot()
  }), envir = list2env(args, parent = getNamespace("ggplot2")))
  save_plot(fig, stem, width = 7, height = 5)
}

#' Scatter plot of PC1 vs PC2 coloured by sample.
plot_pca <- function(srt, ref_assay, subset, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem    <- file.path(out_dir, "pca_pc1_pc2")
  pca_emb <- as.data.frame(
    Seurat::Embeddings(srt, paste0("pca.", ref_assay))[, 1:2])
  colnames(pca_emb) <- c("PC1", "PC2")
  pca_emb$sample    <- srt$sample.name
  fig <- ggplot2::ggplot(pca_emb, ggplot2::aes(x = PC1, y = PC2, colour = sample)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.5) +
    ggplot2::labs(title = paste0("PCA — ", subset), colour = "Sample") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4))) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 7))
  save_plot(fig, stem, width = 7, height = 7)
}

#' Compute average expression and % expressing for canonical markers per cluster,
#' write to CSV, and return the path.
write_canonical_markers_csv <- function(srt, ref_assay, marker_res,
                                        annotation_config, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path      <- file.path(out_dir, "canonical_markers.csv")
  clust_col <- grep(paste0("sNN_pca\\..*_res\\.", gsub("\\.", "\\\\.", marker_res), "$"),
                    colnames(srt@meta.data), value = TRUE)[1]
  canon_genes <- intersect(
    unlist(annotation_config$cell_types, use.names = FALSE),
    rownames(SeuratObject::LayerData(srt, layer = "data"))
  )
  avg_exp  <- Seurat::AverageExpression(srt,
                                        assay    = ref_assay,
                                        layer    = "data",
                                        group.by = clust_col,
                                        features = canon_genes)[[ref_assay]]
  clusters <- srt@meta.data[[clust_col]]
  lvls     <- levels(factor(clusters))
  pct_exp  <- do.call(cbind, lapply(lvls, function(cl) {
    mat <- SeuratObject::LayerData(srt, layer = "data")[
      canon_genes, which(clusters == cl), drop = FALSE]
    Matrix::rowMeans(mat > 0)
  }))
  colnames(pct_exp) <- lvls
  colnames(avg_exp) <- lvls
  avg_long <- reshape2::melt(as.matrix(avg_exp),
                             varnames   = c("gene", "cluster"),
                             value.name = "avg_exp")
  pct_long <- reshape2::melt(pct_exp,
                             varnames   = c("gene", "cluster"),
                             value.name = "pct_exp")
  result         <- merge(avg_long, pct_long, by = c("gene", "cluster"))
  result$gene    <- as.character(result$gene)
  result$cluster <- as.character(result$cluster)
  write.csv(result, path, row.names = FALSE)
  path
}

#' Write cluster size summary (long format, all resolutions) to CSV.
save_cluster_size_csv <- function(srt, out_dir) {
  sizes   <- summarise_cluster_sizes(srt)
  long_df <- do.call(rbind, lapply(names(sizes), function(res) {
    df <- sizes[[res]]; df$resolution <- res; df
  }))
  save_csv(long_df, out_dir, "cluster_size_summary.csv")
}

#' Canonical marker bubble plot saved to file; auto-sizes by gene/cluster count.
plot_bubble <- function(canonical_markers_path, annotation_config, subset, out_dir) {
  gene_order <- unlist(annotation_config$cell_types, use.names = FALSE)
  n_genes    <- length(unique(gene_order))
  n_clust    <- length(unique(read.csv(canonical_markers_path)$cluster))
  stem       <- file.path(out_dir, "canonical_markers_bubble")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fig <- plot_canonical_bubble(canonical_markers_path, gene_order,
                               title = paste0("Canonical markers — ", subset))
  save_plot(fig, stem,
            width  = max(6, n_clust * 0.5 + 3),
            height = max(5, n_genes * 0.25 + 2))
}

#' Per-cluster differential expression between sample groups using DElegate.
#'
#' For each cluster, runs DElegate::FindAllMarkers2 comparing groups defined by
#' `group_column` (default: "sample.name"). Returns a combined data.frame with a
#' `cluster` column prepended.
#'
#' @param srt          Seurat object with Idents set to the cluster column.
#' @param excluded_genes Character vector of genes to exclude (e.g. mt + reporters).
#' @param group_column  Metadata column defining the groups to compare.
find_delegate_markers <- function(srt, excluded_genes, group_column = "sample.name") {
  clust_ids   <- levels(Seurat::Idents(srt))
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]

  l <- lapply(clust_ids, function(clust_x) {
    cells <- colnames(srt)[as.character(Seurat::Idents(srt)) == clust_x]
    DElegate::FindAllMarkers2(
      object       = counts_full[, cells, drop = FALSE],
      meta_data    = srt@meta.data[cells, , drop = FALSE],
      group_column = group_column,
      min_fc       = -Inf,
      min_rate     = 0
    )
  })
  names(l) <- clust_ids
  dplyr::bind_rows(l, .id = "cluster")
}

find_delegate_markers_rest_clusters <- function(srt, excluded_genes, group_column = "sample.name") {
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]

  DElegate::FindAllMarkers2(
    object       = counts_full,
    meta_data    = srt@meta.data,
    group_column = group_column,
    min_fc       = -Inf,
    min_rate     = 0
  )
}

find_delegate_markers_all_clusters <- function(srt, excluded_genes, group_column = "sample.name") {
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]

  DElegate::findDE(
    object       = counts_full,
    meta_data    = srt@meta.data,
    group_column = group_column,
    compare = "all_vs_all"
  )
}

find_delegate_markers_all_clusters2 <- function(srt, excluded_genes, group_column = "sample.name") {
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]
  group_levels <- srt@meta.data[[group_column]] %>% as.character() %>% unique()
  lapply(group_levels, function(group_i){
    cat(group_i,"\n")
    DElegate::findDE(
      object       = counts_full,
      meta_data    = srt@meta.data,
      group_column = group_column,
      compare_is_ref = TRUE,
      compare = group_i)
    }) %>%
    bind_rows()

}

#' Per-cluster all-vs-all pairwise DE between sample groups using DElegate.
#'
#' Like \code{find_delegate_markers} but uses \code{DElegate::findDE} with
#' \code{compare = "all_vs_all"} so every pair of groups is tested (e.g.
#' WT vs HET+, WT vs HET-, HET+ vs HET-). Returns a combined data.frame with a
#' \code{cluster} column prepended.
#'
#' @param srt          Seurat object with Idents set to the cluster column.
#' @param excluded_genes Character vector of genes to exclude (e.g. mt + reporters).
#' @param group_column  Metadata column defining the groups to compare.
find_delegate_markers_all_vs_all <- function(srt, excluded_genes, group_column = "sample.name") {
  clust_ids   <- levels(Seurat::Idents(srt))
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]

  l <- lapply(clust_ids, function(clust_x) {
    cells <- colnames(srt)[as.character(Seurat::Idents(srt)) == clust_x]
    DElegate::findDE(
      object       = counts_full[, cells, drop = FALSE],
      meta_data    = srt@meta.data[cells, , drop = FALSE],
      group_column = group_column,
      compare      = "all_vs_all"
    )
  })
  names(l) <- clust_ids
  dplyr::bind_rows(l, .id = "cluster")
}

# ── Pairwise cluster-wise DElegate DE ────────────────────────────────────────

#' Run DElegate pairwise DE within each cluster across multiple contrasts and
#' clustering levels.
#'
#' For each combination of clustering level and contrast, cells are split by
#' cluster, a minimum-cell guard is applied per group, and
#' \code{DElegate::findDE} is called. Results are returned as a nested list
#' \code{result[[clust_level_name]][[contrast_name]]} — each element is a
#' data.frame with the standard DElegate columns plus a \code{cluster} column.
#'
#' @param srt           Seurat object (SCT assay, SCT@counts used).
#' @param excluded_genes Character vector of genes to drop before testing.
#' @param comparisons   Named list of length-2 character vectors:
#'                      \code{list(name = c("groupA", "groupB"), ...)}
#'                      passed directly to \code{DElegate::findDE(compare=)}.
#' @param cluster_levels Named list of named character vectors mapping cell
#'                      barcodes → cluster labels (one per clustering level).
#' @param group_column  Metadata column that holds the group labels used in
#'                      \code{comparisons} (default \code{"sample.name"}).
#' @param min_cells_per_group Minimum cells required per group in a cluster
#'                      before the cluster is tested (default 3).
#' @return Nested list: \code{[[clust_level]][[contrast]]} → data.frame.
find_delegate_pairwise <- function(srt,
                                   excluded_genes,
                                   comparisons,
                                   cluster_levels,
                                   group_column        = "sample.name",
                                   min_cells_per_group = 3L) {
  feat_use    <- setdiff(rownames(srt), excluded_genes)
  counts_full <- srt@assays$SCT@counts[feat_use, , drop = FALSE]
  meta        <- srt@meta.data

  de_one_cluster <- function(comp, cl, cells) {
    group_labels <- meta[cells, group_column]
    n1 <- sum(group_labels == comp[1], na.rm = TRUE)
    n2 <- sum(group_labels == comp[2], na.rm = TRUE)
    if (isTRUE(n1 < min_cells_per_group) || isTRUE(n2 < min_cells_per_group)) {
      message(sprintf("  skip cluster %s [%s vs %s]: %d / %d cells",
                      cl, comp[1], comp[2], n1, n2))
      return(NULL)
    }
    tryCatch(
      DElegate::findDE(
        object       = counts_full[, cells, drop = FALSE],
        meta_data    = meta[cells, , drop = FALSE],
        group_column = group_column,
        compare      = comp
      ) |> dplyr::mutate(cluster = cl),
      error = function(e) {
        message(sprintf("  DElegate error — cluster %s [%s vs %s]: %s",
                        cl, comp[1], comp[2], conditionMessage(e)))
        NULL
      }
    )
  }

  de_one_contrast <- function(comp, clust_vec) {
    clusters <- unique(clust_vec)
    dplyr::bind_rows(lapply(clusters, function(cl) {
      cells <- names(clust_vec)[clust_vec == cl]
      de_one_cluster(comp, cl, cells)
    }))
  }

  lapply(cluster_levels, function(clust_vec) {
    lapply(comparisons, de_one_contrast, clust_vec = clust_vec)
  })
}

#' Per-cluster DE: Het mCh- vs Het mCh+ using HTO tag as replicate.
#'
#' Subsets to mCh- and mCh+ cells (excludes negative and multiplet dm2_class),
#' then runs \code{DElegate::findDE} with \code{group_column = "sample.name"} and
#' \code{replicate_column = "dm2_class"} so each HTO-demultiplexed tag is treated
#' as a biological replicate.
#'
#' @param srt          Seurat object with Idents set to the cluster column.
#' @param excluded_genes Character vector of genes to exclude.
#' @param mch_neg_pattern Regex matching mCh- sample.name values.
#' @param mch_pos_pattern Regex matching mCh+ sample.name values.
find_delegate_markers_mch <- function(srt, excluded_genes,
                                      mch_neg_pattern = "mCh-",
                                      mch_pos_pattern = "mCh\\+") {
  keep_cells <- colnames(srt)[
    grepl(paste0(mch_neg_pattern, "|", mch_pos_pattern), srt@meta.data$sample.name) &
    !srt@meta.data$dm2_class %in% c("negative", "multiplet")
  ]
  if (length(keep_cells) == 0) {
    message("find_delegate_markers_mch: no mCh- / mCh+ cells found — returning empty table.")
    return(dplyr::tibble())
  }
  srt_sub    <- srt[, keep_cells]
  clust_ids  <- levels(Seurat::Idents(srt_sub))
  feat_use   <- setdiff(rownames(srt_sub), excluded_genes)
  counts_full <- srt_sub@assays$SCT@counts[feat_use, , drop = FALSE]

  l <- lapply(clust_ids, function(clust_x) {
    cells <- colnames(srt_sub)[as.character(Seurat::Idents(srt_sub)) == clust_x]
    if (length(unique(srt_sub@meta.data[cells, "sample.name"])) < 2) return(NULL)
    DElegate::findDE(
      object           = counts_full[, cells, drop = FALSE],
      meta_data        = srt_sub@meta.data[cells, , drop = FALSE],
      group_column     = "sample.name",
      replicate_column = "dm2_class",
      compare          = "all_vs_all"
    )
  })
  names(l) <- clust_ids
  dplyr::bind_rows(l[!sapply(l, is.null)], .id = "cluster")
}


#' Feature plots for the top N marker genes per cluster (by avg_log2FC).
plot_marker_features <- function(srt, cluster_markers, out_dir, top_n = 3, n_cols = 4) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  top_genes <- cluster_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = top_n) |>
    dplyr::ungroup() |>
    dplyr::pull(gene) |>
    unique()
  top_genes <- intersect(top_genes,
                         rownames(SeuratObject::LayerData(srt, layer = "data")))
  plots <- lapply(top_genes, function(g) {
    p        <- Seurat::FeaturePlot(srt, features = g,
                                    reduction = "umap.nn", pt.size = 0.5, order = TRUE)
    col_name <- names(p$data)[sapply(p$data, is.numeric) &
                                !names(p$data) %in% names(p$data)[1:2]][1]
    lean_dimplot(p, colour_col = col_name, legend_title = g, pt_size = 0.5) +
      ggplot2::scale_colour_viridis_c(option = "viridis", na.value = "grey90") +
      ggplot2::labs(title = g)
  })
  stem   <- file.path(out_dir, "marker_feature_plots")
  n_rows <- ceiling(length(top_genes) / n_cols)
  save_plot(patchwork::wrap_plots(plots, ncol = n_cols),
            stem, width = n_cols * 3.5, height = n_rows * 3)
}

# ---------------------------------------------------------------
# Cluster annotation helpers
# ---------------------------------------------------------------

#' Re-cluster a processed Seurat object at multiple resolutions.
#'
#' Uses the existing SNN graph (first matching `^sNN_pca`). Returns a data.frame
#' with one row per cell and one column per resolution (named `res_0_1`, etc.),
#' plus a `cell_id` column.
#'
#' @param srt         Seurat object post-\code{process_all}.
#' @param resolutions Numeric vector of resolutions to test.
#' @return data.frame: cell_id + one column per resolution.
recluster_resolutions_df <- function(srt, resolutions = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  snn_graph <- grep("^sNN_pca", SeuratObject::Graphs(srt), value = TRUE)[1]
  if (is.na(snn_graph))
    stop("No SNN graph found. Expected pattern: sNN_pca.*")
  out <- data.frame(cell_id = colnames(srt), stringsAsFactors = FALSE)
  for (r in resolutions) {
    col <- paste0("res_", gsub("\\.", "_", r))
    srt <- Seurat::FindClusters(srt, graph.name = snn_graph, resolution = r, verbose = FALSE)
    out[[col]] <- as.character(srt$seurat_clusters)
  }
  out
}

# ---------------------------------------------------------------
# Reference datasets for label transfer
# ---------------------------------------------------------------

#' Process a raw census reference Seurat object with SCT normalisation for
#' use with build_azimuth_ref().
#'
#' Runs SCTransform → PCA (pca.sct) → UMAP with return.model = TRUE (umap.sct).
#'
#' @param srt        Seurat object with raw counts in RNA assay (from fetch_ref_census).
#' @param n_features Variable features for SCTransform (default 3000).
#' @param npcs       PCA dimensions (default 30).
#' @param n_neighbors UMAP n.neighbors (default 30).
#' @return Seurat object with SCT assay and pca.sct / umap.sct reductions.
seurat_ref_process <- function(srt, n_features = 3000, npcs = 50, n_neighbors = 30) {
  gene_keep <- Matrix::rowSums(Seurat::GetAssayData(srt, layer = "counts")) > 0
  cell_keep <- Matrix::colSums(Seurat::GetAssayData(srt, layer = "counts")) > 0
  srt <- srt[gene_keep, cell_keep]

  # ── SCT ────────────────────────────────────────────────────────
  # SCTransform uses future internally; raise the global size limit so the
  # large Seurat object doesn't trip the 500 MiB default threshold.
  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 8 * 1024^3)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)

  srt <- Seurat::SCTransform(srt, vst.flavor = "v2", variable.features.n = n_features, new.assay.name = "refAssay",
                              verbose = FALSE)
  srt <- Seurat::RunPCA(srt, npcs = npcs, reduction.name = "refDR", verbose = FALSE)
  srt <- Seurat::RunUMAP(srt, dims = 1:npcs, n.neighbors = n_neighbors, reduction = "refDR", reduction.name = "refUMAP",
                          return.model = TRUE, verbose = FALSE)

  Seurat::DefaultAssay(srt) <- "refAssay"
  srt
}

#' Download a Tabula Muris Senis (10x) dataset from cellxgene.census as raw counts.
#'
#' Uses the high-level cellxgene.census::get_seurat() API. Census version pinned
#' to 2025-11-08. get_seurat() returns Ensembl IDs as rownames; this function
#' remaps them to gene symbols (make.unique to handle any duplicates) so that
#' downstream SCTransform and gene pattern matching work correctly.
#'
#' @param dataset_id  cellxgene dataset UUID.
#' @return Unprocessed Seurat object with raw counts in the RNA assay and
#'   gene symbols as rownames.
fetch_ref_census <- function(dataset_id) {
  census <- cellxgene.census::open_soma(census_version = "2025-11-08")
  on.exit(census$close(), add = TRUE)
  srt <- cellxgene.census::get_seurat(
    census           = census,
    organism         = "Mus musculus",
    obs_value_filter = paste0("dataset_id == '", dataset_id, "'"),
    obs_column_names = c("cell_type", "tissue", "sex", "development_stage")
  )
  # get_seurat() sets rownames to Ensembl IDs; remap to gene symbols.
  # feature_name is stored in the RNA assay's meta.features after get_seurat().
  # rownames(srt) updates both the assay internals and the Seurat object-level
  # feature index — rownames(srt[["RNA"]]) alone would leave a mismatch.
  symbols <- make.unique(as.character(srt[["RNA"]]@meta.features[["feature_name"]]))
  rownames(srt) <- symbols
  srt
}


#' Download the Azimuth Mouse PanSci reference from Zenodo (doi:10.5281/zenodo.14042645).
#'
#' The RDS is a processed Azimuth reference Seurat object. This function:
#'   1. Downloads the file to a session temp dir.
#'   2. Extracts the count matrix and cell type labels into a plain Seurat object
#'      so that seurat_ref_process() can normalise it consistently with other refs.
#'
#' The reference uses the column "celltype.l2" for fine-grained labels and
#' "celltype.l1" for broad labels. We expose "celltype.l2" as "cell_type".
#'
#' @return Unprocessed Seurat object with raw counts in RNA assay and cell_type column.
fetch_ref_pansci <- function() {
  url <- "https://zenodo.org/records/14042645/files/pansci_filtered.rds?download=1"
  tmp <- file.path(tempdir(), "pansci_filtered.rds")
  if (!file.exists(tmp)) {
    message("Downloading Azimuth Mouse PanSci reference (~1.4 GB) ...")
    options(timeout = 3600)
    download.file(url, tmp, mode = "wb", quiet = FALSE)
  }
  message("Loading PanSci RDS ...")
  ref <- readRDS(tmp)

  # Azimuth references may be plain Seurat or wrapped in a list.
  if (is.list(ref) && !inherits(ref, "Seurat")) {
    if ("object" %in% names(ref)) ref <- ref$object else ref <- ref[[1]]
  }

  message("PanSci meta columns: ", paste(colnames(ref@meta.data), collapse = ", "))
  message("PanSci assays: ",       paste(names(ref@assays),       collapse = ", "))
  message("PanSci reductions: ",   paste(names(ref@reductions),   collapse = ", "))
  message("PanSci dims: ", nrow(ref), " genes x ", ncol(ref), " cells")

  # Cell type labels — broad search across common naming conventions
  ct_col <- intersect(
    c("Sub_cell_type", "Main_cell_type",
      "celltype.l2", "celltype.l1", "celltype", "cell_type",
      "annotation", "cluster", "seurat_clusters", "label"),
    colnames(ref@meta.data)
  )[1]
  if (is.na(ct_col))
    stop("PanSci reference has no recognised cell-type column. Columns: ",
         paste(colnames(ref@meta.data), collapse = ", "))

  # Subset to bone marrow cells if Organ_name column exists
  if ("Organ_name" %in% colnames(ref@meta.data)) {
    bm_cells <- colnames(ref)[grepl("(?i)bone.marrow|marrow|BoneMarrow",
                                    ref$Organ_name, perl = TRUE)]
    if (length(bm_cells) > 0) {
      message("Subsetting to bone marrow: ", length(bm_cells), " / ", ncol(ref), " cells")
      ref <- ref[, bm_cells]
    } else {
      message("No bone marrow cells found in Organ_name — using all cells")
    }
  }

  # Convert Ensembl IDs to gene symbols if rownames look like Ensembl IDs
  if (any(grepl("^ENSMUSG", rownames(ref)))) {
    message("Converting Ensembl IDs to gene symbols via biomaRt ...")
    ids <- sub("\\..*", "", rownames(ref))  # strip version suffix
    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    map  <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol"),
      filters    = "ensembl_gene_id",
      values     = ids,
      mart       = mart
    )
    sym <- map$mgi_symbol[match(ids, map$ensembl_gene_id)]
    sym[is.na(sym) | sym == ""] <- ids[is.na(sym) | sym == ""]
    rownames(ref) <- make.unique(sym)
  }

  # Build a plain Seurat with raw counts so seurat_ref_process() can normalise it
  counts_mat <- tryCatch(
    Seurat::GetAssayData(ref, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(ref, assay = "RNA", slot = "counts")
  )
  meta <- data.frame(cell_type = as.character(ref@meta.data[[ct_col]]),
                     row.names = colnames(ref))
  srt <- Seurat::CreateSeuratObject(counts = counts_mat, meta.data = meta)
  rm(ref, counts_mat); gc()
  srt
}

#' Load the Mouse Bone Marrow Cell Atlas Azimuth reference (syn53222728).
#'
#' This is a pre-built SCT reference (refAssay = SCTAssay, data layer only).
#' Reductions: refDR (50-dim PCA-equivalent), refUMAP (with model for ProjectUMAP).
#' Cell type labels: Level1 (broad) and Level2 (fine-grained, used as cell_type).
#'
#' Unlike census/PanSci refs, this object is NOT passed through seurat_ref_process().
#' It is used directly by transfer_labels_df_azimuth() which sets the correct
#' reference.reduction and reduction.model.
#'
#' @param path Path to syn53222728.Rds on disk.
#' @return The reference Seurat object as-is (with cell_type = Level2 added).
load_ref_marrow_atlas <- function(path) {
  ref <- readRDS(path)
  # Expose Level2 as cell_type for compatibility with transfer functions
  ref$cell_type <- as.character(ref$Level2)
  # Fix key/colname mismatch in refDR: feature.loadings columns are named "PC_*"
  # but the reduction key is "refdr_". Seurat v5 validObject() requires them to match.
  if ("refDR" %in% names(ref@reductions)) {
    fl <- ref@reductions[["refDR"]]@feature.loadings
    if (ncol(fl) > 0 && !startsWith(colnames(fl)[1], "refdr_")) {
      colnames(ref@reductions[["refDR"]]@feature.loadings) <-
        paste0("refdr_", seq_len(ncol(fl)))
      colnames(ref@reductions[["refDR"]]@cell.embeddings) <-
        paste0("refdr_", seq_len(ncol(ref@reductions[["refDR"]]@cell.embeddings)))
    }
  }
  ref
}

#' Download the mouse bone marrow niche reference (Seurat v5 RDS) from OSF.
#'
#' Source: "Resolving the Bone Marrow Niche Heterogeneity" (doi:10.3389/fcell.2021.622519)
#' OSF project: https://osf.io/ne9vj/files/wnxc7 — 32,743 BM niche cells, Seurat v5.
#' Direct download URL: https://osf.io/download/wnxc7/
#'
#' @param path Destination file path (e.g. "data/bone-marrow-seurat-v5.rds").
#' @return \code{path} (invisibly), for use as a \code{format = "file"} target.
fetch_ref_bm_niche <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    message("Downloading BM niche reference from OSF (~309 MB)...")
    utils::download.file(
      url     = "https://osf.io/download/wnxc7/",
      destfile = path,
      mode    = "wb",
      method  = "libcurl",
      quiet   = FALSE
    )
  }
  invisible(path)
}

#' Load the mouse bone marrow niche reference for label transfer.
#'
#' Exposes \code{cell_type} metadata column (= \code{CellType} or the first
#' character metadata column that looks like a cell-type label).
#' Also normalises the SCT assay if present so it is compatible with
#' \code{transfer_labels_df_azimuth()}.
#'
#' @param path Path to the downloaded RDS file.
#' @return Seurat reference object with \code{cell_type} column.
load_ref_bm_niche <- function(path) {
  ref <- readRDS(path)
  # Expose cell type label as generic cell_type column
  ct_candidates <- c("CellType", "cell_type", "celltype", "annotation",
                     "label", "cell.type", "Celltype")
  ct_col <- ct_candidates[ct_candidates %in% colnames(ref@meta.data)][1]
  if (!is.na(ct_col)) {
    ref$cell_type <- as.character(ref@meta.data[[ct_col]])
  } else {
    # Fall back to first character column with >1 unique value
    char_cols <- vapply(ref@meta.data, is.character, logical(1))
    cat_cols  <- names(ref@meta.data)[char_cols]
    for (col in cat_cols) {
      n_unique <- length(unique(ref@meta.data[[col]]))
      if (n_unique > 1 && n_unique < nrow(ref@meta.data) * 0.1) {
        ref$cell_type <- as.character(ref@meta.data[[col]])
        message("BM niche: using '", col, "' as cell_type (", n_unique, " values)")
        break
      }
    }
  }
  ref
}

#' Build a Seurat reference from any Mouse Cell Atlas (MCA) DGE dataset.
#'
#' Source: Mouse Cell Atlas (MCA) — https://bis.zju.edu.cn/MCA/
#' All MCA zip files share the same file layout:
#'   <prefix>_dge.csv          — genes × cells raw count matrix (row indices as gene IDs)
#'   <prefix>_gene.csv         — gene name lookup (row index → symbol)
#'   <prefix>_barcodes_anno.csv — cell barcode + Idents.pbmc. cell type label
#'
#' @param dir    Directory containing the extracted MCA files.
#' @param prefix Filename prefix, e.g. "Peripheral-Blood", "Bone-Marrow", etc.
#' @param label  Short label used in log messages, e.g. "MCA PBMC".
#' @return Seurat object (DietSeurat-slimmed) ready for \code{transfer_labels_df}.
load_ref_mca <- function(dir, prefix, label = prefix) {

  message("Building Seurat object (", nrow(mat), " genes × ", ncol(mat), " cells)...")
  srt <- Read10X(dir) |> CreateSeuratObject()
  meta <- read.table(glue::glue("{dir}/metadata.tsv.gz"),
                     sep="\t",
                     header=TRUE,
                     row.names=1)
  srt <- AddMetaData(srt, meta)

  message("SCTransform + PCA + UMAP...")
  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 8 * 1024^3)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)
  srt <- Seurat::SCTransform(srt, vst.flavor = "v2", variable.features.n = 3000,
                              verbose = FALSE,new.assay.name = "refAssay")
  srt <- Seurat::RunPCA(srt, npcs = 50, reduction.name = "refDR", verbose = FALSE)
  srt <- Seurat::RunUMAP(srt, dims = 1:50, reduction = "refDR", reduction.name = "refUMAP",
                          return.model = TRUE, verbose = FALSE)
  srt <- Seurat::DietSeurat(srt, assays = "refAssay", layers = c("counts", "data", "scale.data"),
                            dimreducs = c("refDR", "refUMAP"))
  gc()
  message(label, " reference ready: ", length(unique(srt$cell_type)), " cell types")
  srt
}

load_ref_mca_pbmc        <- function(dir) load_ref_mca(dir, "Peripheral-Blood",        "MCA PBMC")
load_ref_mca_bm          <- function(dir) load_ref_mca(dir, "Bone-Marrow",             "MCA Bone Marrow")
load_ref_mca_bm_ckit     <- function(dir) load_ref_mca(dir, "Bone-Marrow-c-kit",       "MCA BM c-Kit")
load_ref_mca_bm_mesench  <- function(dir) load_ref_mca(dir, "Bone-Marrow-Mesenchyme",  "MCA BM Mesenchyme")

#' Build a Seurat reference from SCP1248 (Baccin et al. 2019 BM niche atlas).
#'
#' SCP1248 files (from the Broad Single Cell Portal):
#'   counts.normalized.csv.gz — genes × cells log-normalised matrix (rows = genes)
#'   metadata.csv             — cell metadata with UMAP coords + Harmonized Label
#'
#' Constructs a Seurat object with:
#'   - RNA assay containing the pre-normalised counts in the \code{data} layer
#'   - UMAP reduction from the X/Y coordinates in metadata
#'   - \code{cell_type} = Harmonized Label
#'   - \code{source_dataset}, \code{source_label} additional metadata columns
#'
#' @param dir Path to the SCP1248 directory containing the two files above.
#' @return Seurat object ready for \code{transfer_labels_df_azimuth()}.
load_ref_scp1248 <- function(dir) {
  counts_path   <- file.path(dir, "counts.normalized.csv.gz")
  metadata_path <- file.path(dir, "metadata.csv")

  message("Reading SCP1248 normalized counts...")
  mat_raw <- data.table::fread(counts_path, data.table = FALSE)
  rownames(mat_raw) <- mat_raw[[1]]
  mat_raw <- mat_raw[, -1, drop = FALSE]
  # Convert to sparse immediately to save memory
  mat <- Matrix::Matrix(as.matrix(mat_raw), sparse = TRUE)
  rm(mat_raw); gc()

  message("Reading SCP1248 metadata...")
  meta <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- meta[meta$NAME != "TYPE", ]
  rownames(meta) <- meta$NAME

  cells <- intersect(colnames(mat), rownames(meta))
  mat   <- mat[, cells, drop = FALSE]
  meta  <- meta[cells, , drop = FALSE]

  message("Building Seurat object (", nrow(mat), " genes × ", ncol(mat), " cells)...")
  # Pass the normalised sparse matrix as counts; CreateSeuratObject stores it
  # in the counts layer. We then copy it to data and treat it as pre-normalised.
  srt <- SeuratObject::CreateSeuratObject(counts = mat, meta.data = meta)
  srt <- SeuratObject::SetAssayData(srt, assay = "RNA", layer = "data", new.data = mat)
  rm(mat); gc()

  colnames(srt@meta.data)[colnames(srt@meta.data) == "Harmonized Label"] <- "cell_type"
  colnames(srt@meta.data)[colnames(srt@meta.data) == "Source Dataset"]   <- "source_dataset"
  colnames(srt@meta.data)[colnames(srt@meta.data) == "Source Label"]     <- "source_label"

  # Add UMAP from X/Y coordinates in metadata
  umap_coords <- matrix(as.numeric(c(meta$X, meta$Y)), ncol = 2,
                        dimnames = list(cells, c("UMAP_1", "UMAP_2")))
  srt[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap_coords, key = "UMAP_", assay = "RNA")

  # PCA for FindTransferAnchors — data layer is already normalised, skip NormalizeData.
  message("Running PCA for label transfer...")
  srt <- Seurat::FindVariableFeatures(srt, nfeatures = 3000, verbose = FALSE)
  var_genes <- Seurat::VariableFeatures(srt)
  srt <- Seurat::ScaleData(srt, features = var_genes, verbose = FALSE)
  srt <- Seurat::RunPCA(srt, features = var_genes, npcs = 30, verbose = FALSE)

  # Slim down: drop counts + scale.data, keep only data layer + PCA + UMAP + metadata
  srt <- Seurat::DietSeurat(srt, assays = "RNA", layers = "data",
                            dimreducs = c("pca", "umap"))
  gc()

  message("SCP1248 reference ready: ", length(unique(srt$cell_type)), " cell types")
  srt
}

#' Load PanSci reference from RDS and prepare for Azimuth.
#'
#' Loads a pre-filtered Seurat object from PanSci (stromal or immune, WT only)
#' and processes it with SCTransform + PCA + UMAP for use as an Azimuth reference.
#'
#' The input RDS must contain:
#'   - \code{Sub_cell_type_organ} column for cell type labels
#'
#' Processing steps:
#'   - Extract cell_type from Sub_cell_type_organ
#'   - SCTransform (v2, 3000 features) → refAssay
#'   - PCA (50 dims) → refDR
#'   - UMAP (dims 1:50, return.model = TRUE) → refUMAP
#'   - DietSeurat to retain only refAssay layers + refDR + refUMAP
#'
#' @param path Path to pansci RDS file (e.g., pansci_stromal_wt.rds or pansci_immune_wt.rds).
#' @param label Label for log messages (e.g., "PanSci Stromal WT" or "PanSci Immune WT").
#' @return Seurat object ready for Azimuth (refAssay/refDR/refUMAP).
load_ref_pansci <- function(path, label = "PanSci") {
  message("Loading ", label, " reference from ", path, "...")
  srt <- readRDS(path)

  message("SCTransform + PCA + UMAP...")
  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 8 * 1024^3)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)

  srt <- Seurat::SCTransform(srt, vst.flavor = "v2", variable.features.n = 3000,
                              verbose = FALSE, new.assay.name = "refAssay")
  srt <- Seurat::RunPCA(srt, npcs = 50, reduction.name = "refDR", verbose = FALSE)
  srt <- Seurat::RunUMAP(srt, dims = 1:50, reduction = "refDR", reduction.name = "refUMAP",
                         return.model = TRUE, verbose = FALSE)
  srt <- Seurat::DietSeurat(srt, assays = "refAssay", layers = c("counts", "data", "scale.data"),
                            dimreducs = c("refDR", "refUMAP"))
  gc()
  message(label, " reference ready: ", length(unique(srt$cell_type)), " cell types")
  srt
}

#' Project query cells onto a reference and transfer labels.
#'
#' Uses FindTransferAnchors (pcaproject onto reference PCA) anchored on the
#' intersection of variable features from both query and reference — focusing
#' anchors on informative genes and improving anchor counts vs. all-gene transfer.
#'
#' MappingScore is computed before TransferData (so it is not overwritten) and
#' stored as \code{ctr.mapping.score} in query metadata. TransferData adds
#' \code{predicted.cell_type} and \code{predicted.cell_type.score}.
#'
#' Output coordinates (\code{refumap_1/2}) are the query cells projected onto
#' the reference UMAP via \code{ProjectUMAP}, stored under the same reduction
#' name as \code{ref_umap} in the query object.
#'
#' All assay and reduction names must be supplied explicitly.
#'
#' Returns a tidy data.frame (one row per query cell):
#'   cell_id, predicted.cell_type, predicted.cell_type.score,
#'   mapping_score, refumap_1, refumap_2
#'
#' @param query         Seurat query object (must have umap.nn reduction).
#' @param reference     Pre-processed Seurat reference with PCA and UMAP.
#' @param ref_name      Short identifier used in column names and messages.
#' @param norm_method   "SCT" or "LogNormalize" — must match both assays.
#' @param ref_assay     Assay name in the reference (e.g. "refAssay", "RNA").
#' @param query_assay   Assay name in the query (e.g. "SCT", "RNA").
#' @param ref_reduction PCA reduction name in the reference (e.g. "refDR", "pca").
#' @param ref_umap      UMAP model reduction name in the reference (e.g. "refUMAP", "umap").
#' @param dims          PCA dimensions to use (e.g. 1:50, 1:30).
#' @return data.frame (one row per query cell).
transfer_labels_df_azimuth <- function(query,
                                       reference,
                                       ref_name,
                                       norm_method   = NULL,
                                       ref_assay     = NULL,
                                       query_assay = NULL,
                                       ref_reduction = NULL,
                                       ref_umap      = NULL,
                                       dims          = NULL) {
  # lowest number of computed dimensions in eithe query or reference
  ndim <- ncol(Seurat::Embeddings(reference, reduction = ref_reduction))
  dims <- 1:min(ndim, 50)

  # set default assay
  Seurat::DefaultAssay(reference) <- ref_assay
  Seurat::DefaultAssay(query) <- query_assay

  # Explicit feature intersection
  query_var_feats <- VariableFeatures(query)
  ref_var_feats <- FindVariableFeatures(reference) %>%
    VariableFeatures()
  features <- intersect(query_var_feats, ref_var_feats)
  message("FindTransferAnchors features: ", length(features), " overlapping genes")
  if (length(features) == 0)
    stop("No overlapping features between reference (", ref_assay, ") and query (", query_assay, ").")

  anchors <- Seurat::FindTransferAnchors(
    reference            = reference,
    query                = query,
    normalization.method = norm_method,
    reference.assay      = ref_assay,
    query.assay          = query_assay,
    features             = features,
    reduction            = "pcaproject",
    reference.reduction  = ref_reduction,
    dims                 = dims,
    mapping.score.k      = 100,
    verbose              = FALSE
  )

  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors,
                            ndim = max(dims)),
    col.name = "ctr.mapping.score"
  )

  query <- TransferData(
    anchorset = anchors,
    refdata = list(cell_type  = "cell_type"),
    reference = reference ,
    query = query,
    verbose = FALSE
  )

  pca_red <- paste0("ref.pca.", ref_name)
  query <- Seurat::IntegrateEmbeddings(
    anchorset          = anchors,
    reference          = reference,
    query              = query,
    reductions         = "pcaproject",
    new.reduction.name = pca_red)

  query <- Seurat::ProjectUMAP(
    query              = query,
    query.reduction    = pca_red,
    reference          = reference,
    reference.reduction = ref_reduction,
    reduction.model    = ref_umap,
    reduction.name = ref_umap,
    verbose            = FALSE)

  out <- Embeddings(query, reduction = ref_umap)[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    left_join(query@meta.data %>% rownames_to_column("cell_id")) %>%
    dplyr::select(cell_id,
                  predicted.cell_type,
                  predicted.cell_type.score,
                  mapping_score = ctr.mapping.score,
                  refumap_1 = matches("_1$"),
                  refumap_2 = matches("_2$"))

  return(out)

}

#' Build and save an Azimuth reference from a processed Seurat object.
#'
#' Wraps AzimuthReference() + SaveAzimuthReference(). The saved folder can then
#' be passed directly to RunAzimuth() / run_azimuth_df().
#'
#' @param srt        Processed Seurat object (SCT assay + PCA + UMAP with return.model=TRUE).
#' @param ref_name   Short identifier, e.g. "tabula_muris_marrow".
#' @param refDR      Name of the PCA reduction in srt (default "pca").
#' @param refUMAP    Name of the UMAP reduction in srt (default "umap").
#' @param refAssay   Name of the SCT assay in srt (default "SCT").
#' @param metadata   Metadata column name for cell type labels (default "cell_type").
#' @param azimuth_dir Base directory where reference folders are saved (default "data/azimuth").
#' @return Path to the saved reference folder (for use as format = "file" target).
build_azimuth_ref <- function(srt, ref_name,
                               refDR      = "pca",
                               refUMAP    = "umap",
                               refAssay   = "SCT",
                               metadata   = "cell_type",
                               azimuth_dir = "data/azimuth") {
  ref_path <- file.path(azimuth_dir, ref_name, "/")
  dims <- length(srt[[refDR]]@stdev)
  col_list <- list(Azimuth::CreateColorMap(srt,unique(srt$cell_type)))
  names(col_list) <- metadata
  dir.create(ref_path, recursive = TRUE, showWarnings = FALSE)
  ref <- Azimuth::AzimuthReference(
    object   = srt,
    refUMAP  = refUMAP,
    plotref  = refUMAP,
    refDR    = refDR,
    metadata = metadata,
    refAssay = refAssay,
    dims = 1:dims,
    colormap = col_list
    )

  Azimuth::ValidateAzimuthReference(ref, ad.name = "AzimuthReference")

  Azimuth::SaveAzimuthReference(ref, folder = ref_path)
  invisible(ref_path)
}

#' Run Azimuth label transfer using a local reference directory and return
#' a data.frame with the same schema as transfer_labels_df_azimuth().
#'
#' RunAzimuth() handles anchor finding, label transfer, mapping score, and
#' ref UMAP projection in one call. This wrapper extracts the results into
#' the standard projection data.frame used by attach_projection() in the notebook.
#'
#' @param query         Seurat query object (SCT-normalised, umap.nn present).
#' @param ref_dir       Path to directory saved by build_azimuth_ref().
#' @param ref_name      Short identifier used in column names.
#' @param annotation.levels  Annotation levels to map (NULL = all available).
#' @return data.frame: cell_id, predicted.cell_type, mapping_score,
#'         predicted.cell_type.score, refumap_1, refumap_2.
#' Rename RNA assay features using a named vector (old_name → new_name).
#'
#' Used for cross-species projection: rename mouse gene names to their human
#' ortholog equivalents before FindTransferAnchors against a human reference.
#' Genes not present in new_feats are dropped.
#'
#' @param srt       Seurat object.
#' @param new_feats Named character vector: names = old gene names, values = new names.
#' @return Seurat object with RNA assay features renamed.
rename_rna_feats <- function(srt, new_feats) {

  #define a new assay with new names and overwrite old RNA assay
  DefaultAssay(srt) <- "SCT"
  mat <- GetAssayData(srt,assay = "RNA",layer = "counts")
  old_feats <- mat@Dimnames[[1]]
  ##limit to genes with human ortho
  feat_filt <- new_feats[names(new_feats) %in% old_feats]
  ##overwrite old RNA assay and rename the features
  filt_mat <- mat[names(feat_filt),]
  rm(mat);gc()
  filt_mat@Dimnames[[1]] <- unname(feat_filt)
  srt[["RNA"]] <- CreateAssay5Object(filt_mat)
  DefaultAssay(srt) <- "RNA"
  srt[["SCT"]] <- NULL
  return(srt)
  }


ms_hs_1_to_1_ortho_vec <- function(mm_hs_ortho,
                                   ref,
                                   query,
                                   assay = "RNA",
                                   layer = "counts"){

  hs_counts <- Matrix::rowSums(GetAssayData(ref,assay = assay, layer = layer) )
  mm_counts <- Matrix::rowSums(GetAssayData(query,assay = assay, layer = layer) )
  rm(list = c("ref", "query")); gc()

  df <- mm_hs_ortho %>%
    filter(nchar(Mouse.gene.name)>0 & nchar(human_ortho)>0) %>%
    left_join(enframe(hs_counts,name = "human_ortho",value = "human_counts")) %>%
    left_join(enframe(mm_counts,name = "Mouse.gene.name", value = "mouse_counts")) %>%
    na.omit() %>%
    group_by(Mouse.gene.name) %>%
    slice_max(order_by = human_counts, n = 1,with_ties = FALSE) %>%
    ungroup() %>%
    group_by(human_ortho) %>%
    slice_max(order_by = mouse_counts, n = 1,with_ties = FALSE) %>%
    ungroup()

  setNames(df$human_ortho, df$Mouse.gene.name)
}


#' Project query cells onto a Seurat reference using anchor-based transfer.
#'
#' Runs FindTransferAnchors → TransferData → IntegrateEmbeddings → ProjectUMAP
#' and returns a flat data.frame matching the run_azimuth_df output format
#' (cell_id, refUMAP_1, refUMAP_2, + all query metadata + predicted labels).
#'
#' @param query            Query Seurat object.
#' @param ref              Reference Seurat object (must have a UMAP model).
#' @param norm_method      Normalisation method (default: "LogNormalize").
#' @param ref_assay        Reference assay (default: "RNA").
#' @param query_assay      Query assay (default: "RNA").
#' @param reduction_method Anchor finding reduction (default: "pcaproject").
#' @param ref_reduction    Reference reduction to use (default: "pca").
#' @param dims             PCA dims to use (default: 1:50).
#' @param mapping.score.k  k for mapping score calculation (default: 50).
#' @param tranfer_list     Named list of reference metadata vectors to transfer.
#' @return data.frame with cell_id, ref.umap coords, and query metadata.
project_seurat_ref <- function(query,
                                ref,
                                norm_method      = "LogNormalize",
                                ref_assay        = "RNA",
                                query_assay      = "RNA",
                                reduction_method = "pcaproject",
                                ref_reduction    = "pca",
                                dims             = 1:50,
                                mapping.score.k  = 50,
                               project_umap,
                                tranfer_list) {
  anchors <- FindTransferAnchors(
    reference         = ref,
    query             = query,
    normalization.method = norm_method,
    reference.assay   = ref_assay,
    query.assay       = query_assay,
    reduction         = reduction_method,
    reference.reduction = ref_reduction,
    dims              = dims,
    mapping.score.k   = mapping.score.k
  )

  query <- AddMetaData(
    object   = query,
    metadata = MappingScore(anchors = anchors, ndim = length(dims)),
    col.name = "mapping.score"
  )

  query <- TransferData(
    anchorset = anchors,
    refdata   = tranfer_list,
    reference = ref,
    query     = query
  )

  query <- IntegrateEmbeddings(
    anchorset      = anchors,
    reference      = ref,
    query          = query,
    reductions     = reduction_method,
    new.reduction.name = "ref.pca"
  )

  if(project_umap){
    query <- ProjectUMAP(
      query            = query,
      query.reduction  = "ref.pca",
      reference        = ref,
      reference.reduction = "pca",
      reduction.model  = "umap",
      reduction.name   = "ref.umap",
      k.param          = 50
    )

    # Build output before rm() — rm() removes local vars which would
    # make query inaccessible if called before Embeddings()
    out <- Embeddings(query, reduction = "ref.umap")[, 1:2] %>%
      as.data.frame() %>%
      rownames_to_column("cell_id") %>%
      left_join(query@meta.data %>% rownames_to_column("cell_id"))

    rm(query, ref, anchors)
    gc()
    return(out)

  }
  # make query inaccessible if called before Embeddings()
  out <- query@meta.data

  rm(query, ref, anchors)
  gc()
  return(out)
}

run_azimuth_df <- function(query,
                           ref_dir) {
  query <- Azimuth::RunAzimuth(
    query             = query,
    reference         = ref_dir,
    verbose           = TRUE
  )

  out <- Embeddings(query, reduction = "ref.umap")[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    left_join(query@meta.data %>% rownames_to_column("cell_id"))

  return(out)
}

project_mm_onto_hs_seurat <- function(ref,
                                      query,
                                      ortho_table,
                                      assay,
                                      tranfer_list,
                                      project_umap = TRUE){

  new_feats <- ms_hs_1_to_1_ortho_vec(ortho_table,
                                      ref,query,assay = "RNA", layer = "counts")
  proj_df <- rename_rna_feats(query,
                            new_feats) %>%
    NormalizeData(assay = assay,
                  verbose = FALSE) %>%
    project_seurat_ref(ref,
                       norm_method = "LogNormalize",
                       ref_assay   = assay,
                       query_assay = assay,
                       reduction_method = "pcaproject",
                       ref_reduction = "pca",
                       dims = 1:50,
                       mapping.score.k=50,
                       tranfer_list = tranfer_list,
                       project_umap = project_umap)

  return(proj_df)
}


#' Project query cells onto a Symphony reference and return predicted labels.
#'
#' Fixes the stale uwot external pointer by resetting save_uwot_path to the
#' correct absolute path before projection. Uses knnPredict for label transfer.
#'
#' @param query      Seurat object (query cells).
#' @param ref_rds    Absolute path to the Symphony reference RDS.
#' @param uwot_path  Absolute path to the matching .uwot model file.
#' @param label_col  Column in ref$meta_data to transfer (default: "CellType_Annotation").
#' @param k          k for kNN label transfer (default: 5).
#' @return data.frame with cell_id, UMAP1, UMAP2, predicted_label, prediction_score.
run_symphony_df <- function(query,
                             ref_rds,
                             uwot_path,
                             label_col = "CellType_Annotation",
                             k         = 5) {
  ref <- readRDS(ref_rds)

  # Fix stale external pointer — Symphony reloads uwot from this path on first use
  ref$save_uwot_path <- uwot_path

  # Variable genes defined in the reference
  vargenes <- ref$vargenes_means_sds$symbol

  # Normalized query expression (log-normalised RNA counts)
  exp_query <- Seurat::GetAssayData(query, assay = "RNA", layer = "data")

  # Keep only genes present in both query and reference
  shared_genes <- intersect(vargenes, rownames(exp_query))
  missing_pct  <- round(100 * (1 - length(shared_genes) / length(vargenes)), 1)
  if (missing_pct > 50)
    warning("run_symphony_df: ", missing_pct, "% of reference var genes missing from query.")

  message("  Mapping ", ncol(exp_query), " query cells onto Symphony reference...")
  result <- symphony::mapQuery(
    exp_query      = exp_query[shared_genes, ],
    metadata_query = query@meta.data,
    ref_obj        = ref,
    vars           = NULL,
    do_umap        = TRUE,
    do_normalize   = FALSE,
    verbose        = FALSE
  )

  # kNN label transfer from reference embeddings to query
  predicted <- symphony::knnPredict(
    query_obj  = result,
    ref_obj    = ref,
    label_col  = label_col,
    k          = k,
    confidence = TRUE
  )

  # Return flat data frame matching run_azimuth_df output format
  umap_df <- as.data.frame(result$umap$embedding)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$cell_id         <- rownames(result$meta_data)
  umap_df$predicted_label <- predicted$cell_type_pred
  umap_df$prediction_score <- predicted$cell_type_pred_score

  umap_df[, c("cell_id", "UMAP1", "UMAP2", "predicted_label", "prediction_score")]
}

# Extract a flat data.frame from a pre-filter (raw) processed Seurat object.
# Includes: umap.nn coords, PC1, seurat_clusters, and all QC metadata columns.
# Also merges in any extra columns from metadata_flagged_all (e.g. pass_* flags)
# that were lost during process_all().
#
# @param srt_raw     mpn_raw_processed_* Seurat object.
# @param flagged_all metadata_flagged_all data.frame (rows = cells, rownames = barcodes).
# @return data.frame with cell_id, UMAP1, UMAP2, PC1, seurat_clusters, and all metadata.
extract_raw_meta <- function(srt_raw, flagged_all = NULL,
                             umap_reduction = "umap.nn",
                             pca_reduction  = "pca.SCT") {
  emb <- as.data.frame(Seurat::Embeddings(srt_raw, umap_reduction))
  colnames(emb) <- c("UMAP1", "UMAP2")

  # First PC — used as a proxy for cell complexity on UMAP
  pca_key <- grep(pca_reduction, Seurat::Reductions(srt_raw), value = TRUE)[1]
  if (!is.na(pca_key)) {
    pc1 <- Seurat::Embeddings(srt_raw, pca_key)[, 1, drop = TRUE]
    emb$PC1 <- pc1[rownames(emb)]
  }

  meta_want <- c("seurat_clusters", "sample.name", "dm2_class", "Phase",
                 "scDblFinder.score", "scDblFinder.class",
                 "nCount_RNA", "nFeature_RNA",
                 "percent.mt", "percent.rp", "percent.hb",
                 "percent.stress", "count_feat_ratio",
                 "pass_cell", "pass_mt", "pass_hb")
  meta_cols <- intersect(meta_want, colnames(srt_raw@meta.data))
  out <- cbind(data.frame(cell_id = rownames(emb)), emb,
               srt_raw@meta.data[rownames(emb), meta_cols, drop = FALSE])

  # Merge in extra columns from metadata_flagged_all (pass_* flags, etc.)
  if (!is.null(flagged_all)) {
    new_cols <- setdiff(colnames(flagged_all), colnames(out))
    shared   <- intersect(rownames(flagged_all), out$cell_id)
    if (length(new_cols) > 0 && length(shared) > 0) {
      extra <- flagged_all[shared, new_cols, drop = FALSE]
      extra$cell_id <- rownames(extra)
      out <- merge(out, extra, by = "cell_id", all.x = TRUE)
    }
  }

  rownames(out) <- NULL
  out
}

# Extract a flat data.frame of query metadata needed for annotation notebooks.
# Includes: umap.nn coords, all res_*/sNN_ cluster columns, and all useful
# metadata columns (QC metrics, cell cycle, doublet scores, HTO, sample info).
# No Seurat object needed in the notebook after this.
extract_anno_meta <- function(srt, umap_reduction = "umap.nn") {
  emb <- as.data.frame(Seurat::Embeddings(srt, umap_reduction))
  colnames(emb) <- c("UMAP1", "UMAP2")
  # cluster columns (main + recluster resolutions)
  res_cols <- grep("^res_|^sNN_", colnames(srt@meta.data), value = TRUE)
  # fixed metadata columns present in processed Seurat objects
  fixed_cols <- c("sample.name", "dm2_class", "genotype",
                  "Phase", "S.Score", "G2M.Score",
                  "scDblFinder.score", "scDblFinder.class",
                  "percent.mt", "percent.rp", "percent.hb",
                  "percent.stress", "nCount_RNA", "nFeature_RNA",
                  "nCount_SCT", "nFeature_SCT", "count_feat_ratio")
  meta_cols <- intersect(fixed_cols, colnames(srt@meta.data))
  cbind(
    data.frame(cell_id = rownames(emb)),
    emb,
    srt@meta.data[rownames(emb), c(res_cols, meta_cols), drop = FALSE]
  )
}

# Extract SCT normalised expression for a set of genes as a flat data.frame.
# Rows = cells, columns = cell_id + one column per gene present in the object.
extract_expr_df <- function(srt, genes, assay = "SCT") {
  genes_ok <- intersect(genes, rownames(srt))
  mat <- Seurat::GetAssayData(srt, assay = assay, layer = "data")[genes_ok, , drop = FALSE]
  df  <- as.data.frame(t(as.matrix(mat)))
  df  <- cbind(data.frame(cell_id = rownames(df)), df)
  rownames(df) <- NULL
  df
}

# Extract UMAP coords + cell type label from a reference Seurat object.
# Used to build the grey silhouette and full reference UMAP panel without
# loading the full Seurat object in notebooks.
extract_ref_umap_df <- function(srt_ref, umap_col = "refUMAP", ct_col = "cell_type") {
  emb <- Seurat::Embeddings(srt_ref, umap_col) %>% as.data.frame()
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb <- emb %>%
    rownames_to_column("cell_id") %>%
    left_join(srt_ref@meta.data %>% rownames_to_column("cell_id"))
  emb
}

extract_azimuth_umap_df <- function(ref_path, umap_col = "refUMAP", ct_col = "cell_type") {
  srt_ref <- readRDS(glue::glue("{ref_path}/ref.Rds"))
  emb <- Seurat::Embeddings(srt_ref, umap_col) %>%
    as.data.frame()
  # Normalise UMAP column names: refUMAP_1/refumap_1 → UMAP1/UMAP2
  names(emb) <- sub("(?i)^ref[uU][mM][aA][pP]_?1$", "UMAP1",
                    sub("(?i)^ref[uU][mM][aA][pP]_?2$", "UMAP2",
                        names(emb), perl = TRUE), perl = TRUE)
  meta <- srt_ref@meta.data %>% rownames_to_column("cell_id")
  # Standardise the cell type column to "cell_type".
  # Try: caller-supplied ct_col, then common aliases, then first non-id column.
  ct_candidates <- c(ct_col, "celltype", "CellType", "cell.type",
                     "annotation", "label", "cluster_label")
  found <- intersect(ct_candidates, colnames(meta))[1]
  if (!is.na(found) && found != "cell_type") {
    meta <- dplyr::rename(meta, cell_type = !!found)
  }
  emb %>%
    rownames_to_column("cell_id") %>%
    left_join(meta)
}

# ── Data-frame based UMAP helpers ────────────────────────────────────────────
# These functions take plain data.frames (from anno_meta_* targets) instead of
# Seurat objects — avoiding large object loads in notebooks.
#
# df: data.frame with at minimum UMAP1, UMAP2 columns and any metadata columns.

#' Colour UMAP by a discrete or continuous metadata column.
#' @param df       anno_meta data.frame (UMAP1, UMAP2, metadata cols).
#' @param col      Column name to colour by.
#' @param cols     Named colour vector (optional; auto-generated if NULL).
#' @param pt_size  Point size (default 0.4).
#' @param title    Plot title.
#' @param label    If TRUE, overlay cluster centroids with ggrepel text.
plot_umap_df <- function(df, col, cols = NULL, pt_size = 0.4,
                         title = NULL, label = FALSE, facet_by = NULL) {
  is_cont <- is.numeric(df[[col]])
  keep_cols <- unique(c("UMAP1", "UMAP2", col, facet_by))
  keep_cols <- intersect(keep_cols, colnames(df))
  d <- df[, keep_cols, drop = FALSE]
  if (!is_cont) d[[col]] <- factor(d[[col]])

  p <- ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2, colour = .data[[col]])) +
    ggplot2::geom_point(size = pt_size, alpha = 0.7) +
    ggplot2::labs(title = title, x = NULL, y = NULL, colour = NULL) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))

  if (is_cont) {
    p <- p + ggplot2::scale_colour_viridis_c(option = "viridis")
  } else if (!is.null(cols)) {
    p <- p + ggplot2::scale_colour_manual(values = cols, na.value = "grey70")
  }

  if (label && !is_cont) {
    d_nona <- d[!is.na(d[[col]]), ]
    cents  <- aggregate(cbind(UMAP1, UMAP2) ~ d_nona[[col]], data = d_nona, FUN = median)
    colnames(cents)[1] <- col
    p <- p + ggrepel::geom_text_repel(
      data = cents, ggplot2::aes(label = .data[[col]]),
      size = 2.8, fontface = "bold", max.overlaps = 20,
      segment.size = 0.2, colour = "grey20")
  }

  if (!is.null(facet_by) && facet_by %in% colnames(d))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)), nrow = 1)

  p
}

#' Feature UMAP from anno_expr + anno_meta — split by a metadata column.
#' @param meta_df  anno_meta data.frame (UMAP1, UMAP2, split_col, ...).
#' @param expr_df  anno_expr data.frame (cell_id + gene columns).
#' @param features Character vector of gene names.
#' @param split_by Metadata column to facet by (NULL = no facet).
#' @param pt_size  Point size (default 0.4).
plot_umap_feature_df <- function(meta_df, expr_df, features, split_by = NULL,
                                 pt_size = 0.4) {
  coords <- merge(meta_df[, c("cell_id", "UMAP1", "UMAP2",
                               intersect(split_by, colnames(meta_df)))],
                  expr_df, by = "cell_id")

  plots <- lapply(features, function(gene) {
    if (!gene %in% colnames(coords)) {
      return(ggplot2::ggplot() +
        ggplot2::annotate("text", x=0.5, y=0.5,
                          label=paste0(gene, "\n(absent)"), size=3) +
        ggplot2::theme_void() + ggplot2::labs(title = gene))
    }
    d <- coords[order(coords[[gene]]), ]
    p <- ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2, colour = .data[[gene]])) +
      ggplot2::geom_point(size = pt_size, alpha = 0.8) +
      ggplot2::scale_colour_viridis_c(option = "mako", name = gene) +
      ggplot2::labs(title = gene, x = NULL, y = NULL) +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                     legend.key.size = ggplot2::unit(0.4, "cm"),
                     legend.text = ggplot2::element_text(size = 7))
    if (!is.null(split_by) && split_by %in% colnames(d))
      p <- p + ggplot2::facet_wrap(~ .data[[split_by]], nrow = 1)
    p
  })
  if (length(plots) == 1) plots[[1]] else patchwork::wrap_plots(plots, ncol = 1)
}

#' Grid UMAP: sample.name × dm2_class — plain and highlighted.
#' @param df  anno_meta data.frame (UMAP1, UMAP2, sample.name, dm2_class).
plot_umap_grid_sample_dm2_df <- function(df) {
  d <- df[, c("UMAP1", "UMAP2", "sample.name", "dm2_class")]
  d$sample.name <- factor(d$sample.name)
  d$dm2_class   <- factor(d$dm2_class)

  base_theme <- ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio  = 1, strip.text = ggplot2::element_text(size = 9),
                   plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
                   panel.spacing = ggplot2::unit(0.2, "lines"))

  p_plain <- ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(colour = "grey40", size = 0.3, alpha = 0.4) +
    ggplot2::facet_grid(sample.name ~ dm2_class) +
    ggplot2::labs(title = "sample \u00d7 dm2_class (plain)") +
    base_theme

  d_bg <- d[, c("UMAP1", "UMAP2")]
  p_highlight <- ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(data = d_bg, colour = "grey85", size = 0.3, alpha = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = dm2_class), size = 0.3) +
    ggplot2::facet_grid(sample.name ~ dm2_class) +
    ggplot2::labs(title = "sample \u00d7 dm2_class (highlight)", colour = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 4, alpha = 1))) +
    base_theme

  list(plain = p_plain, highlight = p_highlight)
}

# ── SCENIC pipeline helpers ───────────────────────────────────────────────────

#' Load pySCENIC AUC matrix from CSV.
#' @param path Path to auc_mtx.csv (cells × regulons, cell IDs in first column).
#' @return Numeric matrix, rownames = cell IDs, colnames = regulon names.
scenic_load_auc <- function(path) {
  dt  <- data.table::fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  mat
}

#' Build slim SCENIC metadata from anno_meta.
#' Keeps cell_id, UMAP coords, sample.name, the requested cluster column, and
#' adds a short \code{condition} label parsed from \code{sample.name}.
#' @param anno_meta  data.frame from the \code{anno_meta_*} pipeline target.
#' @param clust_col  Name of the cluster column to retain.
#' @return data.frame.
scenic_build_meta <- function(anno_meta, clust_col) {
  anno_meta[, c("cell_id", "UMAP1", "UMAP2", "sample.name", clust_col)] |>
    dplyr::mutate(condition = dplyr::case_when(
      grepl("mCh\\+", sample.name) ~ "mCh+",
      grepl("mCh-",   sample.name) ~ "mCh-",
      grepl("WT",     sample.name) ~ "WT",
      TRUE ~ sample.name
    ))
}

#' Aggregate AUC matrix to cluster means.
#' @param auc_mat   cells × regulons numeric matrix.
#' @param meta      data.frame with \code{cell_id} and \code{clust_col}.
#' @param clust_col Metadata column used for grouping.
#' @return clusters × regulons numeric matrix.
scenic_aggregate_auc <- function(auc_mat, meta, clust_col) {
  shared <- intersect(rownames(auc_mat), meta$cell_id)
  auc    <- auc_mat[shared, , drop = FALSE]
  groups <- as.character(meta[[clust_col]][match(shared, meta$cell_id)])
  t(vapply(
    split(seq_len(nrow(auc)), groups),
    function(idx) colMeans(auc[idx, , drop = FALSE]),
    numeric(ncol(auc))
  ))
}

#' Condition proportion matrix for SCENIC heatmap annotation.
#' @param meta        data.frame with \code{clust_col} and \code{condition}.
#' @param clust_col   Grouping column.
#' @param clust_order Character vector giving the desired row order (typically
#'                    \code{rownames(scenic_agg_hsc)}).
#' @param cond_order  Condition levels in display order.
#' @return clusters × conditions numeric matrix (rows sum to 1).
scenic_cond_prop_mat <- function(meta, clust_col, clust_order,
                                 cond_order = c("WT", "mCh-", "mCh+")) {
  meta |>
    dplyr::count(.data[[clust_col]], condition) |>
    dplyr::group_by(.data[[clust_col]]) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup() |>
    tidyr::complete(
      !!rlang::sym(clust_col) := clust_order,
      condition = cond_order,
      fill = list(prop = 0)
    ) |>
    dplyr::select(cluster = .data[[clust_col]], condition, prop) |>
    tidyr::pivot_wider(names_from = condition, values_from = prop,
                       values_fill = 0) |>
    dplyr::select(cluster, dplyr::any_of(cond_order)) |>
    dplyr::arrange(match(cluster, clust_order)) |>
    tibble::column_to_rownames("cluster") |>
    as.matrix()
}

#' Regulon Specificity Score (RSS) — reimplementation of SCENIC::calcRSS().
#'
#' Measures how specifically each regulon's AUC activity is concentrated in
#' each cell group (cluster, condition, etc.) using Jensen-Shannon Divergence.
#'
#' Algorithm (Marín-Béjar & Bhatt et al. 2021):
#'   1. Row-normalise the AUC matrix so each regulon sums to 1 across cells →
#'      treats the regulon's activity as a discrete probability distribution P.
#'   2. For each cell group, build a uniform distribution Q over its cells
#'      (1 / n_cells_in_group for member cells, 0 otherwise).
#'   3. RSS = 1 − sqrt(JSD(P ‖ Q)), where JSD is the Jensen-Shannon divergence.
#'      RSS = 1 means the regulon is perfectly co-localised with the group;
#'      RSS → 0 means it is no more specific to that group than the background.
#'
#' Input convention: cells × regulons (matching scenic_auc_mat_hsc in this
#' project). The original SCENIC package uses regulons × cells — this function
#' transposes internally.
#'
#' @param auc_mat        Numeric matrix, cells × regulons. Accepts the
#'                       \code{scenic_auc_mat_hsc} target directly.
#' @param cell_annotation Named character/factor vector (names = cell IDs) OR
#'                       an unnamed vector aligned to \code{rownames(auc_mat)}.
#'                       Defines the grouping variable (cluster, condition, …).
#' @param cell_types     Optional character vector of group levels to compute.
#'                       Defaults to all unique values in \code{cell_annotation}.
#' @return Numeric matrix, regulons × cell_types, values in [0, 1].
#'   Higher = more specific to that cell type.
calc_rss <- function(auc_mat, cell_annotation, cell_types = NULL) {
  # ── Input checks ─────────────────────────────────────────────────────────────
  if (any(is.na(cell_annotation)))
    stop("NAs in cell_annotation — remove or impute before calling calc_rss().")

  # Align annotation to matrix rows if named
  if (!is.null(names(cell_annotation))) {
    shared <- intersect(names(cell_annotation), rownames(auc_mat))
    if (length(shared) == 0)
      stop("No overlap between names(cell_annotation) and rownames(auc_mat).")
    auc_mat         <- auc_mat[shared, , drop = FALSE]
    cell_annotation <- cell_annotation[shared]
  } else {
    if (length(cell_annotation) != nrow(auc_mat))
      stop("Length of cell_annotation must equal nrow(auc_mat) when unnamed.")
  }

  cell_annotation <- as.character(cell_annotation)
  if (is.null(cell_types)) cell_types <- sort(unique(cell_annotation))

  # ── Normalise AUC: each regulon sums to 1 across cells (row is a regulon) ──
  # Work in regulons × cells orientation to match original SCENIC convention.
  norm_auc <- t(auc_mat)                       # regulons × cells
  reg_sums  <- rowSums(norm_auc)
  zero_regs <- reg_sums == 0
  if (any(zero_regs))
    warning(sum(zero_regs), " regulon(s) have zero total AUC and will yield RSS = 0.")
  norm_auc[!zero_regs, ] <- norm_auc[!zero_regs, ] / reg_sums[!zero_regs]

  # ── Shannon entropy ───────────────────────────────────────────────────────────
  .H <- function(p) {
    p <- p[p > 0]
    -sum(p * log2(p))
  }

  # ── Jensen-Shannon Divergence ─────────────────────────────────────────────────
  .jsd <- function(p, q) {
    m <- (p + q) / 2
    .H(m) - (.H(p) + .H(q)) / 2
  }

  # ── RSS per regulon × cell type ───────────────────────────────────────────────
  regulons <- rownames(norm_auc)

  rss_cols <- lapply(cell_types, function(ct) {
    # Uniform distribution over cells belonging to this type
    is_ct <- as.numeric(cell_annotation == ct)
    n_ct  <- sum(is_ct)
    if (n_ct == 0) return(setNames(rep(NA_real_, length(regulons)), regulons))
    q <- is_ct / n_ct

    vapply(regulons, function(reg) {
      p   <- norm_auc[reg, ]
      jsd <- .jsd(p, q)
      1 - sqrt(max(jsd, 0))   # clamp numerical noise below 0
    }, numeric(1L))
  })

  rss <- do.call(cbind, rss_cols)
  colnames(rss) <- cell_types
  rss
}

#' Plot top-N regulons per cell group as a dot / bubble plot.
#'
#' Takes the matrix returned by \code{calc_rss()} and plots the top \code{top_n}
#' regulons per cell type, sized and coloured by RSS score.
#'
#' @param rss_mat   Regulons × cell_types matrix from \code{calc_rss()}.
#' @param top_n     Number of top regulons to show per cell type. Default 5.
#' @param min_rss   Drop regulon × cell combinations with RSS < this value.
#'                  Default 0.
#' @param strip_dir Logical. Strip the "(+)" / "(-)" suffix from regulon labels.
#'                  Default TRUE.
#' @return A ggplot object.
plot_rss <- function(rss_mat, top_n = 5, min_rss = 0, strip_dir = TRUE) {
  df <- as.data.frame(rss_mat) |>
    tibble::rownames_to_column("regulon") |>
    tidyr::pivot_longer(-regulon, names_to = "cell_type", values_to = "rss") |>
    dplyr::filter(!is.na(rss), rss >= min_rss)

  if (strip_dir)
    df$regulon <- sub("\\([+-]\\)$", "", df$regulon)

  top_df <- df |>
    dplyr::group_by(cell_type) |>
    dplyr::slice_max(rss, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()

  # Regulon order: cluster by cell type of maximum RSS
  reg_order <- top_df |>
    dplyr::group_by(regulon) |>
    dplyr::slice_max(rss, n = 1, with_ties = FALSE) |>
    dplyr::arrange(cell_type, desc(rss)) |>
    dplyr::pull(regulon)

  top_df$regulon   <- factor(top_df$regulon,   levels = rev(unique(reg_order)))
  top_df$cell_type <- factor(top_df$cell_type, levels = colnames(rss_mat))

  ggplot2::ggplot(top_df, ggplot2::aes(cell_type, regulon, size = rss, fill = rss)) +
    ggplot2::geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
    ggplot2::scale_size_continuous(range = c(1, 7), name = "RSS") +
    ggplot2::scale_fill_gradientn(
      colours = c("#f7f7f7", "#fee090", "#d73027"),
      limits  = c(min_rss, 1), name = "RSS"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = ggplot2::element_text(size = 8),
      axis.title   = ggplot2::element_blank(),
      legend.position = "right"
    )
}

# ── Cluster annotation helpers ────────────────────────────────────────────────
# Used by notebooks/cluster-annotation.qmd.
# All functions take plain data.frames — no Seurat objects.

# Convert a resolution number to the column name used in anno_meta / recluster dfs.
anno_res_col <- function(r) paste0("res_", gsub("\\.", "_", r))

# Merge anno_meta with annotation_recluster.
# annotation_recluster adds extra res_* columns from a wider resolution sweep.
anno_merge_recluster <- function(meta_df, recluster_df) {
  new_cols <- setdiff(colnames(recluster_df), colnames(meta_df))
  if (length(new_cols) == 0) return(meta_df)
  extra <- recluster_df[, c("cell_id", new_cols), drop = FALSE]
  merge(meta_df, extra, by = "cell_id", all.x = TRUE)
}

# Merge a meta_df with one projection result (proj_<ref_id>_<subset> target).
anno_merge_proj <- function(meta_df, proj_df) {
  inner_join(meta_df, proj_df)
}

# Load all proj_<ref_id>_<subset> flat tables into a named list.
# Missing targets return NULL (reference not yet computed).
anno_load_proj_list <- function(subset, ref_ids) {
  proj_list <- vector("list", length(ref_ids))
  names(proj_list) <- ref_ids
  for (ref_id in ref_ids) {
    tgt <- paste0("proj_", ref_id, "_", subset)
    proj_list[[ref_id]] <- tryCatch(
      targets::tar_read_raw(tgt),
      error = function(e) { message("Skipping ", tgt, ": ", e$message); NULL }
    )
  }
  proj_list
}

# Load all ref_umap_df_<ref_id> flat tables into a named list.
anno_load_ref_df_list <- function(ref_ids) {
  ref_list <- vector("list", length(ref_ids))
  names(ref_list) <- ref_ids
  for (ref_id in ref_ids) {
    tgt <- paste0("ref_umap_df_", ref_id)
    ref_list[[ref_id]] <- tryCatch(
      targets::tar_read_raw(tgt),
      error = function(e) NULL
    )
  }
  ref_list
}

# UMAP coloured by cluster, one panel per resolution, arranged in a 4-column grid.
plot_resolution_umap_grid <- function(meta_df, resolutions, compartment) {
  panel_list <- vector("list", length(resolutions))
  for (i in seq_along(resolutions)) {
    r   <- resolutions[i]
    col <- anno_res_col(r)
    if (!col %in% colnames(meta_df)) {
      panel_list[[i]] <- ggplot2::ggplot() + ggplot2::theme_void()
      next
    }
    d <- meta_df[, c("UMAP1", "UMAP2", col)]
    d$cluster <- factor(d[[col]])
    cents <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = d, FUN = median)
    panel_list[[i]] <-
      ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2, colour = cluster)) +
      ggplot2::geom_point(size = 0.25, alpha = 0.7) +
      ggrepel::geom_text_repel(data = cents, ggplot2::aes(label = cluster),
        size = 3.2, fontface = "bold", max.overlaps = 20,
        segment.size = 0.2, colour = "grey20") +
      ggplot2::labs(title = paste0("res ", r)) +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(size = 9, face = "bold"))
  }
  patchwork::wrap_plots(panel_list, ncol = 4) +
    patchwork::plot_annotation(
      title = paste0(compartment, " — clustering across resolutions"))
}

# Feature UMAP panels for a set of genes.
# Returns a named list of ggplot objects; absent genes get a placeholder panel.
plot_feature_umap_panels <- function(meta_df, expr_df, genes, compartment) {
  coords   <- merge(meta_df[, c("cell_id", "UMAP1", "UMAP2")], expr_df, by = "cell_id")
  genes_ok <- intersect(genes, colnames(coords))
  genes_na <- setdiff(genes, genes_ok)
  panel_list <- vector("list", length(genes))
  names(panel_list) <- genes
  for (g in genes_ok) {
    d <- coords[order(coords[[g]]), ]
    panel_list[[g]] <-
      ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2, colour = .data[[g]])) +
      ggplot2::geom_point(size = 0.25, alpha = 0.8) +
      ggplot2::scale_colour_viridis_c(option = "mako", name = g) +
      ggplot2::labs(title = paste0(compartment, " — ", g)) +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9),
                     legend.key.size = ggplot2::unit(0.35, "cm"),
                     legend.text = ggplot2::element_text(size = 7))
  }
  for (g in genes_na) {
    panel_list[[g]] <-
      ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = paste0(g, "\n(absent)"), size = 3) +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0(compartment, " — ", g))
  }
  panel_list
}

# Build a colour palette shared between reference and query UMAPs.
# label: cell-type column in ref_df. ref_df may be NULL (cross-species, palette from pred only).
anno_make_type_pal <- function(ref_df, pred_types, label = "cell_type") {
  ref_types   <- if (!is.null(ref_df)) unique(ref_df[[label]]) else character(0)
  query_types <- unique(na.omit(pred_types))
  all_types   <- sort(union(ref_types, query_types))
  if (length(all_types) == 0L) return(setNames(character(0), character(0)))
  setNames(scales::hue_pal()(length(all_types)), all_types)
}

# Reference UMAP coloured by cell type with centroid labels.
# label:   column name in ref_df holding cell type identities.
# ref_key: prefix for ref UMAP axis columns; ref_key+"1" and ref_key+"2" are used.
plot_ref_umap <- function(ref_df, type_pal, ref_name, label = "cell_type",
                          ref_key = "UMAP") {
  df       <- ref_df
  df$.type <- ref_df[[label]]
  df$.x    <- ref_df[[paste0(ref_key, 1)]]
  df$.y    <- ref_df[[paste0(ref_key, 2)]]
  cents    <- aggregate(cbind(.x, .y) ~ .type, data = df, FUN = median)
  ggplot2::ggplot(df, ggplot2::aes(.x, .y, colour = .type)) +
    ggplot2::geom_point(size = 0.2, alpha = 0.5) +
    ggrepel::geom_text_repel(data = cents, ggplot2::aes(label = .type),
      size = 2.8, fontface = "bold", max.overlaps = 20,
      segment.size = 0.2, colour = "grey20") +
    ggplot2::scale_colour_manual(values = type_pal, name = NULL) +
    ggplot2::labs(title = paste0("Reference — ", ref_name)) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))
}

# Query cells projected onto reference UMAP space.
# label:    ref cell-type column; predicted column is paste0("predicted.", label).
# proj_key: prefix for projected UMAP columns in df  (proj_key+"1", proj_key+"2").
# ref_key:  prefix for reference UMAP columns in ref_df (ref_key+"1", ref_key+"2").
plot_ref_proj_umap <- function(df, ref_df, type_pal, label = "cell_type",
                               proj_key = "refumap_", ref_key = "UMAP") {
  pred_col  <- paste0("predicted.", label)
  proj      <- data.frame(.x = df[[paste0(proj_key, 1)]], .y = df[[paste0(proj_key, 2)]],
                           .type = df[[pred_col]])
  cents     <- aggregate(cbind(.x, .y) ~ .type, data = proj[!is.na(proj$.type), ], FUN = median)

  p <- ggplot2::ggplot()
  if (!is.null(ref_df)) {
    ref_bg <- data.frame(.x = ref_df[[paste0(ref_key, 1)]], .y = ref_df[[paste0(ref_key, 2)]])
    p <- p + ggplot2::geom_point(data = ref_bg, ggplot2::aes(.x, .y),
                                 colour = "grey82", size = 0.12, alpha = 0.35)
  }
  p +
    ggplot2::geom_point(data = proj, ggplot2::aes(.x, .y, colour = .type),
                        size = 0.3, alpha = 0.7) +
    ggrepel::geom_text_repel(data = cents, ggplot2::aes(.x, .y, label = .type),
      size = 2.8, fontface = "bold", max.overlaps = 20,
      segment.size = 0.2, colour = "grey20") +
    ggplot2::scale_colour_manual(values = type_pal, na.value = "grey70") +
    ggplot2::labs(title = "Query projected on ref UMAP") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))
}

# Query cell density overlaid on reference UMAP space.
# proj_key: prefix for projected UMAP columns in df.
# ref_key:  prefix for reference UMAP columns in ref_df.
plot_query_density_umap <- function(df, ref_df,
                                    proj_key = "refumap_", ref_key = "UMAP") {
  dens_df <- data.frame(.x = df[[paste0(proj_key, 1)]], .y = df[[paste0(proj_key, 2)]])
  ref_bg  <- data.frame(.x = ref_df[[paste0(ref_key, 1)]], .y = ref_df[[paste0(ref_key, 2)]])
  ggplot2::ggplot() +
    ggplot2::geom_point(data = ref_bg, ggplot2::aes(.x, .y),
                        colour = "grey82", size = 0.12, alpha = 0.35) +
    ggpointdensity::geom_pointdensity(data = dens_df, ggplot2::aes(.x, .y), size = 0.12) +
    ggplot2::scale_colour_viridis_c(option = "inferno", name = "Density") +
    ggplot2::labs(title = glue::glue("{unique(df$sample.name)} density on ref UMAP")) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.35, "cm"),
                   legend.text = ggplot2::element_text(size = 7))
}

# Query UMAP coloured by predicted cell type labels.
# label: ref cell-type column; predicted column is paste0("predicted.", label).
plot_query_label_umap <- function(df, type_pal, label = "cell_type") {
  pred_col <- paste0("predicted.", label)
  q     <- data.frame(UMAP1 = df$UMAP1, UMAP2 = df$UMAP2, .type = df[[pred_col]])
  cents <- aggregate(cbind(UMAP1, UMAP2) ~ .type, data = q[!is.na(q$.type), ], FUN = median)
  ggplot2::ggplot(q, ggplot2::aes(UMAP1, UMAP2, colour = .type)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.6) +
    ggrepel::geom_text_repel(data = cents, ggplot2::aes(label = .type),
      size = 2.8, fontface = "bold", max.overlaps = 20,
      segment.size = 0.2, colour = "grey20") +
    ggplot2::scale_colour_manual(values = type_pal, na.value = "grey70") +
    ggplot2::labs(title = "Query — predicted labels") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))
}

# Mapping score UMAP (turbo palette).
plot_mapping_score_umap <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, colour = mapping.score)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.7) +
    ggplot2::scale_colour_viridis_c(option = "turbo", name = "Mapping\nscore") +
    ggplot2::labs(title = "Mapping score") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.4, "cm"),
                   legend.text = ggplot2::element_text(size = 7))
}

# Prediction score UMAP (plasma palette).
# label: used to derive the score column as paste0("predicted.", label, ".score").
plot_prediction_score_umap <- function(df, label = "cell_type") {
  score_col <- paste0("predicted.", label, ".score")
  ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, colour = .data[[score_col]])) +
    ggplot2::geom_point(size = 0.25, alpha = 0.7) +
    ggplot2::scale_colour_viridis_c(option = "plasma", name = "Prediction\nscore") +
    ggplot2::labs(title = "Prediction score") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                   legend.position = "right",
                   legend.key.height = ggplot2::unit(0.4, "cm"),
                   legend.text = ggplot2::element_text(size = 7))
}

# Print the full 3-row projection panel for one reference.
# Row 1: ref UMAP + query projected on ref + shared legend
# Row 2: query density on ref + query label UMAP
# Row 3: mapping score + prediction score
# label:    ref cell-type column; predicted column is paste0("predicted.", label).
# proj_key: prefix for projected UMAP columns in df (proj_key+"1", proj_key+"2").
# ref_key:  prefix for reference UMAP columns in ref_df (ref_key+"1", ref_key+"2").
print_projection_panel <- function(df, ref_df, ref_name, label = "cell_type",
                                   proj_key = "refumap_", ref_key = "UMAP") {
  pred_col <- paste0("predicted.", label)
  type_pal <- anno_make_type_pal(ref_df, df[[pred_col]], label)
  p_ref    <- plot_ref_umap(ref_df, type_pal, ref_name, label, ref_key)
  p_proj   <- plot_ref_proj_umap(df, ref_df, type_pal, label, proj_key, ref_key)

  p_legend <- p_ref +
    ggplot2::theme(legend.position = "right",
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(0.3, "cm")) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1), ncol = 1))
  p_legend <- cowplot::get_legend(p_legend)

  p_density <- plot_query_density_umap(df, ref_df, proj_key, ref_key)
  p_query   <- plot_query_label_umap(df, type_pal, label)
  p_map     <- plot_mapping_score_umap(df)
  p_pscore  <- plot_prediction_score_umap(df, label)

  row1 <- cowplot::plot_grid(
    cowplot::plot_grid(p_ref, p_proj, nrow = 1),
    p_legend, nrow = 1, rel_widths = c(1, 0.22))
  row2 <- cowplot::plot_grid(p_density, p_query, nrow = 1)
  row3 <- cowplot::plot_grid(p_map, p_pscore, nrow = 1)

  print(cowplot::plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1, 1)))
  invisible(TRUE)
}

# Compute adjusted Rand index between predicted cell type and each cluster resolution.
# label: ref cell-type column; predicted column is paste0("predicted.", label).
compute_ari_table <- function(df, resolutions, label = "cell_type") {
  pred_col <- paste0("predicted.", label)
  pred <- df[[pred_col]]
  rows <- vector("list", length(resolutions))
  for (i in seq_along(resolutions)) {
    r     <- resolutions[i]
    clust <- df[[anno_res_col(r)]]
    keep  <- !is.na(pred) & !is.na(clust)
    ari   <- if (sum(keep) > 1) mclust::adjustedRandIndex(pred[keep], clust[keep]) else NA_real_
    rows[[i]] <- data.frame(resolution = r, ari = ari)
  }
  do.call(rbind, rows)
}

plot_ari_line <- function(ari_df, title = NULL) {
  best_r <- ari_df$resolution[which.max(ari_df$ari)]
  ggplot2::ggplot(ari_df, ggplot2::aes(x = resolution, y = ari)) +
    ggplot2::geom_line(colour = "grey50", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = resolution == best_r), size = 2.5) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey40", "TRUE" = "#d73027"), guide = "none") +
    ggplot2::geom_vline(xintercept = best_r, colour = "#d73027",
                        linewidth = 0.4, linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = ari_df$resolution) +
    ggplot2::labs(x = "Clustering resolution", y = "Adjusted Rand index",
                  title = title,
                  caption = paste0("Best resolution: ", best_r)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# Mapping score bar + violin by a grouping variable (cluster or cell type).
# score_df must have columns: <x_col> (grouping), score, condition.
plot_mapping_score_violin <- function(score_df, x_col, x_label, title) {
  n_conditions <- length(unique(score_df$condition))
  dodge_w      <- 0.85
  box_w        <- dodge_w / n_conditions * 0.45

  p_bar <- ggplot2::ggplot(score_df,
      ggplot2::aes(.data[[x_col]], fill = condition)) +
    ggplot2::geom_bar(position = ggplot2::position_dodge(dodge_w),
                      width = dodge_w * 0.9, alpha = 0.8,
                      colour = "grey20", linewidth = 0.3) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Condition") +
    ggplot2::scale_y_continuous(transform = "log10") +
    ggplot2::labs(y = "# cells", x = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position  = "top",
                   legend.direction = "horizontal",
                   legend.text      = ggplot2::element_text(size = 7),
                   legend.key.size  = ggplot2::unit(0.35, "cm"),
                   panel.grid.major.x = ggplot2::element_blank())

  p_vln <- ggplot2::ggplot(score_df,
      ggplot2::aes(.data[[x_col]], score, fill = condition)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(dodge_w),
                         scale = "width", trim = TRUE,
                         linewidth = 0.25, alpha = 0.8) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(dodge_w),
                          width = box_w, outlier.size = 0.2,
                          colour = "grey20", linewidth = 0.3) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Condition") +
    ggplot2::labs(y = "Mapping score", x = x_label) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11),
                   legend.position = "none",
                   panel.grid.major.x = ggplot2::element_blank())

  p_bar / p_vln + patchwork::plot_layout(heights = c(1, 2.5)) +
    patchwork::plot_annotation(title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10, face = "bold")))
}

# Mapping score violin by cluster.
plot_score_by_cluster <- function(df, ref_name, compartment, resolution = 0.5) {
  clust_col <- anno_res_col(resolution)
  score_df <- df %>%
    dplyr::filter(!is.na(mapping.score), !is.na(.data[[clust_col]])) %>%
    dplyr::transmute(cluster   = factor(.data[[clust_col]]),
                     score     = mapping.score,
                     condition = sample.name)
  title <- paste0(compartment, " — ", ref_name, " (res ", resolution, ")")
  plot_mapping_score_violin(score_df, "cluster", clust_col, title)
}

# Mapping score violin by predicted cell type; ordered by median, rare types removed.
# label: ref cell-type column; predicted column is paste0("predicted.", label).
plot_score_by_celltype <- function(df, ref_name, compartment, label = "cell_type",
                                   min_cells = 10) {
  pred_col <- paste0("predicted.", label)
  score_df <- df %>%
    dplyr::filter(!is.na(mapping.score), !is.na(.data[[pred_col]])) %>%
    dplyr::transmute(cell_type = .data[[pred_col]],
                     score     = mapping.score,
                     condition = sample.name)
  keep     <- score_df$cell_type %>% table() %>% .[. >= min_cells] %>% names()
  score_df <- score_df %>% dplyr::filter(cell_type %in% keep)
  med_order          <- tapply(score_df$score, score_df$cell_type, median) %>%
                        sort(decreasing = TRUE) %>% names()
  score_df$cell_type <- factor(score_df$cell_type, levels = med_order)
  title <- paste0(compartment, " — ", ref_name, " mapping score by cell type")
  plot_mapping_score_violin(score_df, "cell_type", NULL, title)
}

# Resolution hierarchy heatmap: predicted cell type (rows) × finest cluster (cols).
# Column annotation bars show the majority coarser-resolution cluster for each fine cluster.
plot_resolution_heatmap <- function(df, pred_col, resolutions, title = NULL) {
  res_col_map <- sapply(resolutions, anno_res_col)
  names(res_col_map) <- as.character(resolutions)

  fine_col <- res_col_map[as.character(max(resolutions))]
  if (!fine_col %in% colnames(df) || !pred_col %in% colnames(df))
    return(invisible(NULL))

  pred <- df[[pred_col]]
  fine <- df[[fine_col]]
  keep <- !is.na(pred) & !is.na(fine)
  if (sum(keep) == 0) return(invisible(NULL))

  count_tbl    <- table(predicted = pred[keep], cluster = fine[keep])
  col_prop_mat <- t(t(count_tbl) / colSums(count_tbl))
  cell_labels  <- matrix(
    ifelse(count_tbl > 0, as.character(count_tbl), ""),
    nrow = nrow(count_tbl), dimnames = dimnames(count_tbl))

  # Majority coarser-resolution cluster for each fine cluster
  coarser_res  <- sort(resolutions[resolutions < max(resolutions)])
  col_annot_df <- data.frame(row.names = colnames(count_tbl))
  for (r in coarser_res) {
    coarse_col <- res_col_map[as.character(r)]
    coarse     <- df[[coarse_col]]
    majority   <- character(ncol(count_tbl))
    for (j in seq_len(ncol(count_tbl))) {
      fc          <- colnames(count_tbl)[j]
      idx         <- !is.na(fine) & fine == fc & !is.na(coarse)
      majority[j] <- if (any(idx)) names(which.max(table(coarse[idx]))) else NA_character_
    }
    col_annot_df[[paste0("res.", r)]] <- majority
  }

  # Per-column annotation colour palettes
  annot_col_list <- list()
  for (cn in colnames(col_annot_df)) {
    lvls <- sort(unique(na.omit(col_annot_df[[cn]])), method = "radix")
    annot_col_list[[cn]] <- setNames(scales::hue_pal()(length(lvls)), lvls)
  }

  # ComplexHeatmap column annotations (one coloured bar per coarser resolution)
  anno_fun_list <- list()
  for (cn in colnames(col_annot_df)) {
    vals   <- col_annot_df[[cn]]
    pal    <- annot_col_list[[cn]]
    n_cols <- nrow(col_annot_df)
    anno_fun_list[[cn]] <- ComplexHeatmap::AnnotationFunction(
      fun = function(index, k, n) {
        n_col <- length(index)
        for (i in seq_along(index)) {
          v   <- vals[index[i]]
          col <- if (!is.na(v) && v %in% names(pal)) pal[[v]] else "grey90"
          grid::grid.rect(
            x      = grid::unit((i - 0.5) / n_col, "npc"),
            y      = grid::unit(0.5, "npc"),
            width  = grid::unit(1 / n_col, "npc"),
            height = grid::unit(1, "npc"),
            gp     = grid::gpar(fill = col, col = "black", lwd = 0.5))
          if (!is.na(v))
            grid::grid.text(v,
              x  = grid::unit((i - 0.5) / n_col, "npc"),
              y  = grid::unit(0.5, "npc"),
              gp = grid::gpar(fontsize = 6, fontface = "bold"))
        }
      },
      which = "column", n = n_cols, height = grid::unit(5, "mm"))
  }

  col_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
    c(anno_fun_list,
      list(annotation_name_gp = grid::gpar(fontsize = 7), show_legend = FALSE)))

  col_fun <- circlize::colorRamp2(c(0, 1), c("white", "#2166ac"))

  ComplexHeatmap::Heatmap(
    col_prop_mat,
    name              = "Cluster prop.",
    col               = col_fun,
    top_annotation    = col_anno,
    rect_gp           = grid::gpar(col = "grey60", lwd = 0.5),
    cell_fun          = function(j, i, x, y, width, height, fill) {
      lbl <- cell_labels[i, j]
      if (nchar(lbl) > 0)
        grid::grid.text(lbl, x, y, gp = grid::gpar(fontsize = 6, fontface = "bold"))
    },
    cluster_rows            = TRUE,
    cluster_columns         = TRUE,
    clustering_method_rows  = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_row_dend           = TRUE,
    show_column_dend        = TRUE,
    row_names_gp            = grid::gpar(fontsize = 8),
    column_names_gp         = grid::gpar(fontsize = 8),
    column_names_rot        = 45,
    column_title            = title,
    column_title_gp         = grid::gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param    = list(
      title    = "Cluster prop.",
      title_gp = grid::gpar(fontsize = 8)))
}

# ── Per-reference section helpers (used by cluster-annotation.qmd) ───────────
#
# Each function draws one subsection for a single reference × all compartments.
# Called once per reference chunk in the notebook — no loops in the notebook.

# Shared legend strip extracted from a reference UMAP.
ref_make_legend <- function(p_ref) {
  cowplot::get_legend(
    p_ref +
      ggplot2::theme(legend.position  = "right",
                     legend.text      = ggplot2::element_text(size = 7),
                     legend.key.size  = ggplot2::unit(0.3, "cm")) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1), ncol = 1)))
}

# Draw reference UMAP (left) + query cells projected on it (right) + legend.
ref_plot_labels <- function(dfs, ref_df, type_pal, ref_name, subset, label = "cell_type",
                            proj_key = "refumap_", ref_key = "UMAP") {
  df     <- dfs[[subset]]
  p_ref  <- plot_ref_umap(ref_df, type_pal, ref_name, label, ref_key)
  p_proj <- plot_ref_proj_umap(df, ref_df, type_pal, label, proj_key, ref_key)
  cowplot::plot_grid(
    cowplot::plot_grid(p_ref, p_proj, nrow = 1),
    ref_make_legend(p_ref), nrow = 1, rel_widths = c(1, 0.22))
}

# Draw reference UMAP (left) + kernel-density of query cells (right).
ref_plot_density <- function(dfs, ref_df, type_pal, ref_name, subset, label = "cell_type",
                             proj_key = "refumap_", ref_key = "UMAP") {
  df        <- dfs[[subset]]
  p_ref     <- plot_ref_umap(ref_df, type_pal, ref_name, label, ref_key)
  p_density <- plot_query_density_umap(df, ref_df, proj_key, ref_key)
  cowplot::plot_grid(p_ref, p_density, nrow = 1)
}
ref_plot_density2 <- function(dfs, ref_df, type_pal, ref_name, subset, label = "cell_type",
                             proj_key = "refumap_", ref_key = "UMAP") {
  df        <- dfs[[subset]]
  p_density <- df %>%
    split(.,.$sample.name) %>%
    lapply(function(comp) plot_query_density_umap(comp, ref_df, proj_key, ref_key))

  patchwork::wrap_plots(p_density,nrow = 1)
}

# Draw reference UMAP (left) + query coloured by predicted cell type (right) + legend.
ref_plot_celltype <- function(dfs, ref_df, type_pal, ref_name, subset, label = "cell_type",
                              ref_key = "UMAP") {
  df      <- dfs[[subset]]
  p_ref   <- plot_ref_umap(ref_df, type_pal, ref_name, label, ref_key)
  p_query <- plot_query_label_umap(df, type_pal, label)
  cowplot::plot_grid(
    cowplot::plot_grid(p_ref, p_query, nrow = 1),
    ref_make_legend(p_ref), nrow = 1, rel_widths = c(1, 0.22))
}

# Draw query UMAP (left) + mapping-score UMAP (right).
ref_plot_mapping <- function(dfs, type_pal, subset, label = "cell_type") {
  df      <- dfs[[subset]]
  p_query <- plot_query_label_umap(df, type_pal, label)
  p_map   <- plot_mapping_score_umap(df)
  cowplot::plot_grid(p_query, p_map, nrow = 1)
}

# Draw mapping-score violins (by cluster + by cell type) for all compartments stacked.
# Compartments are separated by a bold section label.
# Returns a patchwork; returns NULL if no scores are available.
ref_plot_violins <- function(dfs, ref_name, label = "cell_type",
                              subsets    = c("hsc", "mnc", "stroma"),
                              compartment_map = c(hsc = "HSC", mnc = "MNC", stroma = "Stroma")) {
  blocks <- list()
  for (s in subsets) {
    df <- dfs[[s]]
    comp <- compartment_map[[s]]
    p1   <- plot_score_by_cluster(df, ref_name, comp)
    p2   <- plot_score_by_celltype(df, ref_name, comp, label)
    pair <- Filter(Negate(is.null), list(p1, p2))
    # Section label for this compartment
    label_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = comp,
                        size = 5, fontface = "bold", hjust = 0) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 10, b = 2))
    blocks <- c(blocks, list(label_plot), pair)
  }
  patchwork::wrap_plots(blocks, ncol = 1)
}

# Draw resolution heatmaps for all compartments (one per page).
# Renders directly via ComplexHeatmap::draw() — call inside a chunk, not inline.
# label: ref cell-type column; predicted column is paste0("predicted.", label).
ref_draw_heatmaps <- function(dfs, ref_name, label = "cell_type",
                               resolutions = c(0.1, 0.3, 0.5, 0.7),
                               subsets = c("hsc", "mnc", "stroma"),
                               compartment_map = c(hsc = "HSC", mnc = "MNC", stroma = "Stroma")) {
  pred_col <- paste0("predicted.", label)
  for (s in subsets) {
    df <- dfs[[s]]
    if (is.null(df) || !pred_col %in% colnames(df)) next
    comp <- compartment_map[[s]]
    ht   <- plot_resolution_heatmap(df, pred_col,
              resolutions = resolutions,
              title       = paste0(comp, " \u2014 ", ref_name, " \u00d7 resolution"))
    if (!is.null(ht)) ComplexHeatmap::draw(ht)
  }
}

# Build ARI line plot across compartments.
# Returns a ggplot; returns NULL if no ARI data are available.
# label: ref cell-type column; predicted column is paste0("predicted.", label).
ref_plot_ari <- function(dfs, ref_name, resolutions, label = "cell_type",
                          subsets = c("hsc", "mnc", "stroma"),
                          compartment_map = c(hsc = "HSC", mnc = "MNC", stroma = "Stroma")) {
  pred_col <- paste0("predicted.", label)
  ari_rows <- list()
  for (s in subsets) {
    df <- dfs[[s]]
    if (is.null(df) || !pred_col %in% colnames(df)) next
    ari_df <- compute_ari_table(df, resolutions, label)
    if (!all(is.na(ari_df$ari))) {
      ari_df$compartment <- compartment_map[[s]]
      ari_rows <- c(ari_rows, list(ari_df))
    }
  }
  if (length(ari_rows) == 0) return(invisible(NULL))
  combined <- do.call(rbind, ari_rows)
  ggplot2::ggplot(combined,
    ggplot2::aes(x = resolution, y = ari,
                 colour = compartment, group = compartment)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_brewer(palette = "Set1", name = "Compartment") +
    ggplot2::scale_x_continuous(breaks = resolutions) +
    ggplot2::labs(x = "Clustering resolution", y = "Adjusted Rand index",
                  title = paste0("ARI \u2014 ", ref_name)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# ── Dim-reduction grid helpers ────────────────────────────────────────────────

#' Extract pre-computed SCT dim-reduction data for the parameter grid notebook.
#'
#' Takes an already-processed Seurat object (SCTransform + PCA + umap.nn already
#' run by the pipeline), extracts the elbow, PCA scatter, and pipeline UMAP as
#' plain data frames, then computes one additional Harmony UMAP for the
#' PCA-vs-Harmony visual comparison section.
#'
#' @param srt          Processed Seurat object (mpn_processed_*).
#' @param n_dims       Number of PCA dims shown on the elbow plot.
#' @param k_umap       n.neighbors for the Harmony UMAP.
#' @param pca_name     Name of the PCA reduction (default "pca.SCT").
#' @return Named list: elbow, pca, umap_nn, umap_harmony — all plain data.frames.
compute_dim_grid_sct_pre <- function(srt, n_dims = 50, k_umap = 30,
                                      pca_name = "pca.SCT") {
  # Elbow data
  elbow    <- Seurat::ElbowPlot(srt, ndims = n_dims, reduction = pca_name)$data

  # PCA scatter (PC1 vs PC2, coloured by sample)
  pca_emb           <- as.data.frame(Seurat::Embeddings(srt, pca_name)[, 1:2])
  colnames(pca_emb) <- c("PC1", "PC2")
  pca_emb$sample    <- srt$sample.name

  # Pipeline UMAP (umap.nn — graph-based, identical to the analysis report)
  umap_nn           <- as.data.frame(Seurat::Embeddings(srt, "umap.nn")[, 1:2])
  colnames(umap_nn) <- c("UMAP1", "UMAP2")
  umap_nn$sample    <- srt$sample.name

  # Harmony UMAP — run on top of the existing PCA; for visual comparison only
  suppressMessages(invisible(capture.output({
    srt <- harmony::RunHarmony(srt,
                                group.by.vars  = "sample.name",
                                reduction      = pca_name,
                                reduction.save = "harmony",
                                verbose        = FALSE)
    srt <- Seurat::RunUMAP(srt,
                            dims           = seq(n_dims),
                            reduction      = "harmony",
                            n.neighbors    = k_umap,
                            reduction.name = "umap.harmony",
                            return.model   = FALSE,
                            verbose        = FALSE)
  })))
  umap_harmony           <- as.data.frame(Seurat::Embeddings(srt, "umap.harmony")[, 1:2])
  colnames(umap_harmony) <- c("UMAP1", "UMAP2")
  umap_harmony$sample    <- srt$sample.name

  list(elbow        = elbow,
       pca          = pca_emb,
       umap_nn      = umap_nn,
       umap_harmony = umap_harmony)
}

#' Run the dims × k grid search for one subset (SCT path).
#'
#' SCTransform has already been run on the input object by the pipeline.
#' This function only loops over FindNeighbors + FindClusters + RunUMAP,
#' so SCTransform is NOT repeated per grid cell.
#'
#' @param srt          Clean Seurat object (mpn_clean_* — pre-SCTransform).
#' @param ref_assay    Reference assay, default "SCT".
#' @param dims_grid    Integer vector of PCA dim counts to test.
#' @param k_grid       Integer vector of k_neighbors values to test.
#' @param res_scan     Numeric vector of clustering resolutions to scan.
#' @param n_features   Number of variable features for SCTransform.
#' @param k_umap       n.neighbors for UMAP.
#' @param exclude_genes Character vector of genes to exclude.
#' @return Named list of lightweight data.frames (UMAP + cluster cols + sample),
#'         one entry per "dims<n>_k<k>" combination.
compute_dim_grid_sct_runs <- function(srt, ref_assay = "SCT",
                                       dims_grid, k_grid, res_scan,
                                       n_features, k_umap,
                                       exclude_genes = character()) {
  # Run SCTransform + PCA once — reused across all grid cells
  suppressMessages(invisible(capture.output({
    srt_pre <- preprocess_seurat(srt,
                                  ref_assay     = ref_assay,
                                  n_features    = n_features,
                                  exclude_genes = exclude_genes)
  })))

  pca_name     <- paste0("pca.", ref_assay)
  clust_prefix <- paste0("sNN_pca.", ref_assay, "_")

  runs <- list()
  for (n_dims in dims_grid) {
    for (k in k_grid) {
      label <- paste0("dims", n_dims, "_k", k)
      suppressMessages(invisible(capture.output({
        srt_i <- find_neighbors_clusters(srt_pre,
                                          dims_i    = seq(n_dims),
                                          ref_assay = ref_assay,
                                          k         = k,
                                          res       = res_scan)
        srt_i <- Seurat::RunUMAP(srt_i,
                                  n.neighbors    = k_umap,
                                  nn.name        = Seurat::Neighbors(srt_i),
                                  reduction.name = "umap.nn",
                                  return.model   = FALSE,
                                  verbose        = FALSE)
      })))

      umap_emb           <- as.data.frame(Seurat::Embeddings(srt_i, "umap.nn")[, 1:2])
      colnames(umap_emb) <- c("UMAP1", "UMAP2")
      clust_cols         <- grep(paste0("^", clust_prefix),
                                  colnames(srt_i@meta.data), value = TRUE)
      run_df             <- cbind(umap_emb,
                                   srt_i@meta.data[, c("sample.name", clust_cols),
                                                    drop = FALSE])
      runs[[label]] <- run_df
      rm(srt_i); gc()
    }
  }

  rm(srt_pre); gc()
  runs
}

# ── Milo differential abundance ───────────────────────────────────────────────

#' Add condition label column to a data.frame with a sample.name column.
#'
#' @param df      data.frame with a sample.name column.
#' @param subset  One of "stroma", "hsc", "mnc".
#' @return df with an added `condition` character column.
add_condition_df <- function(df, subset, condition) {
  df$condition <- switch(subset,
    stroma = dplyr::case_when(
      grepl("WT",   df[[condition]]) ~ "WT",
      grepl("Het",  df[[condition]]) ~ "Het",
      TRUE ~ NA_character_),
    hsc = dplyr::case_when(
      grepl("_WT",    df[[condition]]) ~ "WT",
      grepl("mCh-",   df[[condition]]) ~ "Het_mChNeg",
      grepl("mCh\\+", df$sample.name) ~ "Het_mChPos",
      TRUE ~ NA_character_),
    mnc = dplyr::case_when(
      grepl("_WT$",   df[[condition]]) ~ "WT",
      grepl("mCh-",   df[[condition]]) ~ "mChNeg",
      grepl("mCh\\+", df[[condition]]) ~ "mChPos",
      TRUE ~ NA_character_)
  )
  df
}

#' Add condition column to a Seurat object.
#'
#' @param srt    Seurat object with a sample.name metadata column.
#' @param subset One of "stroma", "hsc", "mnc".
#' @return srt with a condition metadata column.
add_condition_seurat <- function(srt, subset, condition) {
  srt$condition <- add_condition_df(
    data.frame(sample.name = srt@meta.data[[condition]]),
    subset,
    condition)$condition
  srt
}

#' Build a Milo object from a Seurat object.
#'
#' Converts to SCE, copies the PCA and UMAP reductions, builds the
#' neighbourhood graph, and counts cells per sample.
#'
#' @param srt   Seurat object with condition and dm2_class columns.
#' @param k     Number of nearest neighbours for graph building.
#' @param d     Number of PCA dimensions to use.
#' @param prop  Proportion of cells to sample as neighbourhood index points.
#' @return Milo object ready for testNhoods.
.build_milo <- function(srt, subset, graph_name, k = 30, d = 30, prop = 0.1, condition, replicates, dim_reduc, umap) {
  #exclude cells with no clear source (can't define replicate)
  srt  <- srt[, !srt@meta.data[[replicates]] %in% c("negative")]
  srt <- add_condition_seurat(srt, subset = subset, condition  = condition)

  srt$milo_sample <- paste(srt@meta.data[[condition]], srt@meta.data[[replicates]], sep = "__")

  pca_mat <- Seurat::Embeddings(srt, dim_reduc)[, seq_len(d)]
  umap_mat <- Seurat::Embeddings(srt, umap)

  sce <- Seurat::as.SingleCellExperiment(srt)
  SingleCellExperiment::reducedDim(sce, dim_reduc) <- pca_mat
  SingleCellExperiment::reducedDim(sce, umap) <- umap_mat

  milo <- miloR::Milo(sce)
  #milo <- miloR::buildGraph(milo, k = k, d = d, reduced.dim = dim_reduc)
  graph_weighted <- ifelse(str_detect(graph_name, "sNN"), TRUE,NULL)
  miloR::graph(milo) <- igraph::graph_from_adjacency_matrix(srt[[graph_name]],
                                                            mode = "undirected",
                                                            weighted = graph_weighted)
  milo <- miloR::makeNhoods(milo, prop = prop, k = k, d = d,
                             refined = TRUE, reduced_dims = dim_reduc)
  milo <- miloR::buildNhoodGraph(milo)
  milo <- miloR::countCells(milo,
                            samples = "milo_sample",
                            meta.data = as.data.frame(SummarizedExperiment::colData(milo)))
  milo
}

#' Run Milo DA for two conditions (cond_test vs cond_ref).
#'
#' Subsets to the two relevant conditions, builds the Milo object, and runs
#' testNhoods with a two-group design.
#'
#' @param srt_all   Seurat object with condition column; singlets only.
#' @param cond_test Test condition label.
#' @param cond_ref  Reference condition label.
#' @param k,d       Neighbours / PCA dims passed to .build_milo.
#' @return List: milo, da_results (data.frame), nhood_sizes (data.frame),
#'         cluster_col (character — name of the sNN cluster column).
.run_milo_pair <- function(srt_all, cond_test, cond_ref,
                            k = 30, d = 30, cluster_col) {
  srt2 <- srt_all[, srt_all$condition %in% c(cond_ref, cond_test)]
  srt2$condition <- factor(srt2$condition, levels = c(cond_ref, cond_test))

  milo <- .build_milo(srt2, k = k, d = d)

  design_df <- as.data.frame(SummarizedExperiment::colData(milo)) %>%
    dplyr::select(milo_sample, condition) %>%
    dplyr::distinct()
  rownames(design_df) <- design_df$milo_sample

  # Check that every condition has >= 2 samples — edgeR requires full rank.
  # This can fail when HTO demultiplexing leaves only one tag for a condition
  # after filtering negatives and multiplets.
  samples_per_cond <- table(design_df$condition)
  insufficient <- names(samples_per_cond)[samples_per_cond < 2]
  if (length(insufficient) > 0) {
    warning(
      "Milo skipped for contrast ", cond_test, " vs ", cond_ref, ": ",
      paste(insufficient, collapse = ", "),
      " has only ", paste(samples_per_cond[insufficient], collapse = ", "),
      " sample(s) — design matrix would be rank-deficient."
    )
    return(NULL)
  }

  da <- miloR::testNhoods(milo,
                           design        = ~ condition,
                           design.df     = design_df,
                           reduced.dim   = "PCA",
                           fdr.weighting = "graph-overlap")
  da <- miloR::annotateNhoods(milo, da, coldata_col = cluster_col)

  nhood_sizes <- data.frame(
    Nhood = seq_len(ncol(miloR::nhoods(milo))),
    Size  = colSums(miloR::nhoods(milo) > 0)
  )

  list(milo        = milo,
       da_results  = da,
       nhood_sizes = nhood_sizes,
       cluster_col = cluster_col)
}

#' Run Milo DA for a single two-condition comparison (e.g. Stroma: WT vs Het).
#'
#' Removes negative/multiplet cells, adds condition labels, builds the Milo
#' object, and runs testNhoods.
#'
#' @param srt        Processed Seurat object.
#' @param subset     One of "stroma", "hsc", "mnc" — passed to add_condition_seurat.
#' @param cluster_col Name of the cluster column in metadata (sNN_pca.* pattern).
#' @param k,d,prop   Milo graph parameters.
#' @return List: milo, da_results, nhood_sizes, cluster_col.
run_milo <- function(srt, subset, cluster_col, k = 30, d = 30, prop = 0.1) {
  srt  <- add_condition_seurat(srt, subset)
  srt  <- srt[, !srt$dm2_class %in% c("negative", "multiplet")]
  milo <- .build_milo(srt, k = k, d = d, prop = prop)

  design_df <- as.data.frame(SummarizedExperiment::colData(milo)) %>%
    dplyr::select(milo_sample, condition) %>%
    dplyr::distinct()
  rownames(design_df) <- design_df$milo_sample

  da <- miloR::testNhoods(milo,
                           design        = ~ condition,
                           design.df     = design_df,
                           reduced.dim   = "PCA",
                           fdr.weighting = "graph-overlap")
  da <- miloR::annotateNhoods(milo, da, coldata_col = cluster_col)

  nhood_sizes <- data.frame(
    Nhood = seq_len(ncol(miloR::nhoods(milo))),
    Size  = colSums(miloR::nhoods(milo) > 0)
  )

  list(milo        = milo,
       da_results  = da,
       nhood_sizes = nhood_sizes,
       cluster_col = cluster_col)
}

#
format_contrast <- function(test_vars){
  attach(test_vars)
  contrast_string <- glue::glue("{test_variable}{test_condition} - {test_variable}{test_control}")
  detach(test_vars)
  return(contrast_string)
}

#
test_milo_pair <- function(test_vars,
                           milo,
                           reduced_dim = "pca.SCT",
                           cluster_col,
                           norm_method = "TMM",
                           condition){
  design_df <- as.data.frame(SummarizedExperiment::colData(milo)) %>%
    dplyr::select(milo_sample, condition) %>%
    dplyr::distinct()
  rownames(design_df) <- design_df$milo_sample
  contrast <- format_contrast(test_vars)
  cat(paste0(contrast, "\n"))

  da <- miloR::testNhoods(milo,
                          design        = ~0+condition,
                          design.df     = design_df,
                          model.contrasts = contrast,
                          reduced.dim   = reduced_dim,
                          fdr.weighting = "graph-overlap",
                          norm.method = norm_method)

  da <- miloR::annotateNhoods(milo, da, coldata_col = cluster_col)

  nhood_sizes <- data.frame(
    Nhood = seq_len(ncol(miloR::nhoods(milo))),
    Size  = colSums(miloR::nhoods(milo) > 0)
  )

  rm(milo); gc()

  out <- list(
       da_results  = da,
       nhood_sizes = nhood_sizes,
       cluster_col = cluster_col,
       model_contrast = contrast)

  return(out)
}


#' Run three pairwise Milo contrasts: neg vs ref, pos vs ref, pos vs neg.
#'
#' Runs .run_milo_pair independently for all three contrasts.
#' NULL results (insufficient replication) are silently dropped.
#'
#' @param srt        Processed Seurat object.
#' @param subset     One of "stroma", "hsc", "mnc".
#' @param cluster_col Name of the cluster column in metadata.
#' @param cond_neg,cond_pos,cond_ref Condition labels.
#' @param k,d        Milo graph parameters.
#' @return List: res_neg, res_pos, res_pos_neg — each the return value of
#'   .run_milo_pair (or absent if that contrast failed).
run_milo_contrasts <- function(srt, subset, cluster_col,
                                cond_neg, cond_pos, cond_ref = "WT",
                                k = 30, d = 30) {
  srt <- add_condition_seurat(srt, subset)
  srt <- srt[, !srt$dm2_class %in% c("negative", "multiplet")]

  res_neg     <- .run_milo_pair(srt, cond_test = cond_neg, cond_ref = cond_ref,
                                k = k, d = d, cluster_col = cluster_col)
  res_pos     <- .run_milo_pair(srt, cond_test = cond_pos, cond_ref = cond_ref,
                                k = k, d = d, cluster_col = cluster_col)
  res_pos_neg <- .run_milo_pair(srt, cond_test = cond_pos, cond_ref = cond_neg,
                                k = k, d = d, cluster_col = cluster_col)

  available <- list(res_neg = res_neg, res_pos = res_pos, res_pos_neg = res_pos_neg)
  available <- Filter(Negate(is.null), available)
  if (length(available) == 0)
    stop("All Milo contrasts failed — no results to return.")
  available
}

# ── Cluster condition purity ────────────────────────────────────────────────

#' Compute per-cluster condition purity across multiple resolutions.
#'
#' For each resolution column in \code{recluster_df}, and for each cluster at
#' that resolution, computes:
#' \itemize{
#'   \item \code{purity}: fraction of cells belonging to the dominant condition.
#'   \item \code{dominant_condition}: the condition with the most cells.
#'   \item \code{n_cells}: cluster size.
#' }
#' Also returns \code{weighted_purity}: the size-weighted mean purity across all
#' clusters at each resolution — useful as a scalar summary to plot vs resolution.
#'
#' @param recluster_df  data.frame from \code{recluster_resolutions_df}:
#'   rows = cells, columns = cell_id + res_* cluster assignments.
#' @param srt           Seurat object with \code{sample.name} metadata.
#' @param subset        One of "stroma", "hsc", "mnc" — passed to
#'   \code{add_condition_df} to derive condition labels.
#' @return data.frame with columns: resolution, cluster, dominant_condition,
#'   purity, n_cells; plus attribute \code{"weighted_purity"} (data.frame:
#'   resolution, weighted_purity).
compute_condition_purity <- function(recluster_df, srt, subset, condition) {
  # Derive condition for each cell
  cond_df <- add_condition_df(
    data.frame(cell_id     = colnames(srt),
               sample.name = srt$sample.name,
               stringsAsFactors = FALSE),
    subset,
    condition
  )

  res_cols <- grep("^res_", colnames(recluster_df), value = TRUE)

  purity_rows <- purrr::map_dfr(res_cols, function(col) {
    res_val <- as.numeric(gsub("_", ".", gsub("^res_", "", col)))

    df <- recluster_df[, c("cell_id", col)]
    colnames(df)[2] <- "cluster"
    df <- dplyr::left_join(df, cond_df[, c("cell_id", "condition")],
                           by = "cell_id")
    df <- df[!is.na(df$condition), ]

    df %>%
      dplyr::group_by(cluster, condition) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        dominant_condition = condition[which.max(n)],
        purity             = max(n) / sum(n),
        n_cells            = sum(n),
        .groups            = "drop"
      ) %>%
      dplyr::mutate(resolution = res_val)
  })

  wp <- purity_rows %>%
    dplyr::group_by(resolution) %>%
    dplyr::summarise(
      weighted_purity = sum(purity * n_cells) / sum(n_cells),
      n_clusters      = dplyr::n(),
      .groups         = "drop"
    )

  attr(purity_rows, "weighted_purity") <- wp
  purity_rows[, c("resolution", "cluster", "dominant_condition", "purity", "n_cells")]
}

#' Assign a condition-purity group label to every cell.
#'
#' At a given resolution, clusters whose purity exceeds \code{purity_cutoff}
#' are labelled by their dominant condition; the rest are labelled "Impure".
#' Returns a data.frame with \code{cell_id}, \code{purity_group}, \code{UMAP1},
#' \code{UMAP2} — ready to load directly in a notebook.
#'
#' @param condition_purity  Output of \code{compute_condition_purity}.
#' @param recluster_df      Output of \code{recluster_resolutions_df} (cell_id + res_* cols).
#' @param srt               Seurat object (for UMAP embedding).
#' @param purity_res        Resolution to evaluate (numeric). Default 10.
#' @param purity_cutoff     Purity threshold above which a cluster is condition-pure. Default 0.6.
#' @param reduction         Name of the UMAP reduction. Default \code{"umap.nn"}.
#' @return data.frame: cell_id, purity_group, UMAP1, UMAP2.
assign_purity_groups <- function(condition_purity, recluster_df, srt,
                                  purity_res    = 10,
                                  purity_cutoff = 0.6,
                                  reduction     = "umap.nn") {
  res_col <- paste0("res_", gsub("\\.", "_", purity_res))

  cluster_label_map <- condition_purity %>%
    dplyr::filter(resolution == purity_res) %>%
    dplyr::mutate(
      purity_group = ifelse(purity >= purity_cutoff, dominant_condition, "Impure")
    ) %>%
    dplyr::select(cluster, purity_group) %>%
    tibble::deframe()

  umap_emb <- as.data.frame(Seurat::Embeddings(srt, reduction))
  colnames(umap_emb) <- c("UMAP1", "UMAP2")
  umap_emb$cell_id <- rownames(umap_emb)

  recluster_df %>%
    dplyr::select(cell_id, dplyr::all_of(res_col)) %>%
    dplyr::mutate(
      purity_group = cluster_label_map[as.character(.data[[res_col]])]
    ) %>%
    dplyr::select(cell_id, purity_group) %>%
    dplyr::left_join(umap_emb, by = "cell_id")
}

# ── SuperCell helpers ───────────────────────────────────────────────────────

#' Compute SuperCell supercells from a processed Seurat object.
#'
#' Runs \code{SuperCell::SCimplify_from_embedding} on the Seurat PCA embedding,
#' splitting cells by \code{sample.name} so no supercell spans multiple samples.
#'
#' @param srt        Processed Seurat object (must have \code{pca_name} reduction
#'                   and \code{sample.name} metadata column).
#' @param gamma      Graining factor (target cells per supercell). Default 10.
#' @param k.knn      k for the kNN graph. Default 5.
#' @param n.pc       Number of PCs to use from the Seurat reduction. Default 30.
#' @param pca_name   Name of the PCA reduction in \code{srt}. Default \code{"pca.SCT"}.
#' @return SCimplify result list (graph.supercells, membership, supercell_size, …).
compute_supercells <- function(srt,
                               gamma      = 20,
                               k.knn      = 5,
                               n.pc       = 30,
                               pca_name   = "pca.SCT",
                               split_by   = "sample.name",
                               fast.merge = FALSE) {
  pca_emb   <- Seurat::Embeddings(srt, pca_name)[, seq_len(n.pc), drop = FALSE]
  split_vec <- if (is.null(split_by)) NULL else srt@meta.data[[split_by]]
  SuperCell::SCimplify_from_embedding(
    X                    = pca_emb,
    k.knn                = k.knn,
    gamma                = gamma,
    fast.merge           = fast.merge,
    cell.split.condition = split_vec
  )
}

#' Extract per-supercell metadata from a SuperCell result.
#'
#' Returns one row per supercell with size, majority cluster, median UMAP
#' coordinates, and sample / condition labels.
#'
#' @param sc          SCimplify result from \code{compute_supercells}.
#' @param srt         Seurat object used to build \code{sc}.
#' @param cluster_col Name of the cluster metadata column in \code{srt}.
#' @return data.frame with columns: supercell, sample, size,
#'   majority_cluster, UMAP1, UMAP2.
extract_supercell_meta <- function(sc, srt, cluster_col) {
  umap_emb   <- Seurat::Embeddings(srt, "umap.nn")
  membership <- sc$membership          # named by cell barcode → supercell id

  cell_df <- data.frame(
    cell_id   = names(membership),
    supercell = as.integer(membership),
    cluster   = srt@meta.data[names(membership), cluster_col],
    sample    = srt@meta.data[names(membership), "sample.name"],
    UMAP1     = umap_emb[names(membership), 1],
    UMAP2     = umap_emb[names(membership), 2],
    stringsAsFactors = FALSE
  )

  # Each supercell is pure (single sample) because of cell.split.condition.
  # group_by supercell + sample captures that without ambiguity.
  cell_df %>%
    dplyr::group_by(supercell, sample) %>%
    dplyr::summarise(
      size             = dplyr::n(),
      majority_cluster = names(which.max(table(cluster))),
      UMAP1            = median(UMAP1),
      UMAP2            = median(UMAP2),
      .groups = "drop"
    )
}

#' Build a supercell gene-expression matrix (genes × supercells).
#'
#' Wraps \code{SuperCell::supercell_GE}, pulling the normalised data layer from
#' the Seurat object so the matrix is on the same scale used for PCA.
#' Column names are formatted as \code{<sample>_SC<id>} for readability.
#'
#' @param sc    SCimplify result from \code{compute_supercells}.
#' @param srt   Seurat object used to build \code{sc}.
#' @param assay Assay to pull the data layer from. Default \code{"SCT"}.
#' @return Sparse matrix (genes × supercells).
extract_supercell_ge <- function(sc, srt, meta, assay = "SCT") {
  ge  <- Seurat::GetAssayData(srt, assay = assay, layer = "data")
  mat <- SuperCell::supercell_GE(ge, sc$membership)
  # Label columns as <sample>_SC<id> using supercell_meta for sample lookup
  sc_sample <- setNames(meta$sample, meta$supercell)
  colnames(mat) <- paste0(sc_sample[colnames(mat)], "_SC", colnames(mat))
  mat
}

#' Plot Milo neighbourhood graph coloured by DA logFC
#'
#' Drop-in replacement for \code{miloR::plotNhoodGraphDA} that merges the
#' internal \code{plotNhoodGraph} call and exposes all visual parameters.
#' Default values reproduce the original package output exactly.
#'
#' @param x          Milo object (with a built neighbourhood graph).
#' @param milo_res   DA results data.frame from \code{miloR::testNhoods}.
#' @param alpha      FDR threshold; nhoods with SpatialFDR > alpha are zeroed. Default 0.05.
#' @param res_column Column in \code{milo_res} used for fill colour. Default "logFC".
#' @param layout     Reduced dimension name to use as node coordinates. Default "UMAP".
#' @param subset.nhoods Integer vector of nhood indices to plot. Default NULL (all).
#' @param size_range \code{c(min, max)} for node point size. Default \code{c(0.5, 3)}.
#' @param node_stroke  Stroke width for node shape 21. Default 0.3.
#' @param node_colour  Stroke colour for node shape 21. Default "black".
#' @param node_alpha  Alpha for significant nodes (|logFC| > 0). Default 1.
#' @param nonsig_alpha Alpha for non-significant nodes (|logFC| == 0). Default 1.
#'   Set < 1 (e.g. 0.1) to visually suppress non-DA nhoods like \code{highlight.da = TRUE}.
#' @param edge_colour Colour of graph edges. Default "grey66".
#' @param edge_alpha  Alpha of graph edges. Default 0.2.
#' @param edge_width_range \code{c(min, max)} for edge width (overlap size). Default \code{c(0.2, 3)}.
#' @param low_col    Fill colour for negative logFC. Default "#2166AC" (blue).
#' @param mid_col    Fill colour for logFC = 0. Default "white".
#' @param high_col   Fill colour for positive logFC. Default "#B2182B" (red).
#' @param base_size  Base font size passed to \code{theme_classic}. Default 14.
#' @return A ggraph/ggplot object.
plot_nhood_graph_da <- function(x,
                                milo_res,
                                alpha            = 0.05,
                                res_column       = "logFC",
                                layout           = "UMAP",
                                subset.nhoods    = NULL,
                                size_range       = c(0.5, 3),
                                node_stroke      = 0.3,
                                node_colour      = "black",
                                node_alpha       = 1,
                                nonsig_alpha     = 1,
                                edge_colour      = "grey66",
                                edge_alpha       = 0.2,
                                edge_width_range = c(0.2, 3),
                                low_col          = "#2166AC",
                                mid_col          = "white",
                                high_col         = "#B2182B",
                                base_size        = 14) {
  # ── 1. Apply FDR threshold (mirrors plotNhoodGraphDA) ──────────────────────
  signif_res <- milo_res
  signif_res$SpatialFDR[is.na(signif_res$SpatialFDR)] <- 1
  signif_res[signif_res$SpatialFDR > alpha, res_column] <- 0

  SummarizedExperiment::colData(x)[res_column] <- NA

  if (!is.null(subset.nhoods)) {
    sub.indices <- miloR::nhoodIndex(x)[subset.nhoods]
    SummarizedExperiment::colData(x)[unlist(sub.indices[signif_res$Nhood]), res_column] <-
      signif_res[, res_column]
  } else {
    SummarizedExperiment::colData(x)[unlist(miloR::nhoodIndex(x)[signif_res$Nhood]), res_column] <-
      signif_res[, res_column]
  }

  # ── 2. Build/subset graph (mirrors plotNhoodGraph is.da = TRUE branch) ─────
  nh_graph <- miloR::nhoodGraph(x)

  # Stamp logFC values onto graph vertices so permute() and colour mapping work
  nh_graph <- igraph::set_vertex_attr(nh_graph, name = res_column,
                                      value = signif_res[, res_column])

  if (!is.null(subset.nhoods)) {
    nh_graph <- igraph::induced_subgraph(
      nh_graph,
      vids = which(as.numeric(igraph::V(nh_graph)$name) %in%
                     unlist(miloR::nhoodIndex(x)[subset.nhoods]))
    )
  }

  # Sort so largest |logFC| nhoods are drawn last (on top)
  nh_graph <- igraph::permute(
    nh_graph,
    order(abs(igraph::vertex_attr(nh_graph)[[res_column]]), decreasing = TRUE)
  )

  # Layout: pull coordinates from the Milo reducedDim
  layout_mat <- SingleCellExperiment::reducedDim(x, layout)[
    as.numeric(igraph::vertex_attr(nh_graph)$name), , drop = FALSE]
  if (!is.matrix(layout_mat)) layout_mat <- as.matrix(layout_mat)

  # Colour column
  col_vals <- SummarizedExperiment::colData(x)[
    as.numeric(igraph::vertex_attr(nh_graph)$name), res_column]
  igraph::V(nh_graph)$colour_by <- col_vals

  # Graph layout data.frame (for filtering significant nodes)
  g.df <- ggraph::create_layout(nh_graph, layout = layout_mat)
  g.df[[res_column]] <- igraph::vertex_attr(nh_graph)[[res_column]]

  # ── 3. Build plot ───────────────────────────────────────────────────────────
  pl <- ggraph::ggraph(igraph::simplify(nh_graph), layout = layout_mat) +
    ggraph::geom_edge_link0(
      ggplot2::aes(width = weight),
      edge_colour = edge_colour,
      edge_alpha  = edge_alpha
    ) +
    # All nodes (non-significant ones at nonsig_alpha)
    ggraph::geom_node_point(
      ggplot2::aes(fill = colour_by, size = size),
      colour = node_colour,
      alpha  = nonsig_alpha,
      shape  = 21,
      stroke = node_stroke
    ) +
    # Significant nodes drawn on top at node_alpha
    ggraph::geom_node_point(
      data    = g.df[abs(g.df[[res_column]]) > 0, , drop = FALSE],
      mapping = ggplot2::aes(fill = colour_by, size = size),
      colour  = node_colour,
      alpha   = node_alpha,
      shape   = 21,
      stroke  = node_stroke
    ) +
    ggraph::scale_edge_width(range = edge_width_range, name = "overlap size") +
    ggplot2::scale_size(range = size_range, name = "Nhood size") +
    ggplot2::scale_fill_gradient2(
      low    = low_col,
      mid    = mid_col,
      high   = high_col,
      name   = res_column
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line  = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  pl
}

plot_nhood_graph_groups <- function(x,
                                    nhood_groups,
                                    colour_by,
                                    layout           = "UMAP",
                                    subset.nhoods    = NULL,
                                    size_range       = c(0.5, 3),
                                    node_stroke      = 0.3,
                                    node_colour      = "black",
                                    node_alpha       = 1,
                                    edge_colour      = "grey66",
                                    edge_alpha       = 0.2,
                                    edge_width_range = c(0.2, 3),
                                    palette          = NULL,
                                    base_size        = 14) {
  # ── 1. Stamp group values onto cell slots via nhoodIndex ────────────────────
  # nhood_groups must have a Nhood column (1-indexed) and a colour_by column.
  # Each nhood's centroid cell gets the group label; unassigned nhoods → NA.
  SummarizedExperiment::colData(x)[[colour_by]] <- NA_character_
  nhood_index <- miloR::nhoodIndex(x)
  cell_idx    <- unlist(nhood_index[nhood_groups$Nhood])
  SummarizedExperiment::colData(x)[cell_idx, colour_by] <-
    as.character(rep(nhood_groups[[colour_by]],
                     lengths(nhood_index[nhood_groups$Nhood])))

  # ── 2. Build / subset graph ─────────────────────────────────────────────────
  nh_graph <- miloR::nhoodGraph(x)

  if (!is.null(subset.nhoods)) {
    nh_graph <- igraph::induced_subgraph(
      nh_graph,
      vids = which(as.numeric(igraph::V(nh_graph)$name) %in%
                     unlist(nhood_index[subset.nhoods]))
    )
  }

  # ── 3. Layout ───────────────────────────────────────────────────────────────
  layout_mat <- SingleCellExperiment::reducedDim(x, layout)[
    as.numeric(igraph::vertex_attr(nh_graph)$name), , drop = FALSE]
  if (!is.matrix(layout_mat)) layout_mat <- as.matrix(layout_mat)

  # ── 4. Attach colour values to vertices ────────────────────────────────────
  igraph::V(nh_graph)$colour_by <- as.character(
    SummarizedExperiment::colData(x)[
      as.numeric(igraph::vertex_attr(nh_graph)$name), colour_by]
  )

  # ── 5. Discrete palette ─────────────────────────────────────────────────────
  grp_levels <- sort(unique(na.omit(igraph::V(nh_graph)$colour_by)))
  if (is.null(palette)) {
    palette <- setNames(scales::hue_pal()(length(grp_levels)), grp_levels)
  }

  # ── 6. Plot ─────────────────────────────────────────────────────────────────
  ggraph::ggraph(igraph::simplify(nh_graph), layout = layout_mat) +
    ggraph::geom_edge_link0(
      ggplot2::aes(width = weight),
      edge_colour = edge_colour,
      edge_alpha  = edge_alpha
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = colour_by, size = size),
      colour = node_colour,
      alpha  = node_alpha,
      shape  = 21,
      stroke = node_stroke
    ) +
    ggraph::scale_edge_width(range = edge_width_range, name = "overlap size") +
    ggplot2::scale_size(range = size_range, name = "Nhood size") +
    ggplot2::scale_fill_manual(values = palette, name = colour_by,
                               na.value = "grey80") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line  = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
}

# ── DE notebook helpers ─────────────────────────────────────────────────────

#' Build a scatter data.frame from DElegate all-vs-all DE results.
#'
#' For two-comparison subsets (HSC, MNC): inner-joins the neg and pos volcano
#' results, retaining all genes (sig and non-sig) for background plotting.
#' Returns columns: feature, cluster, log_fc_neg, log_fc_pos, sig_neg, sig_pos,
#' max_rate_neg, max_rate_pos, max_rate.
#'
#' For single-comparison subsets (Stroma, group_pos = NULL): filters and flags
#' a single volcano.  Returns columns: feature, cluster, log_fc, sig, max_rate.
#'
#' @param markers_df  DElegate all-vs-all output (data.frame).
#' @param group_neg   group1 label for the neg/single comparison.
#' @param group_ref   group2 label (reference).
#' @param group_pos   group1 label for the pos comparison; NULL for single-comparison.
#' @param padj_thresh Adjusted p-value threshold (default 0.001).
#' @param lfc_thresh  log FC threshold; up: log_fc >= lfc_thresh & rate1 >= rate_cutoff,
#'                    down: log_fc <= -lfc_thresh & rate2 >= rate_cutoff (default 1).
#' @param rate_cutoff Minimum expression rate in the higher-expressing group (default 0.2).
#' @return data.frame.
make_scatter_df <- function(markers_df,
                            group_neg,
                            group_ref,
                            group_pos   = NULL,
                            padj_thresh = 0.001,
                            lfc_thresh  = 1,
                            rate_cutoff = 0.2) {
  .sig <- function(padj, log_fc, rate1, rate2)
    (padj <= padj_thresh & log_fc >=  lfc_thresh & rate1 >= rate_cutoff) |
    (padj <= padj_thresh & log_fc <= -lfc_thresh & rate2 >= rate_cutoff)

  if (is.null(group_pos)) {
    markers_df %>%
      dplyr::filter(group1 == group_neg, group2 == group_ref) %>%
      dplyr::mutate(
        sig      = .sig(padj, log_fc, rate1, rate2),
        max_rate = pmax(rate1, rate2)
      )
  } else {
    volcano_neg <- markers_df %>%
      dplyr::filter(group1 == group_neg, group2 == group_ref) %>%
      dplyr::mutate(sig      = .sig(padj, log_fc, rate1, rate2),
                    max_rate = pmax(rate1, rate2))
    volcano_pos <- markers_df %>%
      dplyr::filter(group1 == group_pos, group2 == group_ref) %>%
      dplyr::mutate(sig      = .sig(padj, log_fc, rate1, rate2),
                    max_rate = pmax(rate1, rate2))
    dplyr::inner_join(
      dplyr::select(volcano_neg, feature, cluster,
                    log_fc_neg = log_fc, sig_neg = sig, max_rate_neg = max_rate),
      dplyr::select(volcano_pos, feature, cluster,
                    log_fc_pos = log_fc, sig_pos = sig, max_rate_pos = max_rate),
      by = c("feature", "cluster")
    ) %>%
      dplyr::mutate(max_rate = pmax(max_rate_neg, max_rate_pos))
  }
}


#' Classify scatter genes into named genesets.
#'
#' For two-comparison scatter_df (has log_fc_neg / log_fc_pos): assigns quadrant
#' (RU/RD/LU/LD) and diagonal band (shared / mCh+-biased / mCh--biased);
#' geneset name = "<quadrant>_<band>".
#'
#' For single-comparison scatter_df (has log_fc): assigns "up_Het" / "down_Het".
#'
#' @param scatter_df  Output of make_scatter_df().
#' @param lfc_geneset Minimum |log FC| for geneset membership (default 1).
#' @param min_size    Minimum genes to keep a geneset (default 3).
#' @return List with: genesets_df (data.frame with geneset column),
#'         genesets (named list of character vectors), universe (character vector).
make_genesets <- function(scatter_df, lfc_geneset = 1, min_size = 3) {
  two_comp <- "log_fc_neg" %in% colnames(scatter_df)
  universe <- unique(scatter_df$feature)

  if (two_comp) {
    genesets_df <- scatter_df %>%
      dplyr::filter(abs(log_fc_neg) >= lfc_geneset | abs(log_fc_pos) >= lfc_geneset) %>%
      dplyr::mutate(
        quadrant = paste0(dplyr::if_else(log_fc_neg >= 0, "R", "L"),
                          dplyr::if_else(log_fc_pos >= 0, "U", "D")),
        band = dplyr::case_when(
          log_fc_pos - log_fc_neg >  1 ~ "mCh\u207a-biased",
          log_fc_neg - log_fc_pos >  1 ~ "mCh\u2212-biased",
          TRUE                         ~ "shared"
        ),
        geneset = paste0(quadrant, "_", band)
      )
  } else {
    genesets_df <- scatter_df %>%
      dplyr::filter(sig, abs(log_fc) >= lfc_geneset) %>%
      dplyr::mutate(geneset = dplyr::if_else(log_fc > 0, "up_Het", "down_Het"))
  }

  genesets <- lapply(split(genesets_df$feature, genesets_df$geneset), unique)
  genesets <- genesets[lengths(genesets) >= min_size]

  list(genesets_df = genesets_df, genesets = genesets, universe = universe)
}

#' Run ORA (fgsea::fora) for each geneset against MSigDB collections.
#'
#' Downloads Hallmark, KEGG, Reactome, WikiPathways, GO:BP for the given species
#' and runs fora() per geneset.
#'
#' @param genesets  Named list of gene vectors (from compute_allvsall_genesets).
#' @param universe  Background gene vector (from compute_allvsall_genesets).
#' @param species   Species name passed to msigdbr (default "Mus musculus").
#' @param min_size,max_size fora() pathway size filters.
#' @return Named list of data.frames (one per geneset, with added geneset column).
compute_ora_results <- function(genesets,
                                 universe,
                                 species  = "Mus musculus",
                                 min_size = 5,
                                 max_size = 500) {
  msig_collections <- list(
    list(collection = "H",  subcollection = NULL),
    list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
    list(collection = "C2", subcollection = "CP:REACTOME"),
    list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
    list(collection = "C5", subcollection = "GO:BP")
  )
  msig_parts <- list()
  for (col in msig_collections) {
    msig_parts[[length(msig_parts) + 1]] <- msigdbr::msigdbr(
      species       = species,
      collection    = col$collection,
      subcollection = col$subcollection)
  }
  msig_df        <- do.call(rbind, msig_parts)
  gene_sets_msig <- split(msig_df$gene_symbol, msig_df$gs_name)
  rm(msig_parts, msig_df)

  ora_list <- list()
  for (gs_name in names(genesets)) {
    genes <- genesets[[gs_name]]
    if (length(genes) < min_size) next
    res             <- fgsea::fora(pathways = gene_sets_msig,
                                   genes    = genes,
                                   universe = universe,
                                   minSize  = min_size,
                                   maxSize  = max_size)
    res$geneset     <- gs_name
    ora_list[[gs_name]] <- res
  }
  ora_list
}

#' Winsorise a numeric vector by clipping to quantile bounds.
#'
#' Values below the lower quantile are set to that quantile; values above the
#' upper quantile are set to that quantile. NAs are preserved.
#'
#' Useful before Pearson correlation or mutual information to reduce the
#' influence of extreme outliers without discarding observations.
#'
#' @param x     Numeric vector.
#' @param probs Two-element numeric vector giving the lower and upper quantile
#'              probabilities. Default c(0.025, 0.975).
#' @return Numeric vector the same length as x, with tails clipped.
winsorize <- function(x, probs = c(0.025, 0.975)) {
  q <- quantile(x, probs, na.rm = TRUE)
  pmax(pmin(x, q[2]), q[1])
}

#' Correlate gene expression features to a numeric covariate.
#'
#' Standalone reimplementation of Seurat::GroupCorrelation with configurable
#' association measures. Returns a plain data.frame instead of a mutated Seurat
#' object, and supports classical correlation, rank-based, non-linear, and
#' information-theoretic measures.
#'
#' Cells with NA in the covariate are silently dropped.
#' Features with zero mean across retained cells are excluded.
#'
#' Method groups and required packages:
#'   Classical (fast, vectorised via stats::cor):
#'     "pearson"  — linear correlation
#'     "spearman" — rank-based, monotone relationships
#'     "kendall"  — rank-based, more robust to outliers
#'
#'   Extended (per-feature loop via {correlation}; require: correlation, energy):
#'     "distance"  — distance correlation; detects non-linear / non-monotone
#'                   dependence; 0 iff independent (energy::dcov internally)
#'     "biweight"  — biweight midcorrelation; robust Pearson alternative,
#'                   down-weights extreme cells (useful for sparse scRNA data)
#'     "hoeffding" — Hoeffding's D; powerful omnibus independence test
#'
#'   Information-theoretic (per-feature; requires: infotheo):
#'     "mutual_information" — MI after equal-frequency binning; captures any
#'                            dependency including non-monotone; scale-free
#'
#'   Geometric:
#'     "cosine" — cosine similarity between each feature vector and the
#'                covariate; equivalent to Pearson on uncentred vectors;
#'                reflects both direction and expression magnitude
#'
#' @param expr_mat  Numeric / sparse matrix (features × cells). Accepts
#'                  dgCMatrix. E.g. srt[["SCT"]]$data or GetAssayData(srt).
#' @param covariate A numeric vector (length == ncol(expr_mat)) OR a named list
#'                  of such vectors to compute associations for multiple
#'                  covariates in one call. When a list is supplied the returned
#'                  data.frame gains a `covariate` column with the list names.
#' @param method    Association method (see groups above). Default "pearson".
#' @param features  Optional character vector to restrict to a subset of genes.
#' @param p_value   Logical. Return p-values? Works for all methods:
#'                   - pearson/spearman/kendall: stats::cor.test()
#'                   - distance/biweight/hoeffding: correlation::cor_test()
#'                   - mutual_information: chi-squared test on the
#'                     contingency table of discretised bins
#'                   - cosine: t-test derived from cosine as correlation
#'                  Adds per-feature loop overhead. Default FALSE.
#' @param mi_bins   Number of bins for mutual_information discretisation.
#'                  NULL (default) uses infotheo::discretize() auto-selection
#'                  (equal-frequency, nbins = ceil(n^(1/3))).
#' @return A data.frame sorted by descending |association| with columns:
#'   feature, association, and optionally p_value. When covariate is a list,
#'   an additional `covariate` column identifies the source covariate.
correlate_features_to_covariate <- function(expr_mat,
                                            covariate,
                                            method   = c("pearson", "spearman", "kendall",
                                                         "distance", "biweight", "hoeffding",
                                                         "mutual_information", "cosine"),
                                            features = NULL,
                                            p_value  = FALSE,
                                            mi_bins  = NULL) {
  method <- match.arg(method)

  # -- Multi-covariate dispatch -------------------------------------------------
  if (is.list(covariate)) {
    if (is.null(names(covariate)))
      stop("When covariate is a list it must be named.")
    # Materialise sparse → dense once here so every per-covariate call reuses
    # the same dense matrix instead of triggering a separate coercion each time.
    if (inherits(expr_mat, "sparseMatrix"))
      expr_mat <- as.matrix(expr_mat)
    results <- lapply(names(covariate), function(nm) {
      df <- correlate_features_to_covariate(
        expr_mat  = expr_mat,
        covariate = covariate[[nm]],
        method    = method,
        features  = features,
        p_value   = p_value,
        mi_bins   = mi_bins
      )
      df$covariate <- nm
      df[, c("covariate", setdiff(names(df), "covariate"))]
    })
    return(do.call(rbind, results))
  }

  # -- Feature subsetting -------------------------------------------------------
  if (!is.null(features)) {
    keep <- intersect(features, rownames(expr_mat))
    if (length(keep) == 0) stop("None of the requested features found in expr_mat.")
    expr_mat <- expr_mat[keep, , drop = FALSE]
  }

  # -- Drop NA cells ------------------------------------------------------------
  ok <- !is.na(covariate)
  if (sum(ok) < 3) stop("Fewer than 3 non-NA cells — cannot compute associations.")
  expr_mat  <- expr_mat[, ok, drop = FALSE]
  covariate <- covariate[ok]

  # -- Drop zero-mean features --------------------------------------------------
  row_means <- Matrix::rowMeans(expr_mat)
  expr_mat  <- expr_mat[row_means != 0, , drop = FALSE]
  if (nrow(expr_mat) == 0) stop("All features have zero mean after filtering.")

  feat_names <- rownames(expr_mat)
  n_feat     <- nrow(expr_mat)

  # -- Dispatch -----------------------------------------------------------------

  # ── Classical correlation (fast vectorised path) ──────────────────────────
  if (method %in% c("pearson", "spearman", "kendall")) {
    if (!p_value) {
      vals <- stats::cor(t(as.matrix(expr_mat)), covariate, method = method)[, 1]
      result <- data.frame(feature = feat_names, association = unname(vals),
                           stringsAsFactors = FALSE)
    } else {
      result <- do.call(rbind, lapply(seq_len(n_feat), function(i) {
        ct <- stats::cor.test(as.numeric(expr_mat[i, ]), covariate, method = method)
        data.frame(feature = feat_names[i], association = unname(ct$estimate),
                   p_value = ct$p.value, stringsAsFactors = FALSE)
      }))
    }

  # ── Extended measures via {correlation} ──────────────────────────────────
  } else if (method %in% c("distance", "biweight", "hoeffding")) {
    if (!requireNamespace("correlation", quietly = TRUE))
      stop("Package 'correlation' is required for method '", method,
           "'. Install with: install.packages('correlation')")
    # correlation::cor_test() method names
    cor_method <- switch(method,
      distance  = "distance",
      biweight  = "biweight",
      hoeffding = "hoeffding"
    )
    cov_vec <- as.numeric(covariate)
    result <- do.call(rbind, lapply(seq_len(n_feat), function(i) {
      ct <- correlation::cor_test(as.numeric(expr_mat[i, ]), cov_vec,
                                  method = cor_method)
      row <- data.frame(feature = feat_names[i], association = ct$r,
                        stringsAsFactors = FALSE)
      if (p_value) row$p_value <- ct$p
      row
    }))

  # ── Mutual information ────────────────────────────────────────────────────
  } else if (method == "mutual_information") {
    if (!requireNamespace("infotheo", quietly = TRUE))
      stop("Package 'infotheo' is required for method 'mutual_information'. ",
           "Install with: install.packages('infotheo')")
    disc_cov <- infotheo::discretize(covariate, nbins = mi_bins)
    result <- do.call(rbind, lapply(seq_len(n_feat), function(i) {
      x_disc <- infotheo::discretize(as.numeric(expr_mat[i, ]), nbins = mi_bins)
      mi     <- infotheo::mutinformation(x_disc, disc_cov)
      row    <- data.frame(feature = feat_names[i], association = mi,
                           stringsAsFactors = FALSE)
      if (p_value) {
        # Chi-squared test on the joint frequency table as a proxy p-value
        ct <- suppressWarnings(
          stats::chisq.test(table(x_disc[[1]], disc_cov[[1]]))
        )
        row$p_value <- ct$p.value
      }
      row
    }))

  # ── Cosine similarity ─────────────────────────────────────────────────────
  } else if (method == "cosine") {
    cov_norm <- sqrt(sum(covariate^2))
    if (cov_norm == 0) stop("Covariate has zero norm — cosine similarity undefined.")
    mat_dense  <- as.matrix(expr_mat)                          # features × cells
    row_norms  <- sqrt(Matrix::rowSums(expr_mat^2))
    cosines    <- as.numeric(mat_dense %*% covariate) / (row_norms * cov_norm)
    cosines[row_norms == 0] <- NA_real_
    result <- data.frame(feature = feat_names, association = cosines,
                         stringsAsFactors = FALSE)
    if (p_value) {
      # Treat cosine as Pearson-equivalent r and derive t-stat
      n  <- ncol(expr_mat)
      tt <- cosines * sqrt(n - 2) / sqrt(1 - cosines^2)
      result$p_value <- 2 * stats::pt(-abs(tt), df = n - 2)
    }
  }

  # -- Sort and return ----------------------------------------------------------
  result <- result[order(abs(result$association), decreasing = TRUE, na.last = TRUE), ]
  rownames(result) <- NULL
  result
}

#' Run AddModuleScore on a Seurat object and return a cell × score data.frame.
#'
#' @param srt      Processed Seurat object.
#' @param genesets Named list of gene vectors.
#' @return data.frame with cell_id and one column per geneset (named after the list element).
extract_module_scores <- function(srt, genesets) {
  srt <- Seurat::AddModuleScore(srt, features = genesets)

  rename_map <- setNames(paste0("Cluster", seq_along(genesets)), names(genesets))
  srt@meta.data <- srt@meta.data %>% dplyr::rename(all_of(rename_map))

  score_cols <- names(genesets)
  srt@meta.data[, c("sample.name", score_cols), drop = FALSE] %>%
    tibble::rownames_to_column("cell_id")
}

#' Run UCell on a Seurat object and return a self-contained data.frame.
#'
#' Scores are rank-based (UCell) and do not depend on background cell set,
#' making them more robust than AddModuleScore for cross-sample comparisons.
#' UMAP coords and sample.name are included so downstream notebooks do not
#' need to load the full Seurat object.
#'
#' @param srt       Processed Seurat object.
#' @param genesets  Named list of gene vectors.
#' @param assay     Assay to score (default "SCT").
#' @param layer     Layer to use (default "data").
#' @param umap_name UMAP reduction name (default "umap.nn").
#' @return data.frame: cell_id, UMAP1, UMAP2, sample.name, one column per geneset.
extract_ucell_scores <- function(srt, genesets,
                                 assay     = "SCT",
                                 layer     = "data",
                                 umap_name = "umap.nn") {
  mat    <- Seurat::GetAssayData(srt, assay = assay, layer = layer)
  scores <- UCell::ScoreSignatures_UCell(mat, features = genesets)
  # UCell appends "_UCell" suffix — strip it
  colnames(scores) <- sub("_UCell$", "", colnames(scores))

  umap_emb           <- Seurat::Embeddings(srt, umap_name)[, 1:2]
  colnames(umap_emb) <- c("UMAP1", "UMAP2")

  data.frame(
    cell_id     = colnames(srt),
    UMAP1       = umap_emb[, 1],
    UMAP2       = umap_emb[, 2],
    sample.name = srt$sample.name,
    scores,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Extract UMAP coords + expression values for a list of genes.
#'
#' Thin wrapper for use in feature-expression targets in the DE notebook.
#' Returns a long-format data.frame: cell, UMAP1, UMAP2, sample.name, gene, expr.
#'
#' @param srt       Processed Seurat object.
#' @param genes     Character vector of gene names.
#' @param assay     Assay to pull expression from (default "SCT").
#' @param umap_name UMAP reduction name (default "umap.nn").
#' @return Long data.frame.
extract_feature_expr_long <- function(srt,
                                       genes,
                                       assay     = "SCT",
                                       umap_name = "umap.nn") {
  genes_present <- intersect(genes, rownames(srt))
  umap_emb      <- as.data.frame(Seurat::Embeddings(srt, umap_name)[, 1:2])
  colnames(umap_emb) <- c("UMAP1", "UMAP2")
  umap_emb$cell        <- rownames(umap_emb)
  umap_emb$sample.name <- srt$sample.name

  expr_mat <- as.matrix(Seurat::GetAssayData(srt, assay = assay, layer = "data")[genes_present, , drop = FALSE])
  expr_long <- as.data.frame(t(expr_mat))
  expr_long$cell <- rownames(expr_long)

  merged <- merge(umap_emb, expr_long, by = "cell")
  tidyr::pivot_longer(merged,
                      cols      = dplyr::all_of(genes_present),
                      names_to  = "gene",
                      values_to = "expr")
}

#map mouse to human genes
map_ortho_gene <- function(gene_old){
  gene_new <- stringr::str_to_upper(gene_old)
  names(gene_old) <- gene_new
  return(gene_old)
}

# ---------------------------------------------------------------
# SCENT potency score (expression entropy)
# ---------------------------------------------------------------
#convert gene Entrez ID to symbols
rename_scent_net <- function(){
  # Load SCENT network (human PPI, ENTREZ IDs)
  net <- SCENT::net17Jan16.m

  # Convert ENTREZ IDs to gene symbols (mouse)
  g_syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys      = rownames(net),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )

  # Filter valid symbols and rebuild network
  valid <- !is.na(g_syms)
  net_sym <- net[valid, valid]
  dimnames(net_sym) <- list(unname(g_syms[valid]), unname(g_syms[valid]))

  return(net_sym)
}
# Compute SCENT potency scores for mouse scRNA-seq data
# https://github.com/aet21/SCENT
get_scent_potency <- function(srt_obj, scent_net) {

  # Match network genes to Seurat object
  genes_in_srt <- rownames(srt_obj@assays$SCT@data)
  mouse_to_human <- map_ortho_gene(genes_in_srt)
  human_shared <- intersect(rownames(scent_net),names(mouse_to_human))
  mouse_shared <- mouse_to_human[human_shared]

  # Subset to shared genes
  net_filtered <- scent_net[human_shared, human_shared]
  rm(scent_net);gc()
  expr_mat     <- as.matrix(srt_obj@assays$SCT@data[mouse_shared,])
  rownames(expr_mat) <- human_shared

  # Run CCAT
  potency_scores <- SCENT::CompCCAT(expr_mat, net_filtered)

  # Return data.frame
  data.frame(
    cell.name      = colnames(expr_mat),
    potency        = potency_scores,
    scaled_potency = potency_scores / max(potency_scores),
    row.names      = colnames(expr_mat)
  )
}

# Compute CytoTRACE differentiation scores
# Higher scores = more stem-like, lower scores = more differentiated
# https://github.com/digitalcytometry/cytotrace
get_cytotrace <- function(srt_obj) {
  # Extract counts matrix (CytoTRACE works on counts, not normalized data)
  counts_mat <- as.matrix(srt_obj@assays$SCT@counts)[1:1000,100]

  # Run CytoTRACE
  cytotrace_result <- CytoTRACE2::cytotrace2(counts_mat,
                                             species = "mouse",
                                             ncores = 2)

  # Extract scores
  scores <- cytotrace_result$CytoTRACE

  # Return data.frame
  data.frame(
    cell.name = names(scores),
    cytotrace = scores,
    row.names = names(scores)
  )
}

# ---------------------------------------------------------------
# Alternative dimensionality reduction methods (PHATE, FA2)
# ---------------------------------------------------------------

# Run PHATE on PCA embeddings
# https://github.com/KrishnaswamyLab/phateR
run_phate <- function(srt_obj,
                      reduction = "pca.SCT",
                      dims      = 1:50,
                      knn       = 5,
                      decay     = 40,
                      t         = "auto") {
  set.seed(42)

  # Get PCA embeddings
  input_mat <- Seurat::Embeddings(srt_obj, reduction = reduction)[, dims]

  # Run PHATE
  phate_result <- phateR::phate(
    input_mat,
    ndim       = 2,
    knn        = knn,
    decay      = decay,
    t          = t,
    gamma      = 1,
    n.landmark = 2000
  )

  # Add to Seurat object
  colnames(phate_result$embedding) <- c("PHATE_1", "PHATE_2")
  rownames(phate_result$embedding) <- colnames(srt_obj)

  srt_obj[["phate"]] <- Seurat::CreateDimReducObject(
    embeddings = phate_result$embedding,
    key        = "PHATE_",
    assay      = "SCT"
  )

  srt_obj
}

# Run ForceAtlas2 (Fruchterman-Reingold approximation) on SNN graph
# Uses igraph::layout_with_fr (true FA2 requires python-fa2 or Gephi)
run_fa2 <- function(srt_obj, graph = NULL, niter = 500) {
  set.seed(42)

  # Auto-detect SNN graph
  if (is.null(graph)) {
    snn_graphs <- grep("sNN", names(srt_obj@graphs), value = TRUE)
    graph <- snn_graphs[1]
  }

  # Extract graph as igraph object
  graph_mat <- srt_obj@graphs[[graph]]
  ig <- igraph::graph_from_adjacency_matrix(graph_mat, mode = "undirected", weighted = TRUE)

  # Run Fruchterman-Reingold layout
  layout_result <- igraph::layout_with_fr(ig, niter = niter)

  # Add to Seurat object
  colnames(layout_result) <- c("FA2_1", "FA2_2")
  rownames(layout_result) <- colnames(srt_obj)

  srt_obj[["fa2"]] <- Seurat::CreateDimReducObject(
    embeddings = layout_result,
    key        = "FA2_",
    assay      = "SCT"
  )

  srt_obj
}

# Run PHATE and return the embedding matrix only (cells × 2).
# Used by tar_map targets so the full Seurat object is not cached.
compute_phate_coords <- function(srt, reduction = "pca.SCT", dims = 1:50,
                                 knn = 5, decay = 40, t = "auto") {
  srt <- run_phate(srt, reduction = reduction, dims = dims,
                   knn = knn, decay = decay, t = t)
  Seurat::Embeddings(srt, "phate")
}

# Run FA2 and return the embedding matrix only (cells × 2).
compute_fa2_coords <- function(srt, graph = NULL, niter = 500) {
  srt <- run_fa2(srt, graph = graph, niter = niter)
  Seurat::Embeddings(srt, "fa2")
}

# ── RNA velocity: h5ad export ────────────────────────────────────────────────

#' Export a processed Seurat object to AnnData h5ad for scVelo / CellRank.
#'
#' Converts via SingleCellExperiment (zellkonverter::writeH5AD).
#' Reductions are renamed to scVelo-expected names (X_umap, X_pca).
#' The SCT normalised data layer lands in adata.layers["logcounts"] (scanpy
#' convention); adata.X is set to the same.
#'
#' @param srt            Seurat object (DietSeurat-slimmed is fine).
#' @param path           Output file path (.h5ad).
#' @param assay          Assay to export (default "SCT").
#' @param reduction_umap UMAP reduction name in the Seurat object.
#' @param reduction_pca  PCA reduction name in the Seurat object.
#' @return `path` invisibly (for use as a targets file target).
export_seurat_h5ad <- function(srt,
                                path,
                                assay          = "SCT",
                                reduction_umap = "umap.nn",
                                reduction_pca  = paste0("pca.", assay)) {
  # Normalised expression matrix — cells × genes (dgCMatrix)
  X   <- Matrix::t(Seurat::GetAssayData(srt, assay = assay, layer = "data"))

  # Cell and gene metadata
  obs <- srt@meta.data
  var <- data.frame(row.names = rownames(srt[[assay]]))

  # UMAP and PCA embeddings under scVelo-expected keys
  obsm <- list()
  if (reduction_umap %in% Seurat::Reductions(srt))
    obsm[["X_umap"]] <- Seurat::Embeddings(srt, reduction_umap)
  if (reduction_pca %in% Seurat::Reductions(srt))
    obsm[["X_pca"]] <- Seurat::Embeddings(srt, reduction_pca)

  # SNN and KNN graphs → obsp so scanpy can use the pipeline graph directly
  # SNN: shared-nearest-neighbour weighted adjacency (used for clustering)
  # KNN: binary adjacency (used internally by Seurat FindNeighbors)
  obsp <- list()
  graph_names <- names(srt@graphs)
  snn_name <- grep("^sNN", graph_names, value = TRUE)[1]
  knn_name <- grep("^kNN|^RNA_nn|^SCT_nn|_nn$", graph_names, value = TRUE)[1]
  if (!is.na(snn_name))
    obsp[["connectivities"]] <- as(srt@graphs[[snn_name]], "CsparseMatrix")
  if (!is.na(knn_name))
    obsp[["distances"]] <- as(srt@graphs[[knn_name]], "CsparseMatrix")

  adata <- anndataR::AnnData(X = X, obs = obs, var = var, obsm = obsm,
                              obsp = if (length(obsp)) obsp else NULL)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) file.remove(path)
  adata$write_h5ad(path)
  invisible(path)
}

# stacked violin plot
#https://divingintogeneticsandgenomics.com/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    #scale_x_discrete(position = "top")+
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(2), angle = 0),
          axis.text.y = element_text(size = rel(2)),
          plot.margin = plot.margin )
  return(p)
}

modify_vlnplot_wide <- function(obj,
                                feature,
                                pt.size = 0,
                                plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                                ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ...) +
    coord_flip() +
    xlab("") + ylab(NULL) + ggtitle(feature) +
    theme(plot.title = element_text(size = 20,hjust = 0,face = "plain",angle = 45,  vjust=0.5,margin = margin(20,0,-30,0)),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size =20, angle = 0, hjust = 0),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(size = 20, angle = 45, vjust=1, hjust=1),
          axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackedVlnPlot_wide<- function(obj, features,
                               pt.size = 0,
                               plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                               ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot_wide(obj = obj,feature = x, ...))

  # Add back y-axis title to left-most plot. patchwork is going to support this?
  plot_list[[1]]<- plot_list[[1]] +
    theme(axis.text.y = element_text(size = 20, hjust=1),
          axis.ticks.y = element_blank())

  # change the x-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, nrow = 1)
  return(p)
}

#rename control clusters and plot stacked violion plot
plot_srt_ctr_stacked_vln <- function(srt_ctr, clust_id, feats, cols){
  gg <- StackedVlnPlot_wide(obj = srt_ctr,
                            group.by = clust_id,
                            features = feats
  )&
    scale_fill_manual(values = cols)

  #add x title
  p_lab <- ggplot() +
    annotate(geom = "text", x = 1, y = 1, label = "Normalized expression", angle = 0,size = 7) +
    coord_cartesian(clip = "off")+
    theme_void()

  gg <- (gg/p_lab)+
    plot_layout(heights = c( 1, 0.1))
  rm(srt_ctr)

  return(gg)
}


# ── Pseudobulk expression heatmap ─────────────────────────────────────────────

#' Pseudobulk expression heatmap via ComplexHeatmap.
#'
#' Averages expression (mean of log-normalised data layer) per group, z-scores
#' per gene, clusters rows and columns with Ward.D2 + dendsort, and draws a
#' ComplexHeatmap with an optional column annotation bar coloured by a metadata variable.
#'
#' @param srt              Seurat object.
#' @param genes            Character vector of genes to display.
#' @param group_by         Character vector of metadata column(s) passed to
#'                         AggregateExpression group.by.
#' @param annotate_by      Metadata column used to colour the column annotation
#'                         bar (e.g. "sample.name"). NULL to skip.
#' @param assay            Assay to aggregate (default "SCT").
#' @param layer            Layer to use (default "data").
#' @param scale            Z-score rows before clustering/plotting (default TRUE).
#' @param col_fun          circlize::colorRamp2 colour function. NULL = auto
#'                         (blue–white–red for scaled, white–red for raw).
#' @param title            Heatmap column title string.
#' @param annotation_colors Named list; names match annotate_by values, values
#'                         are named colour vectors. NULL = auto palette.
#' @param ...              Additional arguments forwarded to ComplexHeatmap::Heatmap.
#' @return A ComplexHeatmap Heatmap object.
plot_pseudobulk_heatmap <- function(srt,
                                    genes,
                                    group_by,
                                    annotate_by       = "sample.name",
                                    assay             = "SCT",
                                    layer             = "data",
                                    scale             = TRUE,
                                    col_fun           = NULL,
                                    title             = "",
                                    annotation_colors = NULL,
                                    ...) {
  # ── 1. Average expression (mean of log-normalised data per group) ──────────
  agg <- Seurat::AverageExpression(
    srt,
    assays        = assay,
    features      = genes,
    group.by      = group_by,
    layer         = layer,
    return.seurat = FALSE,
  )[[assay]]

  genes_present <- intersect(genes, rownames(agg))
  if (length(genes_present) == 0)
    stop("None of the requested genes found in the aggregated matrix.")
  mat <- as.matrix(agg[genes_present, , drop = FALSE])

  # ── 2. Row z-score ────────────────────────────────────────────────────────
  if (scale) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }

  # ── 3. Colour function ────────────────────────────────────────────────────
  if (is.null(col_fun)) {
    if (scale) {
      lim     <- max(abs(mat), na.rm = TRUE)
      col_fun <- circlize::colorRamp2(c(-lim, 0, lim),
                                      c("#2166ac", "white", "#d73027"))
    } else {
      col_fun <- circlize::colorRamp2(c(0, max(mat, na.rm = TRUE)),
                                      c("white", "#d73027"))
    }
  }

  # ── 4. Column annotation — stacked proportion bar plot ───────────────────
  col_anno <- NULL
  if (!is.null(annotate_by)) {
    meta <- srt@meta.data

    # Build per-cell group key matching AggregateExpression column names.
    # AggregateExpression prepends "g" to purely numeric group names.
    group_key <- if (length(group_by) == 1) {
      as.character(meta[[group_by]])
    } else {
      do.call(paste, c(lapply(group_by, function(g) as.character(meta[[g]])),
                       list(sep = "_")))
    }
    anno_vals   <- as.character(meta[[annotate_by]])
    anno_levels <- unique(anno_vals)

    # Colour palette for annotate_by levels
    if (!is.null(annotation_colors) && annotate_by %in% names(annotation_colors)) {
      pal <- annotation_colors[[annotate_by]]
      # Ensure all levels are covered
      missing_lvl <- setdiff(anno_levels, names(pal))
      if (length(missing_lvl))
        pal <- c(pal, setNames(scales::hue_pal()(length(missing_lvl)), missing_lvl))
    } else {
      pal <- setNames(scales::hue_pal()(length(anno_levels)), anno_levels)
    }

    # Proportion matrix: one row per heatmap column, one column per annotate_by level.
    # Each row sums to 1.
    prop_mat <- do.call(rbind, lapply(colnames(mat), function(col_name) {
      key <- col_name
      if (!key %in% group_key) key <- sub("^g", "", col_name)   # strip "g" prefix
      mask <- group_key == key
      tab  <- table(factor(anno_vals[mask], levels = anno_levels))
      as.numeric(tab) / max(sum(tab), 1L)
    }))
    rownames(prop_mat) <- colnames(mat)
    colnames(prop_mat) <- anno_levels

    bar_anno <- ComplexHeatmap::anno_barplot(
      prop_mat,
      ylim       = c(0, 1),
      bar_width  = 0.9,
      gp         = grid::gpar(fill = pal[anno_levels], col = NA),
      axis_param = list(side = "left", gp = grid::gpar(fontsize = 7))
    )

    col_anno <- do.call(
      ComplexHeatmap::HeatmapAnnotation,
      c(
        setNames(list(bar_anno), annotate_by),
        list(
          annotation_name_side = "left",
          annotation_name_gp   = grid::gpar(fontsize = 9)
        )
      )
    )
  }

  # ── 5. Dendrogram (Ward.D2 + dendsort) ────────────────────────────────────
  dend_row <- as.dendrogram(
    dendsort::dendsort(hclust(dist(mat),    method = "ward.D2"))
  )
  dend_col <- as.dendrogram(
    dendsort::dendsort(hclust(dist(t(mat)), method = "ward.D2"))
  )

  # ── 6. Draw heatmap ───────────────────────────────────────────────────────
  ComplexHeatmap::Heatmap(
    mat,
    col                  = col_fun,
    name                 = if (scale) "z-score" else "expr",
    cluster_rows         = dend_row,
    cluster_columns      = dend_col,
    show_row_dend        = TRUE,
    show_column_dend     = TRUE,
    row_dend_width       = unit(1.5, "cm"),
    column_dend_height   = unit(1.5, "cm"),
    top_annotation       = col_anno,
    row_names_gp         = grid::gpar(fontsize = 9),
    column_names_gp      = grid::gpar(fontsize = 9),
    column_names_rot     = 45,
    column_title         = title,
    column_title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp  = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8)
    ),
    ...
  )
}

# ── Neighbour / spatial-autocorrelation utilities ─────────────────────────────

#' Moran's I using a sparse weight matrix
#'
#' Computes Moran's I for a single numeric vector using a sparse spatial weight
#' matrix.  Avoids the n² dense allocation that ape::Moran.I and Rfast2::moranI
#' require, making it suitable for large single-cell graphs (~50k cells).
#'
#' @param x Numeric vector of expression values (one gene, one cell per entry).
#' @param W Sparse weight matrix (dgCMatrix, cells × cells). Should have zeros
#'          on the diagonal and be row-standardised (style "W") or at least
#'          non-negative. No p-value is computed.
#'
#' @return Scalar Moran's I statistic (range roughly −1 to 1; higher = more
#'         locally autocorrelated).
moran_sparse <- function(x, W) {
  n     <- length(x)
  x_c   <- x - mean(x)
  W_sum <- sum(W)   # base::sum dispatches on sparseMatrix via S4
  (n / W_sum) * as.numeric(Matrix::crossprod(x_c, W %*% x_c)) / sum(x_c^2)
}

# ── Neighbour utilities ────────────────────────────────────────────────────────

#' Build a symmetric distance matrix from a Seurat Neighbor object
#'
#' Converts the asymmetric KNN index+distance arrays stored in a Seurat
#' Neighbor slot into a symmetric sparse distance matrix.
#'
#' @param nn     A Seurat Neighbor object (e.g. srt@neighbors$RNA.nn).
#' @param method How to symmetrize: "average" (default) takes the mean of the
#'               two directed distances; "max" takes the larger (union of edges);
#'               "min" takes the smaller (intersection of edges).
#'
#' @return A symmetric dgCMatrix (n_cells × n_cells) of distances.
nn_to_sym_dist <- function(nn, method = c("average", "max", "min")) {
  method <- match.arg(method)

  idx  <- nn@nn.idx    # n_cells × k
  dist <- nn@nn.dist   # n_cells × k
  n    <- nrow(idx)
  k    <- ncol(idx)

  i    <- rep(seq_len(n), k)
  j    <- as.vector(idx)
  x    <- as.vector(dist)
  keep <- i != j
  D    <- Matrix::sparseMatrix(
    i        = i[keep],
    j        = j[keep],
    x        = x[keep],
    dims     = c(n, n),
    dimnames = list(nn@cell.names, nn@cell.names)
  )

  D_sym <- switch(method,
    average = (D + Matrix::t(D)) / 2,
    max     = Matrix::pmax(D, Matrix::t(D)),
    min     = Matrix::pmin(D, Matrix::t(D))
  )
  D_sym
}

#' Build protein-coding promoter and gene-body GRanges from an EnsDb
#'
#' Mirrors the promoter-building logic inside hdWGCNA::MotifScan.
#' Returns a named list:
#'   $promoters   — transcript-level promoter windows, UCSC chr names, symbol column
#'   $gene_coords — gene-body coordinates, UCSC chr names, symbol column
#'
#' Both are filtered to protein-coding genes on standard (1–100, X, Y) chromosomes.
.ensdb_add_ucsc_chr <- function(gr, species_genome) {
  gr <- GenomeInfoDb::keepSeqlevels(
    gr, levels(base::droplevels(GenomeInfoDb::seqnames(gr))))
  old <- levels(GenomeInfoDb::seqnames(gr))
  gr  <- GenomeInfoDb::renameSeqlevels(
    gr, stats::setNames(paste0("chr", old), old))
  GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gr)) <- species_genome
  gr
}

#' Protein-coding gene body GRanges from EnsDb (UCSC chr names, symbol column)
build_protein_coding_genes_gr <- function(EnsDb, species_genome = "mm10") {
  std_chr     <- c(as.character(1:100), "X", "Y")
  gene.coords <- ensembldb::genes(EnsDb)
  gene.coords <- gene.coords[gene.coords$gene_biotype == "protein_coding"]
  gene.coords <- gene.coords[
    as.character(GenomeInfoDb::seqnames(gene.coords)) %in% std_chr]
  .ensdb_add_ucsc_chr(gene.coords, species_genome)
}

#' Protein-coding promoter GRanges from EnsDb (UCSC chr names, symbol column)
#' Symbols are resolved via genes_gr (output of build_protein_coding_genes_gr).
build_protein_coding_promoters_gr <- function(EnsDb, genes_gr, species_genome = "mm10") {
  std_chr        <- c(as.character(1:100), "X", "Y")
  gene.promoters <- ensembldb::promoters(EnsDb)
  gene.promoters <- gene.promoters[gene.promoters$tx_biotype == "protein_coding"]
  gene.promoters <- gene.promoters[
    as.character(GenomeInfoDb::seqnames(gene.promoters)) %in% std_chr]
  gene.promoters$symbol <- genes_gr$symbol[
    base::match(gene.promoters$gene_id, names(genes_gr))]
  gene.promoters <- .ensdb_add_ucsc_chr(gene.promoters, species_genome)

  GenomicRanges::GRanges(
    seqnames = base::droplevels(GenomeInfoDb::seqnames(gene.promoters)),
    IRanges::IRanges(
      start = GenomicRanges::start(gene.promoters),
      end   = GenomicRanges::end(gene.promoters)
    ),
    symbol = gene.promoters$symbol,
    genome = species_genome
  )
}

#' Link genomic ranges to protein-coding gene bodies, then reduce.
#' Returns list(gr = reduced GRanges, peak_gene = data.frame(reduced_name, gene)).
link_ranges_to_genes <- function(gr, genes_gr, link_kb = 100) {
  rng_exp <- gr
  GenomicRanges::start(rng_exp) <- pmax(
    1L, GenomicRanges::start(rng_exp) - link_kb * 1000L)
  GenomicRanges::end(rng_exp) <- GenomicRanges::end(rng_exp) + link_kb * 1000L

  hits    <- GenomicRanges::findOverlaps(rng_exp, genes_gr)
  pg_orig <- data.frame(
    rng_idx = S4Vectors::queryHits(hits),
    gene    = genes_gr$symbol[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  pg_orig <- pg_orig[!is.na(pg_orig$gene), ]

  linked_idx <- unique(pg_orig$rng_idx)
  linked_gr  <- gr[linked_idx]
  red_gr     <- plyranges::reduce_ranges(linked_gr)
  names(red_gr) <- paste0("range:", seq_along(red_gr))

  ovlp        <- GenomicRanges::findOverlaps(linked_gr, red_gr)
  orig_to_red <- data.frame(
    rng_idx      = linked_idx[S4Vectors::queryHits(ovlp)],
    reduced_name = names(red_gr)[S4Vectors::subjectHits(ovlp)],
    stringsAsFactors = FALSE
  )
  peak_gene <- dplyr::distinct(
    dplyr::left_join(pg_orig, orig_to_red, by = "rng_idx",
                     relationship = "many-to-many"),
    reduced_name, gene
  )
  peak_gene <- peak_gene[!is.na(peak_gene$reduced_name), ]

  cat(sprintf("[link_ranges_to_genes] %d input → %d gene-linked → %d reduced; %d range–gene links\n",
              length(gr), length(linked_idx), length(red_gr), nrow(peak_gene)))
  list(gr = red_gr, peak_gene = peak_gene)
}

# Build peak-gene links from the TSS and genes.within.100Kb annotation columns
# already present in the Yoshida 2020 ImmGenATAC18 OCR table.
# Returns list(gr, peak_gene) matching the link_ranges_to_genes output format.
parse_yoshida_orig_links <- function(ocr_gr, ocr_df) {
  peak_ids_in_gr <- sub("^atac:", "", names(ocr_gr))
  df_sub         <- ocr_df[match(peak_ids_in_gr, ocr_df$peakID), ]

  # TSS column — single gene per peak
  tss_mask <- !is.na(df_sub$TSS) & nchar(trimws(df_sub$TSS)) > 0
  tss_df <- data.frame(
    reduced_name = names(ocr_gr)[tss_mask],
    gene         = trimws(df_sub$TSS[tss_mask]),
    stringsAsFactors = FALSE
  )

  # genes.within.100Kb column — comma-separated, may be multiple per peak
  kb100_str  <- df_sub[["genes.within.100Kb"]]
  kb100_mask <- !is.na(kb100_str) & nchar(trimws(kb100_str)) > 0
  split_genes <- strsplit(kb100_str[kb100_mask], ",")
  kb100_df <- data.frame(
    reduced_name = rep(names(ocr_gr)[kb100_mask], lengths(split_genes)),
    gene         = trimws(unlist(split_genes)),
    stringsAsFactors = FALSE
  )
  kb100_df <- kb100_df[nchar(kb100_df$gene) > 0, ]

  peak_gene    <- dplyr::distinct(rbind(tss_df, kb100_df))
  linked_names <- unique(peak_gene$reduced_name)
  gr           <- ocr_gr[names(ocr_gr) %in% linked_names]

  cat(sprintf(
    "[parse_yoshida_orig_links] %d peaks → %d gene-linked (orig TSS/100kb cols); %d peak-gene links\n",
    length(ocr_gr), length(gr), nrow(peak_gene)
  ))
  list(gr = gr, peak_gene = peak_gene)
}

# ---------------------------------------------------------------
# Trajectory — binned expression along pseudotime
# ---------------------------------------------------------------

#' Bin gene expression along a single Slingshot pseudotime lineage.
#'
#' Cells with NA pseudotime (not assigned to this lineage) are excluded.
#' Bins with no cells return NA columns, which are then linearly interpolated.
#' Optional per-gene LOESS smoothing is applied across bins.
#'
#' @param srt          Seurat object containing the expression data.
#' @param sling_df     Data frame returned by \code{sling_cell_df()};
#'                     rows = cells, columns include pseudotime per lineage.
#' @param lineage_col  Name of the raw pseudotime column (e.g. "Lineage1").
#' @param n_bins       Number of equal-width pseudotime bins (default 50).
#' @param assay        Seurat assay to pull expression from (default "SCT").
#' @param layer        Layer within that assay (default "data").
#' @param genes        Genes to include.  NULL → HVGs from the Seurat object.
#' @param smooth       Apply per-gene LOESS smoothing across bins (default TRUE).
#' @param smooth_span  LOESS span (default 0.3).
#' @return A \code{genes × n_bins} numeric matrix.
#'   Column names are the midpoint of each pseudotime bin (0–100 scale).
bin_expr_along_pseudotime <- function(srt, sling_df, lineage_col,
                                       n_bins      = 50L,
                                       assay       = "SCT",
                                       layer       = "data",
                                       genes       = NULL,
                                       smooth      = TRUE,
                                       smooth_span = 0.3) {
  # Cells assigned to this lineage
  pt_raw    <- sling_df[[lineage_col]]
  use_cells <- rownames(sling_df)[!is.na(pt_raw)]
  pt_use    <- pt_raw[!is.na(pt_raw)]

  if (length(use_cells) < 2L)
    stop(sprintf("Lineage '%s' has fewer than 2 cells with non-NA pseudotime.", lineage_col))

  if (is.null(genes)) genes <- Seurat::VariableFeatures(srt)
  genes <- intersect(genes, rownames(srt))
  if (length(genes) == 0L) stop("No matching genes found in the Seurat object.")

  # Dense sub-matrix (genes × lineage-cells)
  expr <- as.matrix(
    Seurat::GetAssayData(srt, assay = assay, layer = layer)[genes, use_cells, drop = FALSE])

  # Normalise pseudotime to 0–100 for interpretable column labels
  rng    <- range(pt_use)
  pt_norm <- (pt_use - rng[1]) / (rng[2] - rng[1]) * 100

  breaks  <- seq(0, 100, length.out = n_bins + 1L)
  bin_ids <- cut(pt_norm, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  bin_mids <- round((breaks[-length(breaks)] + breaks[-1L]) / 2, 1)

  # Mean expression per bin
  mat <- do.call(cbind, lapply(seq_len(n_bins), function(b) {
    idx <- which(bin_ids == b)
    if (length(idx) == 0L) return(rep(NA_real_, nrow(expr)))
    if (length(idx) == 1L) return(as.numeric(expr[, idx]))
    rowMeans(expr[, idx, drop = FALSE])
  }))
  rownames(mat) <- genes
  colnames(mat) <- bin_mids

  # Linearly interpolate empty bins
  mat <- t(apply(mat, 1L, function(x) {
    if (sum(!is.na(x)) < 2L) return(x)
    stats::approx(seq_along(x), x, xout = seq_along(x), rule = 2L)$y
  }))
  rownames(mat) <- genes

  # Optional per-gene LOESS smoothing
  if (smooth) {
    bin_seq <- seq_len(n_bins)
    mat <- t(apply(mat, 1L, function(x) {
      tryCatch(
        stats::predict(stats::loess(x ~ bin_seq, span = smooth_span, degree = 1L), bin_seq),
        error = function(e) x
      )
    }))
    rownames(mat) <- genes
  }

  mat
}


#' Plot binned expression heatmaps for one or more Slingshot lineages side by side.
#'
#' Lineages share the same row (gene) order, enabling direct visual comparison
#' of expression dynamics across trajectories.  Requires \pkg{ComplexHeatmap}
#' and \pkg{circlize}.
#'
#' @param bin_mats   Named list of \code{genes × n_bins} matrices as returned
#'                   by \code{bin_expr_along_pseudotime()}.
#' @param lineages   Character vector: which list elements to display.
#'                   NULL → all.
#' @param genes      Restrict to these genes (rows).  NULL → all shared rows.
#' @param n_top      Keep the \code{n_top} genes with the highest bin-variance
#'                   (computed on the first selected lineage).  NULL → keep all.
#' @param gene_order How to order rows:
#'                   \code{"peak"} – by pseudotime bin of maximum expression
#'                     (first lineage); shows expression waves;
#'                   \code{"clust"} – Ward.D2 hierarchical clustering on the
#'                     first lineage;
#'                   character vector – explicit gene order.
#' @param scale      Z-score each gene across bins before plotting (default TRUE).
#' @param col_fun    circlize colour function.  NULL → auto (RdBu for scaled,
#'                   viridis for raw).
#' @param show_row_names  Logical.  NULL → auto (TRUE when ≤ 60 genes).
#' @param row_anno   Optional \code{ComplexHeatmap::rowAnnotation()} added to
#'                   the left of the first panel.
#' @param ...        Additional arguments forwarded to each \code{Heatmap()} call.
#' @return A \code{HeatmapList} object (use \code{ComplexHeatmap::draw()} or
#'         just print to display).
plot_lineage_heatmap <- function(bin_mats,
                                  lineages        = NULL,
                                  genes           = NULL,
                                  n_top           = NULL,
                                  gene_order      = "peak",
                                  scale           = TRUE,
                                  col_fun         = NULL,
                                  show_row_names  = NULL,
                                  row_names_gp    = NULL,
                                  row_anno        = NULL,
                                  ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Please install ComplexHeatmap.")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Please install circlize.")

  if (is.null(lineages)) lineages <- names(bin_mats)
  mats <- bin_mats[lineages]

  # Restrict to genes present in all selected lineages
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  if (!is.null(genes)) common_genes <- intersect(genes, common_genes)
  if (length(common_genes) == 0L) stop("No genes remain after intersecting lineages.")
  mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])

  # Top-N by bin variance (first lineage)
  if (!is.null(n_top) && n_top < length(common_genes)) {
    gvar     <- apply(mats[[1L]], 1L, stats::var, na.rm = TRUE)
    top_keep <- names(sort(gvar, decreasing = TRUE))[seq_len(n_top)]
    mats         <- lapply(mats, function(m) m[top_keep, , drop = FALSE])
    common_genes <- top_keep
  }

  # Z-score each gene across bins (within each lineage independently)
  if (scale) {
    mats <- lapply(mats, function(m) {
      m_scaled <- t(scale(t(m)))
      # Replace NaN (zero-variance genes) with 0
      m_scaled[!is.finite(m_scaled)] <- 0
      m_scaled
    })
  }

  # Gene ordering
  ord <- switch(
    if (is.character(gene_order) && length(gene_order) == 1L) gene_order else "explicit",
    peak = order(apply(mats[[1L]], 1L, which.max)),
    clust = {
      d <- stats::dist(mats[[1L]])
      stats::hclust(d, method = "ward.D2")$order
    },
    "explicit" = {
      if (length(gene_order) > 1L) match(gene_order, common_genes) else seq_along(common_genes)
    }
  )
  ord  <- ord[!is.na(ord)]
  mats <- lapply(mats, function(m) m[ord, , drop = FALSE])

  # Auto colour function
  if (is.null(col_fun)) {
    all_vals <- unlist(lapply(mats, function(m) m[is.finite(m)]))
    if (scale) {
      lim     <- min(stats::quantile(abs(all_vals), 0.99, na.rm = TRUE), 3)
      col_fun <- circlize::colorRamp2(c(-lim, 0, lim),
                                       c("#2166ac", "white", "#b2182b"))
    } else {
      hi      <- stats::quantile(all_vals, 0.99, na.rm = TRUE)
      col_fun <- circlize::colorRamp2(c(0, hi / 2, hi),
                                       c("#440154", "#21908c", "#fde725"))
    }
  }

  n_genes  <- nrow(mats[[1L]])
  n_bins   <- ncol(mats[[1L]])
  show_rn  <- if (!is.null(show_row_names)) show_row_names else (n_genes <= 60L)
  rn_gp    <- if (!is.null(row_names_gp)) row_names_gp else
                grid::gpar(fontsize = max(4, min(9, floor(200 / n_genes))))

  # Column ticks: label every 25 % of pseudotime
  tick_pos   <- round(seq(1, n_bins, length.out = 5))
  col_labels <- rep("", n_bins)
  col_labels[tick_pos] <- paste0(seq(0, 100, by = 25), "%")

  # Build HeatmapList
  ht_list <- NULL
  for (i in seq_along(mats)) {
    lin_nm   <- names(mats)[i]
    is_first <- i == 1L
    is_last  <- i == length(mats)

    ht <- ComplexHeatmap::Heatmap(
      mats[[i]],
      name                = if (is_first) (if (scale) "Z-score" else "Expr") else paste0(" ", i),
      col                 = col_fun,
      cluster_rows        = FALSE,
      cluster_columns     = FALSE,
      show_row_names      = show_rn && is_first,
      row_names_gp        = rn_gp,
      column_title        = lin_nm,
      column_title_gp     = grid::gpar(fontsize = 10, fontface = "bold"),
      column_labels       = col_labels,
      column_names_gp     = grid::gpar(fontsize = 7),
      column_names_rot    = 0,
      show_heatmap_legend = is_last,
      use_raster          = TRUE,
      raster_quality      = 3,
      left_annotation     = if (is_first && !is.null(row_anno)) row_anno else NULL,
      ...
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }

  ht_list
}

# ── Colour scale helpers ──────────────────────────────────────────────────────

expr_scale_pal_gradientn <- function(pal, floor = 0.03) {
  ggplot2::scale_colour_gradientn(
    colours  = pal,
    values   = c(0, seq(floor, 1, length.out = length(pal) - 1)),
    na.value = pal[1]
  )
}

# ── MAGIC imputation ──────────────────────────────────────────────────────────

#' Run MAGIC on SCT normalised data from a Seurat object.
#'
#' @param srt        Seurat object with an SCT assay and `data` layer.
#' @param genes      Character vector of genes, or "all_genes" / "pca_only".
#' @param t          MAGIC diffusion time (default 3).
#' @param knn        k-NN for graph construction (default 5).
#' @param python_path Optional path to a Python executable for reticulate.
#' @return Long data.frame with columns: cell_id, gene, raw, magic, plus all
#'         metadata columns from srt@meta.data, and UMAP coordinates if present.
compute_magic <- function(srt,
                          genes        = "all_genes",
                          t            = 3,
                          knn          = 5,
                          python_path  = NULL) {
  if (!requireNamespace("Rmagic", quietly = TRUE))
    stop("Rmagic is not installed. Run: remotes::install_github('KrishnaswamyLab/MAGIC', subdir='Rmagic')")

  if (!is.null(python_path))
    reticulate::use_python(python_path, required = TRUE)

  # cells × genes normalised matrix (dense; MAGIC input)
  mat_raw <- Matrix::t(
    Seurat::GetAssayData(srt, assay = "SCT", layer = "data")
  )

  gene_set <- if (identical(genes, "all_genes")) rownames(srt) else genes
  gene_set <- intersect(gene_set, colnames(mat_raw))

  mat_sub <- as.matrix(mat_raw[, gene_set, drop = FALSE])

  magic_out <- Rmagic::magic(mat_sub, genes = gene_set, t = t, knn = knn)

  list(
    raw   = mat_sub,
    magic = magic_out$result,
    cells = rownames(mat_sub)
  )
}


