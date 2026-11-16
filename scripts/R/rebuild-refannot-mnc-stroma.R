# Rebuild only the reference-annotation targets of MNC and Stroma (and what they need).
# From the project root: env/bin/Rscript scripts/R/rebuild-refannot-mnc-stroma.R
targets::tar_make(
  names = c(
    # shared / upstream
    "ferchen_26_mm_bm_umap", "hematlas_mm_bm", "hematlas_mm_fl",
    "mca_mm_bm_ckit_umap", "mca_mm_bm_mesench_umap", "mca_mm_pbmc_umap",
    "zeng_25_hs_bm_umap", "mm_hs_ortho_long",
    # MNC
    "anno_meta_mnc", "annotation_recluster_mnc", "anno_expr_mnc",
    "proj_ferchen_26_mm_bm_mnc", "proj_hematlas_mm_bm_mnc", "proj_hematlas_mm_fl_mnc",
    "proj_mca_mm_bm_ckit_mnc", "proj_mca_mm_bm_mesench_mnc", "proj_mca_mm_pbmc_mnc",
    "proj_zeng_25_hs_bm_mnc",
    # Stroma
    "anno_meta_stroma", "annotation_recluster_stroma", "anno_expr_stroma",
    "proj_ferchen_26_mm_bm_stroma", "proj_hematlas_mm_bm_stroma", "proj_hematlas_mm_fl_stroma",
    "proj_mca_mm_bm_ckit_stroma", "proj_mca_mm_bm_mesench_stroma", "proj_mca_mm_pbmc_stroma",
    "proj_zeng_25_hs_bm_stroma"
  )
)
