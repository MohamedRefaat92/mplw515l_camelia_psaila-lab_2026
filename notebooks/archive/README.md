# notebooks/archive/

Reports that can't render as they are: they read targets the pipeline no longer defines (renamed or
removed during the analysis). They are not registered in `_targets.R` and Quarto ignores this folder
(`_quarto.yml`). To bring one back: point its `tar_read()` / `tar_load()` calls at the current
targets, move it to `notebooks/<name>/<name>.qmd` with the provenance and session-info includes
(like the other reports), and add a `quarto_report()` line to `_targets.R`.

| Notebook | Reads, no longer in the pipeline | Current targets closest to it |
| --- | --- | --- |
| `compositional-analysis.qmd` | `milo_mnc`, `milo_stroma` | `milo_mnc_obj` + `milo_mnc_da`, `milo_stroma_obj` + `milo_stroma_da` |
| `differential-expression-hsc.qmd` | `cluster_markers_res_lst_hsc`, `delegate_pairwise_markers_hsc` | `cluster_de_each_vs_rest_hsc_*`, `condition_de_per_cluster_hsc_*` |
| `differential-expression-mnc.qmd` | `cluster_markers_res_lst_mnc`, `delegate_pairwise_markers_mnc` | `cluster_de_each_vs_rest_mnc_0_9`, `condition_de_per_cluster_mnc_0_9` |
| `differential-expression-stroma.qmd` | `cluster_markers_res_lst_stroma`, `delegate_pairwise_markers_stroma` | `cluster_de_each_vs_rest_stroma_0_9`, `condition_de_per_cluster_stroma_0_9` |
| `degs_gsea.qmd` | the above for HSC, plus `module_scores_mnc`, `module_scores_stroma` | as above; MNC/Stroma ORA needs the retired `scatter_df_*` targets (see `_targets.R`) |
| `fibrosis-myelofibrosis.qmd` | `delegate_allvsall_markers_hsc` | `condition_de_per_cluster_hsc_*` |
| `stress-response.qmd` | `delegate_allvsall_markers_hsc` | `condition_de_per_cluster_hsc_*` |
| `pca-gsea.qmd` | `mpn_processed_hsc_mnc` (the combined HSC + MNC object) | none: the HSC + MNC subset was removed |
| `differential-expression.qmd`, `index.qmd` | (the old book's index pages) | |
| `analysis-book/` | read CSV / PNG files of per-run output folders (`run_label`) that the pipeline no longer writes; `first-analysis-milestone-slides.qmd` reads `delegate_samples_markers_*` | `first-analysis-milestone-report` replaces them |

The "closest" column is a starting point: the replacements return differently shaped results.
