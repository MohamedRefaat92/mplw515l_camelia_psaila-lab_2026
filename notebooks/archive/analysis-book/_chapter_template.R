# _chapter_template.R
# Sourced by each analysis chapter to render standard sections.
# Expects these variables to exist in the calling environment:
#   cfg        — list from yaml::read_yaml("metadata/project_config.yaml")
#   output_dir — path to output/<run_label>/[subset/]
#   subset_label — human-readable name, e.g. "Full cohort" or "Stroma"

library(ggplot2)

chapter_figure <- function(...) file.path(output_dir, "figure", ...)
chapter_table  <- function(...) file.path(output_dir, "table",  ...)

marker_res <- cfg$processing$marker_resolution
