#!/usr/bin/env Rscript
# Driver for the integration-site analysis. The work is in src/R/integration-sites.R; this file only
# decides which tubes to read and where the tables go, in two stages with bwa in between.
#
# Usage, from the project root (scripts/slurm/analysis/integration-sites.sh does this for you):
#   env/bin/Rscript scripts/R/call-integration-sites.R extract <out_dir> <bam_dir> <tube>...
#   env/bin/Rscript scripts/R/call-integration-sites.R call    <out_dir> <genes.gtf.gz>
#
# extract writes work/discovery.fa and work/support.fa for bwa, plus the two tables the report plots
# (anchor_ends.csv, clip_waterfall.csv). call reads work/discovery.sam and work/support.sam back and
# writes integration_sites.csv, cell_site_assignments.csv and site_calling_summary.csv.

source("src/R/integration-sites.R")   # functions only; attaches nothing

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
out   <- args[2]
stopifnot(stage %in% c("extract", "call"), !is.na(out))
work <- file.path(out, "work")
dir.create(work, recursive = TRUE, showWarnings = FALSE)

if (stage == "extract") {
  bam_dir <- args[3]
  tubes   <- args[-(1:3)]
  stopifnot(length(tubes) > 0L)

  classified <- lapply(tubes, function(tube) {
    bam <- file.path(bam_dir, paste0("cellranger_", tube), "sample_alignments.bam")
    message(">>> ", tube)
    clips <- read_junction_clips(bam, tube)
    message("    ", nrow(clips), " soft-clipped reads in the scan window")
    classify_clips(clips)
  })
  classified <- dplyr::bind_rows(classified)

  readr::write_csv(clip_waterfall(classified),       file.path(out, "clip_waterfall.csv"))
  readr::write_csv(anchor_end_histogram(classified), file.path(out, "anchor_ends.csv"))
  writeLines(tubes, file.path(work, "tubes.txt"))

  n_disc <- write_clip_fasta(classified, "discovery", file.path(work, "discovery.fa"))
  n_supp <- write_clip_fasta(classified, "support",   file.path(work, "support.fa"))
  message("=== discovery clips: ", n_disc, "   support clips: ", n_supp)
  print(as.data.frame(clip_waterfall(classified)), row.names = FALSE)

} else {
  gtf   <- args[3]
  tubes <- readLines(file.path(work, "tubes.txt"))
  stopifnot(file.exists(gtf))

  discovery <- read_bwa_hits(file.path(work, "discovery.sam"))
  support   <- read_bwa_hits(file.path(work, "support.sam"))
  message(">>> ", nrow(discovery), " discovery and ", nrow(support),
          " support clips placed uniquely (MAPQ >= ", IS_DEFAULTS$min_mapq, ")")

  sites <- call_sites(discovery)
  hits  <- label_sites(attach_support(sites, support))

  site_table <- summarise_sites(hits, gtf, tubes)
  cell_table <- summarise_cells(hits)

  readr::write_csv(site_table, file.path(out, "integration_sites.csv"))
  readr::write_csv(cell_table, file.path(out, "cell_site_assignments.csv"))
  readr::write_csv(
    dplyr::tibble(
      stage = c("discovery clips placed uniquely",
                "discovery clips in a site that passed the cell and UMI floor",
                "support clips placed uniquely",
                "support clips joined to a discovered site",
                "sites called",
                "cell-site assignments",
                "cells at more than one site"),
      value = c(nrow(discovery), nrow(sites), nrow(support),
                sum(hits$tier == "support"), nrow(site_table), nrow(cell_table),
                dplyr::n_distinct(cell_table$barcode[cell_table$ambiguous]))
    ),
    file.path(out, "site_calling_summary.csv")
  )

  message("=== ", nrow(site_table), " sites, ", nrow(cell_table), " cell-site assignments")
  if (nrow(site_table)) {
    print(as.data.frame(site_table[, c("site_id", "chrom", "pos", "cells", "gene", "context")]),
          row.names = FALSE)
  }
}
