#!/usr/bin/env Rscript
# Driver for the integration-site analysis. The work is in src/R/integration-sites.R; this file only
# decides which tubes to read and where the tables go.
#
# Three stages, because bwa runs twice and the second run needs the first's output:
#
#   extract <out> <bam_dir> <tube>...            clips out of the BAMs -> discovery.fa, assign.fa,
#                                                five_prime.fa, and the tables the report plots
#     bwa  -> discovery.sam, five_prime.sam      against the genome
#   sites   <out> <genes.gtf.gz> <genome.fa>     call the sites, write the mini-reference
#     bwa  -> assign.sam                         against the mini-reference
#   assign  <out>                                cell assignments, TSD check, summary
#
# Splitting sites from assign is not bookkeeping: the mini-reference is built from the called sites,
# so the clips can only be assigned once those exist. Discovery searches the genome, assignment
# searches the flanks - see read_mini_hits().

source("src/R/integration-sites.R")   # functions only; attaches nothing

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
out   <- args[2]
stopifnot(stage %in% c("extract", "sites", "assign"), !is.na(out))
work <- file.path(out, "work")
dir.create(work, recursive = TRUE, showWarnings = FALSE)

if (stage == "extract") {
  bam_dir <- args[3]
  tubes   <- args[-(1:3)]
  stopifnot(length(tubes) > 0L)

  classified <- dplyr::bind_rows(lapply(tubes, function(tube) {
    bam <- file.path(bam_dir, paste0("cellranger_", tube), "sample_alignments.bam")
    message(">>> ", tube)
    clips <- read_junction_clips(bam, tube)
    message("    ", nrow(clips), " soft-clipped reads in the scan window")
    classify_clips(clips)
  }))

  readr::write_csv(clip_waterfall(classified),       file.path(out, "clip_waterfall.csv"))
  readr::write_csv(anchor_end_histogram(classified), file.path(out, "anchor_ends.csv"))
  writeLines(tubes, file.path(work, "tubes.txt"))

  # discovery.fa goes to the genome and defines the sites; assign.fa is every placeable clip,
  # discovery ones included, and goes to the mini-reference built from those sites.
  n_disc <- write_clip_fasta(classified, "discovery", file.path(work, "discovery.fa"))
  n_asn  <- write_clip_fasta(classified, c("discovery", "assign"), file.path(work, "assign.fa"))
  message("=== discovery clips: ", n_disc, "   assignable clips: ", n_asn)

  # Second channel: the provirus's 5' junction, read off the same BAMs.
  five <- dplyr::bind_rows(lapply(tubes, function(tube) {
    bam <- file.path(bam_dir, paste0("cellranger_", tube), "sample_alignments.bam")
    classify_five_prime_clips(read_five_prime_clips(bam, tube))
  }))
  readr::write_csv(clip_waterfall(five), file.path(out, "clip_waterfall_5p.csv"))
  message("=== 5-prime clips: ",
          write_clip_fasta(five, "discovery", file.path(work, "five_prime.fa")))
  print(as.data.frame(clip_waterfall(classified)), row.names = FALSE)

} else if (stage == "sites") {
  gtf       <- args[3]
  genome_fa <- args[4]
  tubes     <- readLines(file.path(work, "tubes.txt"))
  stopifnot(file.exists(gtf), file.exists(genome_fa))

  discovery <- read_bwa_hits(file.path(work, "discovery.sam"), channel = "three_prime")
  message(">>> ", nrow(discovery), " discovery clips placed uniquely (MAPQ >= ",
          IS_DEFAULTS$min_mapq, ")")

  sites      <- label_sites(dplyr::mutate(call_sites(discovery), tier = "discovery"))
  site_table <- summarise_sites(sites, gtf, tubes)
  readr::write_csv(site_table, file.path(out, "integration_sites.csv"))
  saveRDS(list(n_discovery = nrow(discovery), n_in_sites = nrow(sites)),
          file.path(work, "site_counts.rds"))

  write_site_flanks(site_table, genome_fa, file.path(work, "site_flanks.fa"))
  message("=== ", nrow(site_table), " sites; mini-reference written")
  if (nrow(site_table)) {
    print(as.data.frame(site_table[, c("site_id", "chrom", "pos", "cells", "gene", "context")]),
          row.names = FALSE)
  }

} else {
  site_table <- readr::read_csv(file.path(out, "integration_sites.csv"), show_col_types = FALSE)
  counts     <- readRDS(file.path(work, "site_counts.rds"))

  mini <- read_mini_hits(file.path(work, "assign.sam"))
  mini <- dplyr::semi_join(mini, dplyr::select(site_table, site_id), by = "site_id")
  cell_table <- summarise_cells(dplyr::mutate(mini, tier = "assign"))
  readr::write_csv(cell_table, file.path(out, "cell_site_assignments.csv"))

  # The counts in the site table must come from the same place the cell table does.
  site_table <- update_site_counts(site_table, mini, readLines(file.path(work, "tubes.txt")))
  readr::write_csv(site_table, file.path(out, "integration_sites.csv"))

  five <- read_bwa_hits(file.path(work, "five_prime.sam"), channel = "five_prime")
  tsd  <- tsd_check(site_table, five)
  readr::write_csv(tsd, file.path(out, "tsd_check.csv"))

  # Per-base pile-up either side of each junction, for the report's browser-style figure.
  readr::write_csv(
    junction_coverage(file.path(work, "assign.sam"), file.path(work, "five_prime.sam"), site_table),
    file.path(out, "junction_coverage.csv"))

  readr::write_csv(
    dplyr::tibble(
      stage = c("discovery clips placed uniquely",
                "discovery clips in a site that passed the cell and UMI floor",
                "sites called",
                "clips assigned against the mini-reference",
                "cell-site assignments",
                "distinct cells labelled",
                "cells at more than one site"),
      value = c(counts$n_discovery, counts$n_in_sites, nrow(site_table), nrow(mini),
                nrow(cell_table),
                dplyr::n_distinct(paste(cell_table$tube, cell_table$barcode)),
                dplyr::n_distinct(cell_table$barcode[cell_table$ambiguous]))
    ),
    file.path(out, "site_calling_summary.csv")
  )

  if (nrow(tsd)) {
    message("=== target-site duplication (expect delta = 4)")
    print(as.data.frame(tsd), row.names = FALSE)
  }
  message("=== ", nrow(cell_table), " cell-site assignments over ",
          dplyr::n_distinct(paste(cell_table$tube, cell_table$barcode)), " cells")
}
