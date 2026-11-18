#!/usr/bin/env Rscript
# Export the hashtag demultiplex calls from the targets store to a flat file.
#
#   Rscript scripts/R/export-hashtag-calls.R <out_dir>
#
# Why this exists. Replicate attribution is load-bearing in three places - the cross-mouse
# co-occurrence control, the mono-replicate assignment check, and every B0304 caveat - and it lived
# only inside the targets store, which meant it was computed at RENDER time from a 2 GB Seurat
# dependency. A number that only exists while a report is rendering cannot be checked, cannot be
# diffed between runs, and silently changes if an upstream target is rebuilt.
#
# This makes it a file: one row per cell, the demultiplexer's own call, written once per analysis run
# and read by the report like any other input. The report keeps working if the targets store is
# rebuilt or absent.
#
# It is deliberately a separate script rather than a pipeline stage: the integration-site pipeline
# reads BAMs and knows nothing about targets, and giving it a dependency on the targets store to
# fetch a metadata column would be the wrong coupling.
suppressMessages({
  library(targets)
  library(dplyr)
})
source("src/R/integration-sites.R")   # for write_csv_checked()

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1L) args[1] else "output/objects/integration-sites"
stopifnot(dir.exists(out))

grab <- function(target, panel) {
  d <- targets::tar_read_raw(target, store = "output/_targets")
  dplyr::bind_rows(lapply(names(d), function(tube) {
    fa <- d[[tube]]$final_assign
    dplyr::tibble(panel = panel, tube = tube, cell = names(fa), call = as.character(fa))
  }))
}

calls <- dplyr::bind_rows(
  grab("demux_res_hsc",    "HSC"),
  grab("demux_res_mnc",    "MNC"),
  grab("demux_res_stroma", "stroma")
) |>
  # "multiplet" and "negative" are failures of the CALL, not animals. Kept, because the
  # cross-mouse control needs them: a droplet holding two mice's cells is exactly what a doublet
  # looks like, and dropping those rows would bias that control along the axis it tests.
  dplyr::mutate(replicate = sub("^TotalSeq-", "", call),
                is_singlet = grepl("^TotalSeq-B03", call))

write_csv_checked(calls, file.path(out, "work", "hashtag_calls.csv.gz"))
message("=== hashtag calls: ", nrow(calls), " cells over ",
        dplyr::n_distinct(calls$tube), " tubes; ",
        sum(calls$is_singlet), " singlet calls, ",
        sum(!calls$is_singlet), " multiplet/negative")
print(as.data.frame(dplyr::count(calls, tube, replicate, name = "cells") |>
                      tidyr::pivot_wider(names_from = replicate, values_from = cells,
                                         values_fill = 0L)), row.names = FALSE)
