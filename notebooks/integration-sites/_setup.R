# Shared setup for the four scarab reports.
#
# Split out of integration-sites.qmd when that report became four. Every report reads the same
# pipeline outputs and must colour a clone the same way, so this is one file rather than four
# copies - four copies of a palette is how IS01 ends up a different colour in two reports.
#
# Sourced from the project root, because that is where Quarto executes (execute_dir = getwd()
# in quarto_report()).

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(scales)
library(DT)
library(targets)
# The graph and counts targets are dgCMatrix. rownames() and `[` on them are S4 methods, so without
# Matrix's methods registered rownames() can come back NULL and an order check against it compares a
# character vector to NULL - which is FALSE, and reads as a cell-order mismatch rather than as a
# missing method. Attached rather than left to lazy loading, because whether a namespace happens to
# be loaded by that point depends on which chunks ran before.
library(Matrix)

library(scarab)   # the pipeline, as a package: /well/ludwig/users/nrs939/projects/scarab
# Was source("src/R/integration-sites.R"). That file was the only definition of this method and is
# now scarab, installed like any other dependency. A second copy in the project would guarantee
# drift between them - which happened twice while packaging, in both directions, and is why
# dev/compare-to-project.R exists in the package.
#
# IS_DEFAULTS is scarab_opts() now: a validated constructor rather than a bare list, so a
# misspelled override is an error instead of a silent no-op. The alias is kept because the report
# reads thresholds from it in a dozen places and `opt` is the pipeline-wide name for the object.
IS_DEFAULTS <- scarab_opts()

sites_dir <- "output/objects/integration-sites"
stopifnot(dir.exists(sites_dir))

anchor_ends <- readr::read_csv(file.path(sites_dir, "anchor_ends.csv"), show_col_types = FALSE)
waterfall   <- readr::read_csv(file.path(sites_dir, "clip_waterfall.csv"), show_col_types = FALSE)
sites       <- readr::read_csv(file.path(sites_dir, "integration_sites.csv"), show_col_types = FALSE)
assign      <- readr::read_csv(file.path(sites_dir, "cell_site_assignments.csv"), show_col_types = FALSE)
summary_tbl <- readr::read_csv(file.path(sites_dir, "site_calling_summary.csv"), show_col_types = FALSE)
tsd         <- readr::read_csv(file.path(sites_dir, "tsd_check.csv"), show_col_types = FALSE)
# The per-site evidence tier and `defensible`. Here rather than in one report's headline chunk,
# because two of the four need it: scarab-qc builds the headline from it and scarab-clones asks which
# clones are defensible. It was read in scarab-qc only, which made `ev` resolve while that report was
# the one being rendered and fail with "object 'ev' not found" in the other - and, worse, meant
# scarab-clones read a table it did not declare, so the provenance banner verified one fewer file
# than the report actually depends on.
ev          <- readr::read_csv(file.path(sites_dir, "site_evidence.csv"), show_col_types = FALSE)
coverage    <- readr::read_csv(file.path(sites_dir, "junction_coverage.csv"), show_col_types = FALSE)

# The molecule-level table and its companions. clip_tbl is ONE ROW PER MOLECULE with channel
# membership as logicals, which is why nothing below adds umis_three to umis_genomic - that sum
# double-counts the 83 molecules present in both channels.
clip_tbl        <- readr::read_csv(file.path(sites_dir, "clip_table.csv.gz"), show_col_types = FALSE)
detect_tbl      <- readr::read_csv(file.path(sites_dir, "detection_model.csv"), show_col_types = FALSE)
# Read here rather than in the completeness section: the clustering test needs it ~500 lines earlier,
# and a chunk that defines a variable below its first use renders only until someone reorders the
# document.
detr            <- readr::read_csv(file.path(sites_dir, "detection_regression.csv"),
                                   show_col_types = FALSE)
# The reference GTF, for the cis-loading check. One place, so the report and the pipeline cannot
# point at different annotations; the stopifnot is what makes that a check rather than a hope.
gtf_path <- file.path(Sys.getenv(
  "CR_REF",
  "/well/ludwig/users/nrs939/data_hub/10x_resources/reference_genomes/modified/GRCm39-2024-A_LeGO_SFFV_mCherry_DRE"),
  "genes", "genes.gtf.gz")
stopifnot(file.exists(gtf_path))
called_barcodes <- readr::read_csv(file.path(sites_dir, "work", "called_barcodes.csv.gz"),
                                   show_col_types = FALSE)
# Tube -> what it is, so figures can say "HSC Het mCh+" rather than "Rx_4".
tube_map <- readr::read_csv("metadata/tube_sample_map.csv", show_col_types = FALSE)

# One factor with every site as a level, so a clone keeps the same colour in every panel even where
# it has no cells. As a bare character it would be re-levelled per panel and the colours would shift.
assign$site_id <- factor(assign$site_id, levels = sites$site_id)

# Clone colours. These are nominal categories with no order, so a sampled continuous ramp (viridis,
# turbo) is the wrong tool: neighbouring clones come out as neighbouring hues. glasbey is built by
# maximising distance in CIELAB, and measurably wins here - the closest pair of colours is 67 dE
# apart against 28 for the project's .discrete_colours() at n = 9 and 24 for turbo. Its first entry
# is black, dropped so no clone is confused with text or an axis.
#
# The trade-off, stated because it is real: glasbey is not optimised for colour vision deficiency
# (closest pair under deuteranopia ~9 dE). Okabe-Ito and Tol muted do better there (~17) but are far
# worse for normal vision (~26 and ~23), and their low chroma washes out at this point size. If a
# reader needs CVD-safe output, swap this one line.
# Keyed on the site LABEL, not on its position in the table. The labels are pinned to loci by the
# registry, but the table is sorted by size, so keying the palette on position would move every
# clone's colour the moment one changed rank - the same failure the registry fixed for the labels,
# and harder to notice in a figure. Sorted ids mean IS01 keeps its colour and a new site appends.
pal_ids      <- sort(levels(assign$site_id))
site_palette <- stats::setNames(as.character(pals::glasbey(length(pal_ids) + 1))[-1], pal_ids)

# Taken from the pipeline rather than restated, so the figures cannot drift from the run.
JUNCTION <- IS_DEFAULTS$junction
MIN_DISC <- IS_DEFAULTS$min_discovery
MIN_ASSIGN <- IS_DEFAULTS$assign_floor$min_len
CLEAVAGE <- 3152L   # measured from the polyA clips; see reference_additions/README.md

theme_set(theme_minimal(base_size = 11))

# The hashtag levels, in panel order then the two non-calls. Written once because it appeared three
# times inline and a fourth copy is how one of them ends up missing a tag.
HASHTAGS <- c("TotalSeq-B0301", "TotalSeq-B0302", "TotalSeq-B0303", "TotalSeq-B0304",
              "negative", "multiplet")

# ── The cell tables ─────────────────────────────────────────────────────────────────────────────
# Three of the four reports draw clones on cells, so the join between the assignment table and the
# cell embedding lives here rather than in whichever report happened to need it first. It was in
# scarab-clones until scarab-qc and scarab-exploration both turned out to need `embedding` and
# `labelled` too.
#
# `build_cell_tables()` takes the embedding as an ARGUMENT instead of calling tar_read() itself.
# That is deliberate and not a style choice: tar_knitr_deps() scans the .qmd only, never a file it
# sources, so a tar_read() buried in here would be invisible to the pipeline and the report would
# build against a stale target without complaint. Each report keeps a literal
# `tar_read(is_cell_embedding)` in its own setup chunk, and passes the result in.

#' Cell names, sort label and UMAP coordinates for one processed compartment.
#'
#' @param emb_all the `is_cell_embedding` target: all four objects' coordinates in one tibble.
#' @return a small tibble, so the Seurat objects need never be read by a report at all.
cell_embedding <- function(emb_all, panel_name, stage_name = "post-QC") {
  emb_all |>
    dplyr::filter(panel == panel_name, stage == stage_name) |>
    dplyr::transmute(
      cell, panel,
      tube        = as.character(tube),
      sample_name = as.character(`sample.name`),
      sample      = as.character(`sample.name.long`),
      # the hashtag demultiplex call: which replicate mouse this cell came from. "multiplet" and
      # "negative" are failures of the call, not animals, and appear mostly before QC.
      # NOT the authoritative source. Replicate attribution is load-bearing in three places - the
      # cross-mouse co-occurrence control, the mono-replicate assignment check, and every B0304
      # caveat - so it is joined from hashtag_calls.csv.gz below instead. This column is kept so the
      # join has something to check against, which is also a check on the file.
      replicate_srt = factor(as.character(dm2_class), levels = HASHTAGS),
      compartment = as.character(compartment),
      pass_cell   = as.logical(pass_cell),
      UMAP_1, UMAP_2)
}

#' Attach `replicate` from the calls file, checking it against the object's own column.
#'
#' One helper for every embedding, because `cell_embedding()` has two callers and an earlier version
#' of this fixed up only the post-QC one - leaving the pre-QC embedding with no `replicate` column at
#' all, which the facets then failed on.
#'
#' The file is checked against the object in two parts rather than with one identical(). A first
#' version demanded identical() and halted the render: of 50,088 HSC cells present in both, ZERO
#' disagree, and a single cell is absent from the file - so the test was failing on coverage while
#' the content was perfect. The two failures mean different things and deserve different treatment.
#'
#'  - a DISAGREEMENT means one source is stale, most likely the file after an upstream target was
#'    rebuilt, and nothing downstream can be trusted: halt.
#'  - a GAP means the demultiplexer's output does not cover a cell the processed object kept. One
#'    such cell exists (its call is "negative"). Fill it from the object and say how many, but halt
#'    if it is more than a handful, because that would mean the file is from a different run.
#'
#' `strict` is TRUE for the called-cell embedding, where the demultiplexer saw essentially every cell
#' and a gap means the file is from another run. It is FALSE for the pre-QC embedding, which contains
#' barcodes that were never called cells and which the demultiplexer therefore never saw: a large gap
#' there is expected, not a fault.
attach_replicate <- function(emb, label, strict = TRUE) {
  both <- !is.na(emb$replicate) & !is.na(emb$replicate_srt)
  bad  <- sum(as.character(emb$replicate[both]) != as.character(emb$replicate_srt[both]))
  if (bad > 0L) {
    stop("hashtag calls disagree with the Seurat objects for ", bad, " of ", sum(both),
         " cells in ", label, ". One is stale; re-run scripts/R/export-hashtag-calls.R.",
         call. = FALSE)
  }
  gap <- sum(is.na(emb$replicate) & !is.na(emb$replicate_srt))
  if (strict && gap > 10L) {
    stop("hashtag_calls.csv.gz is missing ", gap, " of ", label, "'s cells, too many to be a ",
         "demultiplexer edge case; the file is probably from another run.", call. = FALSE)
  }
  if (gap > 0L) {
    message("note: ", gap, " cell(s) of ", label,
            " absent from hashtag_calls.csv.gz, filled from the processed object")
  }
  emb |>
    dplyr::mutate(replicate = dplyr::coalesce(replicate, replicate_srt)) |>
    dplyr::select(-replicate_srt)
}

#' The four cell tables every report that touches cells needs.
#'
#' @param emb_all the `is_cell_embedding` target, read literally by the calling .qmd.
#' @return a list of `embedding`, `labelled`, `embedding_raw`, `labelled_raw`, `hashtag_calls`.
#'   `list2env()` it into the report's environment.
build_cell_tables <- function(emb_all) {
  # The demultiplex calls, from a file rather than from the objects. The assertion in
  # attach_replicate() is the point: if the file and the Seurat column ever disagree, one of them is
  # stale and the report should say so rather than quietly pick one.
  hashtag_calls <- readr::read_csv(file.path(sites_dir, "work", "hashtag_calls.csv.gz"),
                                   show_col_types = FALSE)
  calls <- dplyr::transmute(hashtag_calls, cell,
                            replicate = factor(call, levels = HASHTAGS))

  # One read of one small target, then filtered per panel and stage. This used to be four literal
  # tar_read() calls on 2 GB objects, taken sequentially so they were never resident together - a
  # constraint that disappears when the thing being read is a few megabytes.
  embedding <- dplyr::bind_rows(cell_embedding(emb_all, "HSC"),
                                cell_embedding(emb_all, "MNC")) |>
    dplyr::left_join(calls, by = "cell") |>
    # A factor with all three levels, not a character: the wild-type samples contain no labelled
    # cells, so as a character a facet on this would silently drop the column that matters most.
    dplyr::mutate(sorted = factor(
      dplyr::case_when(
        grepl("_positive$", sample) ~ "FACS mCherry+",
        grepl("_negative$", sample) ~ "FACS mCherry−",
        grepl("wildtype",   sample) ~ "wild-type (no vector)"),
      levels = c("FACS mCherry+", "FACS mCherry−", "wild-type (no vector)"))) |>
    attach_replicate("the called-cell embedding", strict = TRUE)

  # The pre-QC embedding. The pre-QC UMAP figure that used to accompany it was dropped as a
  # near-duplicate of the post-QC one, but these two tables are not duplicates - three tables in
  # scarab-qc count labelled barcodes against `pass_cell`, and removing the figure took the data
  # with it until a chunk-level smoke test caught the orphan.
  embedding_raw <- dplyr::bind_rows(cell_embedding(emb_all, "HSC", "pre-QC"),
                                    cell_embedding(emb_all, "MNC", "pre-QC")) |>
    dplyr::left_join(calls, by = "cell") |>
    attach_replicate("the pre-QC embedding", strict = FALSE)

  # tube + barcode is the Seurat cell name, e.g. Rx_4 + AAACCAAAGAATCATC-1 -> Rx_4_AAACCAAAGAATCATC-1
  join_assign <- function(emb) {
    assign |>
      dplyr::mutate(cell = paste0(tube, "_", barcode)) |>
      dplyr::inner_join(dplyr::select(emb, -tube), by = "cell")
  }

  list(hashtag_calls = hashtag_calls,
       embedding     = embedding,
       labelled      = join_assign(embedding),
       embedding_raw = embedding_raw,
       labelled_raw  = join_assign(embedding_raw))
}
