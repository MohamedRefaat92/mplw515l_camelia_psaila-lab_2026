# Lentiviral integration sites from read-through transcripts.
#
# Transcription of the cassette does not stop cleanly at the polyA site in the 3' LTR's R region: a
# minority of transcripts read through to the end of the provirus and on into flanking host DNA. Such
# a read aligns to the reporter contig as far as its last base and soft-clips the rest, and that
# clipped tail is host sequence naming where the vector integrated. The read still carries its cell
# barcode, so the locus arrives attached to a cell.
#
# Functions only, no side effects on load: source() this and step through it, or let
# scripts/R/call-integration-sites.R drive it. The alignment of the clips happens between the two
# stages, in scripts/slurm/analysis/integration-sites.sh, because it needs bwa.
#
#   source("src/R/integration-sites.R")
#   clips <- read_junction_clips("data/cellranger_Rx_4/sample_alignments.bam", "Rx_4")
#   table(classify_clips(clips)$class)

# Attaches nothing: every call below is namespace-qualified and the pipes are base |>, so sourcing
# this file has no side effects and the caller decides what to put on its search path. The packages
# it reaches into, all in renv:
#
#   Rsamtools GenomicAlignments GenomicRanges IRanges S4Vectors   reading the BAM, clustering
#   rtracklayer                                                   reading the gene GTF
#   dplyr tidyr                                                   the tables
#
# This matters for the report, which sources the file only for IS_DEFAULTS: attaching here would
# drag the whole Bioconductor stack into a session that never touches a BAM.

# The contig and the one position on it that matters. The junction is the contig's LAST base: the
# contig ends ...TCTCTAGCA, and that CA is what integrase leaves after 3' processing trims the
# terminal GT from the blunt reverse-transcription product. Every provirus is joined to the host
# genome at that CA, so a clip anchored anywhere else cannot be an integration junction. Re-cut the
# cassette with a different 3' boundary and this number must move with it.
IS_DEFAULTS <- list(
  contig         = "LeGO_SFFV_mCherry_DRE",
  junction       = 3236L,   # the processed CA; the integration junction itself
  scan_from      = 3150L,   # anchors can only start upstream of it
  min_discovery  = 25L,     # clip long enough to place on its own
  min_support    = 15L,     # clip that may only add cells to a site already found
  min_mapq       = 30L,     # bwa: a unique placement
  cluster_window = 20L,     # hits this close are one site
  min_cells      = 3L,
  min_umis       = 3L,
  polya_fraction = 0.8,
  # The plasmid's two LTRs are byte-identical over all 181 bases and the contig keeps only one copy,
  # so a molecule spanning the *5'* LTR's junction with psi matches the contig's 3' terminus
  # perfectly: it is placed at the junction with MAPQ 255 and clips into vector sequence,
  # structurally identical to a host junction. No MAPQ or coverage filter sees it, because the anchor
  # genuinely is a unique alignment. Only the clip's sequence gives it away.
  psi_prefix     = "GTGGCGCCCGAACAGG"
)


# ── Stage 1: pull the clipped tails out of one BAM ────────────────────────────────────────────────

#' Every right-soft-clipped read whose anchor lies in the scan window of the reporter contig.
#'
#' One pass over a ~90-base window of one contig, so this is seconds despite the BAM's size.
#' Returns one row per read, with the anchor's reference end and the clipped tail, and no filtering
#' yet: classify_clips() decides what each row is, so the discards stay countable.
read_junction_clips <- function(bam, tube, opt = IS_DEFAULTS) {
  stopifnot(file.exists(bam))
  # readGAlignments does the CIGAR arithmetic itself, so end() is the anchor's reference end and no
  # hand-rolled parser is needed. seq and the barcode tags ride along in mcols().
  ga <- GenomicAlignments::readGAlignments(
    bam,
    param = Rsamtools::ScanBamParam(
      which = GenomicRanges::GRanges(opt$contig,
                                     IRanges::IRanges(opt$scan_from, opt$junction)),
      what  = "seq",
      tag   = c("CB", "UB")
    )
  )
  if (length(ga) == 0L) return(empty_clip_table(tube))

  cigar <- GenomicAlignments::cigar(ga)
  keep  <- !is.na(cigar) & grepl("[0-9]+S$", cigar)      # a right soft clip, i.e. a tail to recover
  if (!any(keep)) return(empty_clip_table(tube))

  ga    <- ga[keep]
  cigar <- cigar[keep]
  meta  <- S4Vectors::mcols(ga)
  # Soft-clipped bases are kept in SEQ (unlike hard clips), which is why the tail is recoverable.
  seqs     <- as.character(meta$seq)
  clip_len <- as.integer(sub(".*?([0-9]+)S$", "\\1", cigar))
  clip     <- substring(seqs, nchar(seqs) - clip_len + 1L, nchar(seqs))

  dplyr::tibble(
    tube       = tube,
    anchor_end = as.integer(GenomicRanges::end(ga)),
    clip_len   = clip_len,
    clip       = clip,
    # the clipped tail's adenine fraction: the transcript's own polyA tail, where this is high
    polya_frac = nchar(gsub("[^A]", "", clip)) / pmax(clip_len, 1L),
    barcode    = na_if_empty(meta$CB),
    umi        = na_if_empty(meta$UB)
  )
}

empty_clip_table <- function(tube) {
  dplyr::tibble(tube = character(), anchor_end = integer(), clip_len = integer(),
                clip = character(), polya_frac = numeric(),
                barcode = character(), umi = character())[0, ] |>
    dplyr::mutate(tube = as.character(tube))
}

na_if_empty <- function(x) {
  if (is.null(x)) return(NA_character_)
  x[!nzchar(x)] <- NA_character_
  x
}


#' Label each clip with the one reason it is kept or dropped.
#'
#' The order is deliberate and the labels are exclusive, so the waterfall adds up and every discard
#' is attributable. Anchor position comes first because it is the only criterion that is a statement
#' about enzymology rather than a threshold: it removes the large polyA population at the cleavage
#' site, which carries good MAPQ, real barcodes and real UMIs and so would otherwise survive
#' everything downstream and be counted as sites.
classify_clips <- function(clips, opt = IS_DEFAULTS) {
  clips |> dplyr::mutate(class = dplyr::case_when(
    anchor_end != opt$junction            ~ "anchor not at the junction",
    startsWith(clip, opt$psi_prefix)      ~ "vector (5' LTR to psi)",
    polya_frac >= opt$polya_fraction      ~ "polyA tail",
    is.na(barcode) | is.na(umi)           ~ "no cell barcode or UMI",
    clip_len >= opt$min_discovery         ~ "discovery",
    clip_len >= opt$min_support           ~ "support",
    TRUE                                  ~ "clip too short to place"
  ))
}

#' Per-tube counts of every class, in a fixed order, for the report's waterfall.
clip_waterfall <- function(classified) {
  levels <- c("anchor not at the junction", "vector (5' LTR to psi)", "polyA tail",
              "no cell barcode or UMI", "clip too short to place", "support", "discovery")
  classified |>
    dplyr::count(tube, class, name = "reads") |>
    dplyr::mutate(class = factor(class, levels = levels)) |>
    dplyr::arrange(tube, class)
}

#' Anchor end against clip length, for the report's coverage view.
#'
#' Every clip, whatever its class: the figure's whole point is to show the polyA pile at the cleavage
#' site next to the junction pile, so filtering here would hide what it is meant to demonstrate.
anchor_end_histogram <- function(classified) {
  classified |>
    dplyr::mutate(polya = polya_frac >= IS_DEFAULTS$polya_fraction) |>
    dplyr::count(tube, anchor_end, clip_len, polya, name = "reads")
}

#' Write one tier's clips as FASTA for bwa.
#'
#' The read name carries provenance the SAM cannot: tube|barcode|umi|serial. The serial only makes
#' names unique, since one cell can contribute several clips of the same UMI.
write_clip_fasta <- function(classified, tier, path) {
  d <- dplyr::filter(classified, class == tier)
  if (nrow(d) == 0L) {
    file.create(path)
    return(invisible(0L))
  }
  names <- sprintf(">%s|%s|%s|%d", d$tube, d$barcode, d$umi, seq_len(nrow(d)))
  writeLines(as.vector(rbind(names, d$clip)), path)
  invisible(nrow(d))
}


# ── Stage 2: place the clips and call sites ───────────────────────────────────────────────────────

#' Primary, uniquely-placed alignments from a bwa SAM, with the cell each clip came from.
#'
#' Secondary (0x100) and supplementary (0x800) records are dropped so each clip contributes once.
read_bwa_hits <- function(sam, opt = IS_DEFAULTS) {
  if (!file.exists(sam) || file.size(sam) == 0L) return(empty_hit_table())
  lines <- readLines(sam, warn = FALSE)
  lines <- lines[!startsWith(lines, "@")]
  if (length(lines) == 0L) return(empty_hit_table())

  f     <- strsplit(lines, "\t", fixed = TRUE)
  field <- function(i) vapply(f, function(x) if (length(x) >= i) x[[i]] else NA_character_, "")
  flag  <- as.integer(field(2))
  mapq  <- as.integer(field(5))
  chrom <- field(3)

  keep <- bitwAnd(flag, 0x100) == 0L & bitwAnd(flag, 0x800) == 0L &
    chrom != "*" & !is.na(mapq) & mapq >= opt$min_mapq
  if (!any(keep)) return(empty_hit_table())

  provenance <- do.call(rbind, strsplit(field(1)[keep], "|", fixed = TRUE))
  rev <- bitwAnd(flag[keep], 0x10) != 0L
  bwa_pos <- as.integer(field(4)[keep])

  # The junction is NOT bwa's POS on the minus strand.
  #
  # A 3' clip is the host sequence FOLLOWING the provirus, so in the read it runs
  # [provirus ... 3236][h1 h2 ... hL] and the base abutting the integration junction is h1, the
  # clip's *first* base. bwa reports the leftmost reference coordinate of the alignment, which is
  # h1 only when the clip aligned forward; when it aligned reverse-complemented, h1 sits at the
  # alignment's far end instead. Taking POS regardless puts minus-strand sites a whole clip length
  # away from their true junction -- about 48 bp here, and silently, because every clip at one site
  # shares a strand so the cluster still looks tight.
  #
  # Caught by the 5' target-site-duplication check: plus-strand sites came out at the expected 4 bp
  # offset while minus-strand ones sat at 45. Do not "simplify" this back to POS.
  junction <- ifelse(rev, bwa_pos + ref_span(field(6)[keep]) - 1L, bwa_pos)

  dplyr::tibble(
    chrom    = chrom[keep],
    pos      = as.integer(junction),
    strand   = ifelse(rev, "-", "+"),
    mapq     = mapq[keep],
    tube     = provenance[, 1],
    barcode  = provenance[, 2],
    umi      = provenance[, 3],
    clip_len = nchar(field(10)[keep])
  )
}

#' Bases of the reference an alignment spans, from its CIGAR.
#'
#' M/D/N/=/X consume reference, I/S/H/P do not. Written out rather than calling
#' GenomicAlignments::cigarWidthAlongReferenceSpace(), which is deprecated.
ref_span <- function(cigar) {
  ops <- regmatches(cigar, gregexpr("[0-9]+[MDN=X]", cigar))
  vapply(ops, function(x) sum(as.integer(sub("[MDN=X]$", "", x))), numeric(1))
}

empty_hit_table <- function() {
  dplyr::tibble(chrom = character(), pos = integer(), strand = character(), mapq = integer(),
                tube = character(), barcode = character(), umi = character(), clip_len = integer())
}

#' Group hits into sites, adding a site_id column.
#'
#' GenomicRanges::reduce does the single-linkage merge. Junctions of one provirus scatter over a few
#' bases, because a clip starts wherever its read happened to begin, so a fixed-width bin would split
#' one site in two; reduce with min.gapwidth merges by distance instead.
cluster_hits <- function(hits, opt = IS_DEFAULTS) {
  if (nrow(hits) == 0L) return(dplyr::mutate(hits, site_id = character()))
  gr      <- GenomicRanges::GRanges(hits$chrom, IRanges::IRanges(hits$pos, width = 1L))
  merged  <- GenomicRanges::reduce(gr, min.gapwidth = opt$cluster_window + 1L)
  cluster <- GenomicRanges::findOverlaps(gr, merged, select = "first")
  hits |> dplyr::mutate(site_id = sprintf("c%04d", cluster))
}

#' Sites supported by enough independent cells and UMIs, from the discovery tier alone.
call_sites <- function(discovery_hits, opt = IS_DEFAULTS) {
  cluster_hits(discovery_hits, opt) |>
    dplyr::group_by(site_id) |>
    dplyr::filter(dplyr::n_distinct(paste(tube, barcode)) >= opt$min_cells,
                  dplyr::n_distinct(paste(tube, barcode, umi)) >= opt$min_umis) |>
    dplyr::ungroup()
}

#' Join the support tier onto sites discovery already found; discard the rest.
#'
#' A support clip is too short to place on its own, so it may add a cell to a known site but never
#' define one. Measured on the two large transduced tubes, every support clip that passed the MAPQ
#' filter landed on a site the long clips had already found, which is what makes this safe: it
#' recovers roughly a third more cells at the large sites with no way to invent a locus.
attach_support <- function(sites, support_hits, opt = IS_DEFAULTS) {
  if (nrow(support_hits) == 0L || nrow(sites) == 0L) {
    return(dplyr::mutate(sites, tier = "discovery"))
  }
  anchors <- sites |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(chrom = dplyr::first(chrom), pos = min(pos), .groups = "drop")

  near <- GenomicRanges::distanceToNearest(
    GenomicRanges::GRanges(support_hits$chrom, IRanges::IRanges(support_hits$pos, width = 1L)),
    GenomicRanges::GRanges(anchors$chrom, IRanges::IRanges(anchors$pos, width = 1L))
  )
  joined <- support_hits[S4Vectors::queryHits(near), ] |>
    dplyr::mutate(site_id  = anchors$site_id[S4Vectors::subjectHits(near)],
                  distance = S4Vectors::mcols(near)$distance) |>
    dplyr::filter(distance <= opt$cluster_window) |>
    dplyr::select(-distance)

  dplyr::bind_rows(dplyr::mutate(sites, tier = "discovery"),
                   dplyr::mutate(joined, tier = "support"))
}


# ── Annotation and output ─────────────────────────────────────────────────────────────────────────

#' Gene the junction falls in, else the nearest one.
#'
#' Where genes overlap, a protein-coding one wins over a lncRNA and then the shorter one, so a
#' junction inside a gene nested in a long lncRNA is reported against the gene.
annotate_positions <- function(chrom, pos, gtf) {
  genes <- rtracklayer::import(gtf, format = "gtf", feature.type = "gene")
  query <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, width = 1L))

  ov <- GenomicRanges::findOverlaps(query, genes)
  best <- dplyr::tibble(
    i      = S4Vectors::queryHits(ov),
    j      = S4Vectors::subjectHits(ov),
    coding = genes$gene_type[S4Vectors::subjectHits(ov)] == "protein_coding",
    width  = GenomicRanges::width(genes)[S4Vectors::subjectHits(ov)]
  ) |>
    dplyr::arrange(i, dplyr::desc(coding), width) |>
    dplyr::distinct(i, .keep_all = TRUE)

  out <- dplyr::tibble(gene = NA_character_, gene_type = NA_character_,
                       gene_strand = NA_character_, context = NA_character_,
                       .rows = length(query))
  out$gene[best$i]        <- genes$gene_name[best$j]
  out$gene_type[best$i]   <- genes$gene_type[best$j]
  out$gene_strand[best$i] <- as.character(GenomicRanges::strand(genes))[best$j]
  out$context[best$i]     <- "intragenic"

  missing <- which(is.na(out$context))
  if (length(missing)) {
    near <- GenomicRanges::distanceToNearest(query[missing], genes)
    k    <- S4Vectors::subjectHits(near)
    kb   <- S4Vectors::mcols(near)$distance %/% 1000L
    rows <- missing[S4Vectors::queryHits(near)]
    out$gene[rows]        <- genes$gene_name[k]
    out$gene_type[rows]   <- genes$gene_type[k]
    out$gene_strand[rows] <- as.character(GenomicRanges::strand(genes))[k]
    out$context[rows]     <- sprintf("intergenic, %d kb away", kb)
  }
  out
}

#' Replace the clustering's internal ids with stable IS01.. labels, largest site first.
#'
#' Done once, before either summary, so the site table and the per-cell table cannot disagree about
#' which site is which.
label_sites <- function(hits) {
  if (nrow(hits) == 0L) return(hits)
  order <- hits |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(cells = dplyr::n_distinct(paste(tube, barcode)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(cells)) |>
    dplyr::mutate(label = sprintf("IS%02d", dplyr::row_number()))
  hits |>
    dplyr::left_join(dplyr::select(order, site_id, label), by = "site_id") |>
    dplyr::mutate(site_id = label) |>
    dplyr::select(-label)
}

#' One row per site: position, host gene, and cell counts overall and per tube.
summarise_sites <- function(hits, gtf, tubes) {
  if (nrow(hits) == 0L) return(dplyr::tibble())
  per_site <- hits |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(
      chrom  = dplyr::first(chrom),
      pos    = min(pos),
      strand = paste(sort(unique(strand)), collapse = "/"),
      reads  = dplyr::n(),
      cells  = dplyr::n_distinct(paste(tube, barcode)),
      umis   = dplyr::n_distinct(paste(tube, barcode, umi)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(cells))

  per_tube <- hits |>
    dplyr::group_by(site_id, tube) |>
    dplyr::summarise(cells = dplyr::n_distinct(barcode), .groups = "drop") |>
    tidyr::pivot_wider(names_from = tube, values_from = cells,
                       values_fill = 0L, names_prefix = "cells_")
  for (t in setdiff(paste0("cells_", tubes), names(per_tube))) per_tube[[t]] <- 0L

  dplyr::bind_cols(per_site, annotate_positions(per_site$chrom, per_site$pos, gtf)) |>
    dplyr::left_join(per_tube, by = "site_id") |>
    dplyr::select(site_id, chrom, pos, strand, gene, gene_type, gene_strand, context,
                  reads, cells, umis, dplyr::starts_with("cells_"))
}

#' One row per cell per site: the table that joins onto the Seurat metadata.
#'
#' A cell at more than one site is flagged, not resolved: a genuine multi-copy integrant and a
#' doublet are indistinguishable here, so the choice belongs downstream.
summarise_cells <- function(hits) {
  if (nrow(hits) == 0L) return(dplyr::tibble())
  hits |>
    dplyr::group_by(site_id, tube, barcode) |>
    dplyr::summarise(reads = dplyr::n(),
                     umis = dplyr::n_distinct(umi),
                     max_clip_len = max(clip_len),
                     tiers = paste(sort(unique(tier)), collapse = "+"),
                     .groups = "drop") |>
    dplyr::group_by(tube, barcode) |>
    dplyr::mutate(n_sites_in_cell = dplyr::n_distinct(site_id),
                  ambiguous = n_sites_in_cell > 1L) |>
    dplyr::ungroup() |>
    dplyr::arrange(tube, barcode, site_id)
}
