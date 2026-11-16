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
  # The provirus's OTHER end. Reverse transcription makes both LTRs identical, and the contig keeps
  # one copy, so a read crossing the host -> 5' LTR junction aligns here and clips host to its LEFT.
  # 2736 is the processed 5' terminus: the contig reads ...GGGGGG AC TGGAAGG..., and integrase
  # removes that AC, leaving TG - the partner of the CA that ends the contig at 3236.
  att            = 2736L,
  # Unlike the 3' end, anchors here SMEAR. The contig continues upstream of 2736, so host bases at
  # the clip's 3' end that happen to match it are absorbed into the alignment by greedy extension and
  # the reported start slides left. A GGGGGG run at 2728-2733 makes this common: most clips lose one
  # base, some five. Accept [att - att_window, att] and add the absorbed bases back; see
  # read_five_prime_clips(). 2736 is a hard upper bound - a later start is not this junction.
  att_window     = 10L,
  min_discovery  = 25L,     # clip long enough to place on its own, genome-wide
  # Assignment is a different question from discovery and gets a different threshold. Discovery asks
  # "is there an integration anywhere in 2.7 Gb?"; assignment asks "does this clip match one of the
  # handful of flanks already called?". The second search space is ~3.6 kb, and the clip must land at
  # one exact offset in it, so a far shorter clip is decisive. 8 bp placed against the mini-reference
  # is safer than 25 bp placed genome-wide: with it, both wild-type tubes still yield zero cells.
  min_assign     = 8L,
  flank          = 200L,    # bases either side of a called site in the mini-reference
  # A flank containing a long homopolymer attracts clips that belong nowhere. JH584304.1:35,417 has a
  # 23-base run against 5-6 for every genuine site, and the mini-reference inflated it 3 -> 14 cells,
  # all in mCherry-NEGATIVE tubes and none in the two high-signal ones. Such a site is still reported
  # from discovery, but it is kept out of the mini-reference so it cannot collect assignments.
  max_homopolymer = 12L,
  max_softclip   = 2L,      # an assignment clip must align end to end, not partially
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


#' Every left-soft-clipped read anchored at the provirus's 5' terminus, with the clip repaired.
#'
#' The second detection channel. A host transcript reading into the provirus, or the provirus's own
#' upstream flank, gives a read of the form [host][5' LTR...] which aligns from the att and clips the
#' host to its left. Because both proviral LTRs are identical and the contig keeps one copy, that
#' read lands on the contig's 3' LTR - the same trick that makes the psi artefact look like a
#' junction, used deliberately here.
#'
#' Repairing the clip is the whole difficulty. Greedy extension absorbs host bases that happen to
#' match the contig upstream of the att, so the reported start slides left by `absorbed` bases and
#' the clip is that many bases too short. Those bases are genuinely ambiguous to the aligner - they
#' match both sequences - and are resolved by biology instead: the vector-side boundary is known to
#' be `att`, so everything before it is host. Skipping this step puts junctions tens of bases out and
#' makes the TSD check meaningless.
read_five_prime_clips <- function(bam, tube, opt = IS_DEFAULTS) {
  stopifnot(file.exists(bam))
  # Known cosmetic noise: this call emits "Warning: stack imbalance in '<-'" a few times per run.
  # It is R's evaluation-stack bookkeeping, raised at C level inside the Rsamtools/GenomicAlignments
  # read, not a condition this code can catch or a sign of bad output. It does not reproduce when
  # readGAlignments is called with the same arguments outside this file, and the results are
  # verified: 113 clips on Rx_4, 71 of them discovery-length, and the TSD check passes at 4 bp for
  # IS01, IS03, IS04 and IS08. Do not spend time on it.
  lo  <- opt$att - opt$att_window
  win <- GenomicRanges::GRanges(opt$contig, IRanges::IRanges(lo, opt$att))
  ga  <- GenomicAlignments::readGAlignments(
    bam,
    param = Rsamtools::ScanBamParam(which = win, what = "seq", tag = c("CB", "UB"))
  )
  if (length(ga) == 0L) return(empty_clip_table(tube))

  cigar <- GenomicAlignments::cigar(ga)
  start <- GenomicRanges::start(ga)
  keep  <- !is.na(cigar) & grepl("^[0-9]+S", cigar) &
    start <= opt$att & start >= opt$att - opt$att_window
  if (!any(keep)) return(empty_clip_table(tube))

  ga    <- ga[keep]
  cigar <- cigar[keep]
  meta  <- S4Vectors::mcols(ga)
  seqs  <- as.character(meta$seq)

  anchor   <- GenomicRanges::start(ga)
  absorbed <- opt$att - anchor                        # host bases swallowed by greedy extension
  clip_len <- as.integer(sub("^([0-9]+)S.*", "\\1", cigar)) + absorbed

  # Drop before building the table, not after. Filtering the tibble on `seqs` would make dplyr's
  # data mask evaluate a vector from this function's scope, which works but leaves R's evaluation
  # stack unbalanced and emits "stack imbalance in '<-'" on every call.
  ok <- clip_len <= nchar(seqs)                       # repair cannot exceed the read
  if (!any(ok)) return(empty_clip_table(tube))
  anchor <- anchor[ok]; clip_len <- clip_len[ok]; meta <- meta[ok, ]
  clip   <- substring(seqs[ok], 1L, clip_len)

  dplyr::tibble(
    tube       = tube,
    anchor_end = as.integer(anchor),                  # where bwa put it, before repair
    clip_len   = as.integer(clip_len),
    clip       = clip,
    polya_frac = nchar(gsub("[^A]", "", clip)) / pmax(clip_len, 1L),
    barcode    = na_if_empty(meta$CB),
    umi        = na_if_empty(meta$UB)
  )
}

#' Label each 5' clip with the one reason it is kept or dropped.
#'
#' Mirrors classify_clips(), minus the anchor test (the scan window already enforces it) and the psi
#' screen (that artefact presents at the 3' end only).
classify_five_prime_clips <- function(clips, opt = IS_DEFAULTS) {
  clips |> dplyr::mutate(class = dplyr::case_when(
    polya_frac >= opt$polya_fraction ~ "polyA tail",
    is.na(barcode) | is.na(umi)      ~ "no cell barcode or UMI",
    clip_len >= opt$min_discovery    ~ "discovery",
    clip_len >= opt$min_assign       ~ "assign",
    TRUE                             ~ "clip too short to place"
  ))
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
    clip_len >= opt$min_assign            ~ "assign",
    TRUE                                  ~ "clip too short to place"
  ))
}

#' Per-tube counts of every class, in a fixed order, for the report's waterfall.
clip_waterfall <- function(classified) {
  levels <- c("anchor not at the junction", "vector (5' LTR to psi)", "polyA tail",
              "no cell barcode or UMI", "clip too short to place", "assign", "discovery")
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
  d <- dplyr::filter(classified, class %in% tier)
  if (nrow(d) == 0L) {
    file.create(path)
    return(invisible(0L))
  }
  names <- sprintf(">%s|%s|%s|%d|%d", d$tube, d$barcode, d$umi, seq_len(nrow(d)), d$clip_len)
  writeLines(as.vector(rbind(names, d$clip)), path)
  invisible(nrow(d))
}


# ── Stage 2: place the clips and call sites ───────────────────────────────────────────────────────

#' Primary, uniquely-placed alignments from a bwa SAM, with the cell each clip came from.
#'
#' Secondary (0x100) and supplementary (0x800) records are dropped so each clip contributes once.
read_bwa_hits <- function(sam, opt = IS_DEFAULTS, channel = c("three_prime", "five_prime")) {
  channel <- match.arg(channel)
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
  #
  # The 5' channel is the mirror image. There the host PRECEDES the provirus, so the
  # junction-adjacent base is the clip's LAST one, and the two strand cases swap over.
  far_end <- bwa_pos + ref_span(field(6)[keep]) - 1L
  junction <- if (channel == "three_prime") ifelse(rev, far_end, bwa_pos)
              else                          ifelse(rev, bwa_pos, far_end)

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

#' Write the mini-reference: the flanks of the called sites, from the genome FASTA.
#'
#' Phase two of a two-question design. Discovery searched 2.7 Gb and needed a clip long enough to be
#' unique there. Assignment asks only whether a clip matches one of the flanks already called, so the
#' search space collapses to a few kb and the required evidence collapses with it.
#'
#' Built from DISCOVERY-tier sites only. Seeding it with anything the assignment step produced would
#' make the method confirm itself.
write_site_flanks <- function(sites, genome_fa, path, opt = IS_DEFAULTS) {
  stopifnot(file.exists(genome_fa), nrow(sites) > 0L)
  fa  <- Rsamtools::FaFile(genome_fa)
  gr  <- GenomicRanges::GRanges(
    sites$chrom,
    IRanges::IRanges(pmax(sites$pos - opt$flank, 1L), sites$pos + opt$flank))
  seqs <- Rsamtools::scanFa(fa, gr)

  hp   <- max_homopolymer(as.character(seqs))
  drop <- hp > opt$max_homopolymer
  if (any(drop)) {
    message("    excluded from the mini-reference (homopolymer > ", opt$max_homopolymer, "): ",
            paste(sprintf("%s (%d)", sites$site_id[drop], hp[drop]), collapse = ", "))
    seqs  <- seqs[!drop]
    sites <- sites[!drop, , drop = FALSE]
  }
  if (length(seqs) == 0L) stop("every called site's flank is low-complexity; nothing to assign to")
  # name carries where the junction sits inside the flank, so read_mini_hits() need not recompute it
  names(seqs) <- sprintf("%s|%s|%d", sites$site_id, sites$strand,
                         sites$pos - pmax(sites$pos - opt$flank, 1L) + 1L)
  Biostrings::writeXStringSet(seqs, path)
  invisible(length(seqs))
}

#' Longest run of a single base in each sequence.
max_homopolymer <- function(x) {
  vapply(regmatches(x, gregexpr("(.)\\1*", x)),
         function(runs) max(nchar(runs)), numeric(1))
}

#' Assign clips to sites by alignment against the mini-reference.
#'
#' The discriminator is POSITION, not alignment score. A clip belongs to a site only if its
#' junction-adjacent base lands on that site's junction base inside the flank - offset `j` in the
#' sequence name. The clip's first base is the junction one (the 3' channel), so a forward alignment
#' must START at j and a reverse one must END there. That single constraint is what makes an 8 bp
#' clip safe: it has 400 places to land in each flank and only one is accepted, and in practice both
#' wild-type tubes return zero cells under it.
#'
#' Alignment quality is then only a sanity filter: end-to-end, few mismatches.
read_mini_hits <- function(sam, opt = IS_DEFAULTS) {
  if (!file.exists(sam) || file.size(sam) == 0L) return(empty_hit_table())
  lines <- readLines(sam, warn = FALSE)
  lines <- lines[!startsWith(lines, "@")]
  if (length(lines) == 0L) return(empty_hit_table())

  f     <- strsplit(lines, "\t", fixed = TRUE)
  field <- function(i) vapply(f, function(x) if (length(x) >= i) x[[i]] else NA_character_, "")
  flag  <- as.integer(field(2))
  rname <- field(3)

  keep <- bitwAnd(flag, 0x100) == 0L & bitwAnd(flag, 0x800) == 0L & rname != "*"
  if (!any(keep)) return(empty_hit_table())

  f <- f[keep]; flag <- flag[keep]
  field <- function(i) vapply(f, function(x) if (length(x) >= i) x[[i]] else NA_character_, "")
  cigar <- field(6)
  pos   <- as.integer(field(4))
  rev   <- bitwAnd(flag, 0x10) != 0L
  ref   <- do.call(rbind, strsplit(rname[keep], "|", fixed = TRUE))
  junc  <- as.integer(ref[, 3])

  # NM:i: and the total soft-clipped length, from the optional fields
  tags <- vapply(f, function(x) paste(x[-(1:11)], collapse = "\t"), "")
  nm   <- suppressWarnings(as.integer(sub(".*NM:i:([0-9]+).*", "\\1", tags)))
  nm[is.na(nm)] <- 99L
  soft <- vapply(regmatches(cigar, gregexpr("[0-9]+S", cigar)),
                 function(x) sum(as.integer(sub("S$", "", x))), numeric(1))

  prov <- do.call(rbind, strsplit(field(1), "|", fixed = TRUE))
  clip_len <- as.integer(prov[, 5])
  span     <- ref_span(cigar)

  at_junction <- ifelse(rev, pos + span - 1L == junc, pos == junc)
  good <- at_junction & soft <= opt$max_softclip & nm <= pmax(1L, clip_len %/% 10L)
  if (!any(good)) return(empty_hit_table())

  dplyr::tibble(
    chrom    = ref[good, 1],          # the site_id; the genomic coordinate comes from the site table
    pos      = junc[good],
    strand   = ref[good, 2],
    mapq     = as.integer(field(5))[good],
    tube     = prov[good, 1],
    barcode  = prov[good, 2],
    umi      = prov[good, 3],
    clip_len = clip_len[good]
  ) |>
    dplyr::rename(site_id = chrom)
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
  # Ties are broken by locus, not left to grouping order. Three sites sat at 3 cells and the labels
  # reshuffled between runs, silently repointing IS07-IS09 at different loci and invalidating every
  # document that named them.
  order <- hits |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(cells = dplyr::n_distinct(paste(tube, barcode)),
                     chrom = dplyr::first(chrom), pos = min(pos), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(cells), chrom, pos) |>
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

#' Replace a site table's read/cell/UMI counts with the ones from the assignment stage.
#'
#' The site LIST and its coordinates come from discovery; the cell COUNTS come from assignment
#' against the mini-reference, which reaches far more cells. Without this the site table and
#' cell_site_assignments.csv disagree about the same quantity - the table reporting only the cells
#' that happened to carry a 25 bp clip.
update_site_counts <- function(site_table, mini, tubes) {
  if (nrow(site_table) == 0L) return(site_table)
  per_site <- mini |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(reads = dplyr::n(),
                     cells = dplyr::n_distinct(paste(tube, barcode)),
                     umis  = dplyr::n_distinct(paste(tube, barcode, umi)), .groups = "drop")
  per_tube <- mini |>
    dplyr::group_by(site_id, tube) |>
    dplyr::summarise(n = dplyr::n_distinct(barcode), .groups = "drop") |>
    tidyr::pivot_wider(names_from = tube, values_from = n, values_fill = 0L,
                       names_prefix = "cells_")
  for (t in setdiff(paste0("cells_", tubes), names(per_tube))) per_tube[[t]] <- 0L

  site_table |>
    dplyr::select(-dplyr::any_of(c("reads", "cells", "umis")),
                  -dplyr::starts_with("cells_")) |>
    dplyr::left_join(per_site, by = "site_id") |>
    dplyr::left_join(per_tube, by = "site_id") |>
    dplyr::mutate(dplyr::across(c(reads, cells, umis, dplyr::starts_with("cells_")),
                                ~ tidyr::replace_na(.x, 0L))) |>
    dplyr::arrange(dplyr::desc(cells), chrom, pos)
}

#' Match 5' junctions to called sites and test the target-site duplication.
#'
#' Integrase cuts the two host strands with a 5 bp stagger, and repair fills the overhang, so the
#' duplicated bases appear on both flanks of every provirus. The two junctions therefore share a 5 bp
#' span, which means their coordinates differ by 4 - mind that off-by-one before judging a site.
#'
#' This is the strongest check available: two independent read populations, from opposite ends of the
#' provirus, agreeing to the base. Nothing that merely resembles a junction can produce it. It has
#' already earned its place by catching the minus-strand coordinate bug, which was invisible to every
#' filter in the 3' channel.
#'
#' A 5' junction with no site within `search` is reported with site_id NA: it may be a real site the
#' 3' channel missed, so it is surfaced rather than dropped.
tsd_check <- function(sites, five_hits, opt = IS_DEFAULTS, search = 2000L, expect = 4L) {
  if (nrow(five_hits) == 0L || nrow(sites) == 0L) return(dplyr::tibble())
  five <- cluster_hits(five_hits, opt) |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(chrom = dplyr::first(chrom), pos_5p = min(pos),
                     reads_5p = dplyr::n(),
                     cells_5p = dplyr::n_distinct(paste(tube, barcode)), .groups = "drop") |>
    dplyr::select(-site_id)

  # Nearest called site to each 5' junction, as a plain loop over a handful of rows. rowwise() with a
  # list-column is the idiomatic-looking alternative and it errors inside mutate(); this is clearer
  # and cannot.
  n       <- nrow(five)
  site_id <- rep(NA_character_, n)
  pos_3p  <- rep(NA_integer_,   n)
  delta   <- rep(NA_integer_,   n)
  for (i in seq_len(n)) {
    near <- which(sites$chrom == five$chrom[i] & abs(sites$pos - five$pos_5p[i]) <= search)
    if (!length(near)) next
    k <- near[which.min(abs(sites$pos[near] - five$pos_5p[i]))]
    site_id[i] <- sites$site_id[k]
    pos_3p[i]  <- as.integer(sites$pos[k])
    delta[i]   <- as.integer(abs(sites$pos[k] - five$pos_5p[i]))
  }

  five |>
    dplyr::mutate(site_id = site_id, pos_3p = pos_3p, delta = delta) |>
    # A TSD is 4 bp, so only a junction within a few bases of a site is a claim about that site.
    # Anything further is a distinct feature that merely happens to be nearby, and must not be
    # labelled FAIL against a site that passed on its own junction - chr19:38,932,385 sits 1.4 kb
    # from IS03 and would otherwise read as IS03 failing.
    dplyr::mutate(tsd = dplyr::case_when(
      is.na(delta)         ~ "no 3' site nearby",
      delta == expect      ~ "PASS",
      delta <= expect + 2L ~ "near",
      delta <= 50L         ~ "FAIL",
      TRUE                 ~ "distinct locus"
    )) |>
    dplyr::mutate(site_id = ifelse(tsd == "distinct locus", NA_character_, site_id)) |>
    dplyr::arrange(dplyr::desc(reads_5p)) |>
    dplyr::select(site_id, chrom, pos_3p, pos_5p, delta, tsd, reads_5p, cells_5p)
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
