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
  # one exact offset in it, so a far shorter clip is decisive: 10 bp against the mini-reference is a
  # stricter test than 25 bp placed genome-wide, not a looser one.
  #
  # Set to 6 by the decoy measurement, one base at a time, and only in company with exact_below.
  # The clip table carries every clip against every flank at every mismatch count with no floor, so
  # the 100 decoy flanks answer this directly: a decoy has no junction, so a decoy match at a given
  # length IS the chance rate. Under exact matching they are clean from 6 bp up - 28 molecules at
  # 6 bp, all matching exactly one real site, none matching any of the 100 decoys - and dirty below:
  # 1 of 11 molecules hits a decoy at 5 bp, 10 of 17 at 4 bp, and at 1-2 bp a clip matches most of
  # the decoy set. Worth 73 assignments over 47 cells that had none, 4% on 1186, with no wild-type
  # tube assignment added.
  #
  # The history, because this number has been wrong twice in opposite directions. 8 came from the
  # wild-type tubes being empty, which bounds cross-tube hopping rather than within-tube chance
  # matching - the wrong instrument. 10 came from decoy hits in coarse bands under the mismatch
  # budget, which is the right instrument read at the wrong resolution: in bands, and with one
  # mismatch allowed, 6-9 bp looks dirty. Raise or lower this only against a fresh per-base decoy
  # curve, and never without deciding exact_below in the same breath.
  # ONE OBJECT, NOT TWO CONSTANTS, and that is the guard rather than a tidiness preference.
  # min_len = 6 is only defensible with exact matching beneath exact_below: under the L %/% 10
  # mismatch budget a 6-base clip is allowed one mismatch and 21 of the 35 molecules at that length
  # match a decoy flank - a 60% error rate. Asking a reader not to separate the two numbers is not a
  # guard; making them one value is. Anything that wants to move the floor has to move the policy.
  assign_floor   = list(min_len = 6L, exact_below = 15L),
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
  # Decoys measure the assignment tier's error rate instead of arguing about it. They are flanks cut
  # from random positions carrying no provirus, passed through the identical positional test, so
  # under the null a decoy is exchangeable with a real flank and its hit rate IS epsilon:
  #   expected spurious assignments to real sites ~= n_real * decoy_hits / n_decoys
  # They also audit the positional filter itself. A decoy has no junction, so ANY decoy hit that
  # passes anchoring means the anchoring is mis-implemented. That is a stronger statement than the
  # wild-type tubes can make: those bound cross-tube hopping and carry-over, not within-tube chance
  # matching, and conflating the two is a reasoning error this analysis has already made once.
  # They must share the FASTA with the real flanks - a separate run would not share the aligner's
  # view of the search space, which is the thing being measured.
  # 1000, not 100. With 100 the decoy curve has no power where it matters: observed hits track
  # expected-under-chance at every length (10 against 9.24 at 4 bp, 0 against 1.67 at 6 bp), so a
  # zero at the floor carries p = 0.19 and resolves nothing. At 1000 flanks chance predicts ~13 hits
  # at 6 bp, so a zero there is a measurement rather than an absence of evidence. The cost is the
  # matcher, which runs over every flank: ~9x more comparisons in the assignment step.
  n_decoy        = 1000L,
  decoy_seed     = 20261117L,
  decoy_prefix   = "DECOY",
  # The genome-mapped channel: reads STAR placed on the host genome with the vector soft-clipped,
  # the mirror image of the 3' channel. Discovery needs the clip to match the provirus terminus
  # exactly over at least this many bases, measured AT the clip boundary. 13 is set by specificity,
  # not taste: boundary-constrained on chr19's 5.34M soft-clipped reads, 13 returns 29 hits at one
  # position (IS03) and nothing else, and the yield saturates by 12 (32 hits), so there is nothing
  # below it to gain. Scale any background claim genome-wide before quoting it - one chromosome is
  # 2.3% of the genome and cannot see a rate this low.
  min_genomic    = 13L,
  # MANDATORY, not a refinement. The provirus's terminal A is absorbed into the host alignment
  # whenever it matches the host's first base, so the terminus overhangs the clip boundary by 1-2
  # bases. Without this tolerance the boundary test returns ZERO hits at a site known to be there,
  # and the failure presents as "the channel does not work" rather than as a bug.
  absorb_window  = 2L,
  # The provirus 3' terminus, as it appears at the clip boundary. Forward for a plus-strand site;
  # read_genomic_clips() reverse-complements it for the minus-strand case, which is the orientation
  # trap that made the first probe of this channel return zero at every site.
  terminus       = "TGTGGAAAATCTCTAGCA",
  # The plasmid's two LTRs are byte-identical over all 181 bases and the contig keeps only one copy,
  # so a molecule spanning the *5'* LTR's junction with psi matches the contig's 3' terminus
  # perfectly: it is placed at the junction with MAPQ 255 and clips into vector sequence,
  # structurally identical to a host junction. No MAPQ or coverage filter sees it, because the anchor
  # genuinely is a unique alignment. Only the clip's sequence gives it away.
  psi_prefix     = "GTGGCGCCCGAACAGG"
)

# Checked when this file is sourced, not left to a reader's discretion. The floor and the exact-match
# threshold are one policy (see assign_floor): a clip admitted below ~15 bases MUST be matched
# exactly, because the L %/% 10 budget allows it one mismatch and a 6-base clip with one mismatch
# matches 60% of a decoy set. If someone lowers the floor without the exact-match rule reaching it,
# this stops the run at source() rather than producing a quietly worse table.
local({
  fl <- IS_DEFAULTS$assign_floor
  if (!(fl$min_len >= 15L || fl$exact_below >= 15L)) {
    stop("assign_floor is incoherent: min_len = ", fl$min_len, " admits clips below 15 bases while ",
         "exact matching is only required below ", fl$exact_below, ". A short clip under the ",
         "mismatch budget matches most of a decoy set; see clip_floor_curve.csv.", call. = FALSE)
  }
  if (fl$exact_below < fl$min_len) {
    stop("assign_floor is incoherent: exact_below (", fl$exact_below, ") is below min_len (",
         fl$min_len, "), so no admitted clip is ever matched exactly.", call. = FALSE)
  }
})


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
    clip_len >= opt$assign_floor$min_len       ~ "assign",
    TRUE                             ~ "clip too short to place"
  ))
}

# ── The genome-mapped channel ─────────────────────────────────────────────────────────────────────

#' Junction reads STAR placed on the HOST genome, with the vector soft-clipped.
#'
#' The mirror image of read_junction_clips(). A read spanning an integration junction has a vector
#' part and a host part, and STAR keeps whichever alignment scores better. When the host part is the
#' longer one the read is reported against the genome, with the vector clipped - so the junction's
#' genomic coordinate is already resolved, at MAPQ 255, with 120-140 bases of host behind it. That is
#' an order of magnitude better placement than the 25-50 bp clips the 3' channel has to align.
#'
#' Those reads were invisible to this analysis because nothing looked for them: the 3' channel reads
#' one window of one contig, and these never touch the contig. Intersected against the cells the
#' mini-reference already labels, the overlap is low - 12 of 67 at IS01, 2 of 20 at IS03 - so the two
#' channels are close to independent captures of the same clone, which is what makes a
#' capture-recapture clone size estimable at all.
#'
#' Two uses, with different evidence requirements, exactly as discovery and assignment differ:
#'   discovery  - no prior site needed. The clip must match the terminus over >= min_genomic bases
#'                AT its boundary, which on chr19's 5.34M soft-clipped reads returns hits at one
#'                position and nothing else.
#'   assignment - a shorter match, but the junction coordinate must already be a called site.
#'
#' `region` restricts the scan; NULL walks the whole BAM, which is ~17 minutes a tube.
read_genomic_clips <- function(bam, tube, region = NULL, opt = IS_DEFAULTS, k = opt$min_genomic) {
  what <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq")
  tags <- c("CB", "UB")
  p <- if (is.null(region)) {
    Rsamtools::ScanBamParam(what = what, tag = tags)
  } else {
    Rsamtools::ScanBamParam(what = what, tag = tags, which = region)
  }
  # Whole-BAM scan: stream it. These BAMs hold ~765M alignments and their SEQ alone would be over
  # 100 GB in memory, so with no `which` the file is read in yields and each yield reduced to the
  # handful of junction reads it contains before the next is fetched. The recursion is one level
  # deep - each call below passes a materialised chunk, never NULL - and keeps the matching logic in
  # exactly one place rather than duplicating it for the two access patterns.
  if (is.null(region)) {
    # 5e5 rather than 1e6: eight of these run forked in parallel, one per tube, and each yield holds
    # that many 150-base SEQs plus R's overhead on them.
    bf <- Rsamtools::BamFile(bam, yieldSize = 5e5L)
    open(bf); on.exit(close(bf), add = TRUE)
    acc <- list()
    repeat {
      got <- Rsamtools::scanBam(bf, param = p)[[1]]
      if (!length(got$cigar)) break
      acc[[length(acc) + 1L]] <- read_genomic_clips_chunk(got, tube, opt, k)
    }
    return(dplyr::bind_rows(acc))
  }

  chunks <- Rsamtools::scanBam(Rsamtools::BamFile(bam), param = p)
  dplyr::bind_rows(lapply(chunks, read_genomic_clips_chunk, tube = tube, opt = opt, k = k))
}

#' The matching itself, over one scanBam chunk. Shared by the windowed and whole-BAM paths.
#'
#' Flattening scanBam's chunks is where this went wrong twice, so nothing is flattened here: the
#' conversion happens inside the chunk. unlist() on `rname` returns the FACTOR CODES rather than the
#' chromosome names, so every coordinate lands on a chromosome called "19"; and
#' do.call(c, lapply(...)) hands the DNAStringSet to c() as a NAMED argument - the name being the
#' range - which returns a list that as.character() then deparses into one 80,000-character string.
#' Both failures produce an object of plausible class and length, so neither errors and the function
#' simply returns nothing.
read_genomic_clips_chunk <- function(x, tube, opt = IS_DEFAULTS, k = opt$min_genomic) {
  if (!length(x$cigar)) return(empty_genomic_table(tube))
  b <- dplyr::tibble(chrom = as.character(x$rname), pos = as.integer(x$pos),
                     mapq = as.integer(x$mapq), cigar = x$cigar,
                     seq = as.character(x$seq),
                     barcode = if (is.null(x$tag$CB)) NA_character_ else x$tag$CB,
                     umi     = if (is.null(x$tag$UB)) NA_character_ else x$tag$UB)
  if (!nrow(b)) return(empty_genomic_table(tube))

  # Cheap prefilter before any per-base work: a junction read must carry the terminus somewhere and
  # must have a soft clip to carry it in. This is what makes a 765M-alignment scan affordable.
  term0  <- substring(opt$terminus, nchar(opt$terminus) - k + 1L)
  rterm0 <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(term0)))
  pre <- grepl("S", b$cigar, fixed = TRUE) &
    (grepl(term0, b$seq, fixed = TRUE) | grepl(rterm0, b$seq, fixed = TRUE))
  if (!any(pre)) return(empty_genomic_table(tube))
  b <- b[pre, , drop = FALSE]

  cig <- b$cigar
  # suppressWarnings because as.integer() is evaluated for every element, including the CIGARs with
  # no clip on that side, where sub() returns the CIGAR unchanged and coercion gives NA. ifelse()
  # discards those, but the warning is noise that would otherwise fire on every call.
  lead  <- suppressWarnings(
    ifelse(grepl("^[0-9]+S", cig), as.integer(sub("^([0-9]+)S.*", "\\1", cig)), 0L))
  trail <- suppressWarnings(
    ifelse(grepl("[0-9]+S$", cig), as.integer(sub(".*[^0-9]([0-9]+)S$", "\\1", cig)), 0L))
  seqs <- b$seq
  nseq <- nchar(seqs)

  # The terminus as it appears at the clip boundary. A plus-strand provirus puts host AFTER it, so
  # the read is [vector][host]: the clip leads and ends with the terminus. A minus-strand provirus
  # is stored reverse-complemented, so the clip trails and STARTS with the complement. Searching only
  # the forward form is why the first probe of this channel returned zero at every minus-strand site
  # while working at the plus-strand ones - a failure that looks like absence, not like a bug.
  term  <- substring(opt$terminus, nchar(opt$terminus) - k + 1L)
  rterm <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(term)))

  # ABSORPTION. The provirus's terminal base is absorbed into the host alignment whenever it matches
  # the host's first base, so the terminus overhangs the clip boundary by 0 to absorb_window bases.
  # Without this the test returns zero hits at a site known to be there. `absorbed` is also the
  # correction the junction coordinate needs, so it is carried out of here rather than recomputed.
  w <- opt$absorb_window
  hit <- rep(FALSE, length(cig)); absorbed <- rep(NA_integer_, length(cig)); side <- rep("", length(cig))
  for (a in 0:w) {
    # leading clip: the terminus ends `a` bases past the clip's end
    e  <- lead + a
    ok <- !hit & lead >= k & e <= nseq & (e - k + 1L) >= 1L &
      substr(seqs, e - k + 1L, e) == term
    hit[ok] <- TRUE; absorbed[ok] <- a; side[ok] <- "lead"
    # trailing clip: the complement starts `a` bases before the clip's start
    s  <- nseq - trail + 1L - a
    ok <- !hit & trail >= k & s >= 1L & (s + k - 1L) <= nseq &
      substr(seqs, s, s + k - 1L) == rterm
    hit[ok] <- TRUE; absorbed[ok] <- a; side[ok] <- "trail"
  }
  if (!any(hit)) return(empty_genomic_table(tube))

  i <- which(hit)
  span <- ref_span(cig[i])
  # The junction base is the host base abutting the provirus, which is the NEAR end of the mapped
  # segment - the end the clip is on - not the far one. A leading clip puts the vector before the
  # host, so the junction is the segment's START; a trailing clip puts it after, so the junction is
  # the segment's END. Taking the far end instead lands a whole host-length away, which presents as
  # junctions scattered 100-136 bp downstream of the real site rather than as an obvious error.
  #
  # Absorption then moves it back toward the vector, because the absorbed bases were counted as host
  # by the aligner but belong to the provirus: a leading clip's start is pushed left, so add; a
  # trailing clip's end is pushed right, so subtract. Both forms are checked against the hand
  # measurement - IS02 (lead) 53,079,914 + 2 = 53,079,916 and IS01 (trail) 83,681,312 - 1 =
  # 83,681,311, each the called site exactly.
  junction <- ifelse(side[i] == "lead",
                     b$pos[i] + absorbed[i],
                     b$pos[i] + span - 1L - absorbed[i])

  dplyr::tibble(
    chrom     = b$chrom[i],
    pos       = as.integer(junction),
    strand    = ifelse(side[i] == "lead", "+", "-"),
    mapq      = b$mapq[i],
    tube      = tube,
    barcode   = b$barcode[i],
    umi       = b$umi[i],
    clip_len  = ifelse(side[i] == "lead", lead[i], trail[i]) + absorbed[i],
    host_len  = as.integer(span),
    absorbed  = absorbed[i]
  ) |>
    dplyr::filter(!is.na(chrom), !is.na(barcode), !is.na(umi))
}

empty_genomic_table <- function(tube) {
  dplyr::tibble(chrom = character(), pos = integer(), strand = character(), mapq = integer(),
                tube = character(), barcode = character(), umi = character(),
                clip_len = integer(), host_len = integer(), absorbed = integer())
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
    clip_len >= opt$assign_floor$min_len            ~ "assign",
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
    clip_len = nchar(field(10)[keep]),
    # The clip's own bases. Carried so tsd_check() can test the duplicated span against OBSERVED
    # sequence rather than only against the reference, which cannot falsify a coordinate convention.
    seq      = field(10)[keep]
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
                tube = character(), barcode = character(), umi = character(), clip_len = integer(),
                seq = character())
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

  # The inclusion gate. Assignment recruits cells to a flank, so a flank for a site that is not
  # really there manufactures a clone - which is how JH584304.1 reached 14 cells from 3, every one of
  # them in a mCherry-negative tube.
  #
  # Gating on the TSD alone was considered and rejected: only four of nine sites have any 5' reads at
  # all, so it would delete IS02 with its 145 cells. Gating on the bare floor is what already failed.
  # The gate is therefore TSD-verified OR comfortably clear of the floor, with sites sitting at
  # EXACTLY the floor admitted only if the 5' channel vouches for them. Certainty can come from
  # weight of evidence as well as from a second channel.
  if (!is.null(sites$tsd_pass) && nrow(sites)) {
    at_floor <- sites$cells <= opt$min_cells
    gated    <- at_floor & !isTRUE_vec(sites$tsd_pass)
    if (any(gated)) {
      message("    held out of the mini-reference (at the cell floor, no 5' confirmation): ",
              paste(sites$site_id[gated], collapse = ", "))
      seqs  <- seqs[!gated]
      sites <- sites[!gated, , drop = FALSE]
    }
  }
  if (length(seqs) == 0L) stop("every called site's flank is low-complexity; nothing to assign to")
  # name carries where the junction sits inside the flank, so assignment need not recompute it
  names(seqs) <- sprintf("%s|%s|%d", sites$site_id, sites$strand,
                         sites$pos - pmax(sites$pos - opt$flank, 1L) + 1L)

  decoys <- draw_decoy_flanks(sites, fa, opt)
  if (length(decoys)) {
    message("    plus ", length(decoys), " decoy flanks, to measure the assignment error rate")
    seqs <- c(seqs, decoys)
    # The decoys set the assignment floor, so the floor can move with nothing else changing if they
    # are redrawn differently. They are already seeded (opt$decoy_seed), but a seed is only a promise
    # until something checks it: the digest of the drawn set is recorded next to the flanks, and a
    # change in it is visible in the diff of a run rather than only in the numbers it produces.
    dg <- paste(substr(openssl::md5(paste(names(decoys), as.character(decoys), collapse = "|")),
                       1, 12), collapse = "")
    writeLines(c(paste0("n_decoy\t", length(decoys)),
                 paste0("decoy_seed\t", opt$decoy_seed),
                 paste0("flank\t", opt$flank),
                 paste0("digest\t", dg)),
               file.path(dirname(path), "decoy_provenance.tsv"))
    message("    decoy set digest ", dg, " (seed ", opt$decoy_seed, ")")
  }
  Biostrings::writeXStringSet(seqs, path)
  invisible(length(seqs))
}

#' Flanks from random positions with no provirus, named so they pass through the identical test.
#'
#' Matched to the real flanks on the three things that decide whether a clip can land: the same
#' chromosomes (sampled in proportion to length, so the scaffolds a real site sits on are in play),
#' the same width, and the same low-complexity filter. That is mappability matching by proxy rather
#' than from a mappability track - enough for the hit rate to be interpretable, and the one
#' approximation in this measurement worth stating out loud.
#'
#' The junction offset in the name is the centre, exactly as for a real flank, so the positional test
#' in assign_clips_exact() treats a decoy and a site identically. A decoy has no junction, so any hit
#' that survives is either chance matching at the real rate or a bug in the anchoring.
draw_decoy_flanks <- function(sites, fa, opt = IS_DEFAULTS) {
  if (opt$n_decoy <= 0L) return(Biostrings::DNAStringSet())
  idx  <- Rsamtools::scanFaIndex(fa)
  keep <- as.character(GenomicRanges::seqnames(idx)) %in% unique(sites$chrom)
  idx  <- idx[keep]
  if (!length(idx)) return(Biostrings::DNAStringSet())

  width <- 2L * opt$flank + 1L
  set.seed(opt$decoy_seed)                       # reproducible: the rate must not move between runs
  # Oversample, then filter. Random draws land in N-runs and homopolymers often enough that a single
  # pass would return fewer than asked for, and silently.
  want <- opt$n_decoy
  out  <- Biostrings::DNAStringSet()
  for (attempt in seq_len(10L)) {
    if (length(out) >= want) break
    n   <- 3L * (want - length(out))
    chr <- sample(seq_along(idx), n, replace = TRUE,
                  prob = as.numeric(GenomicRanges::width(idx)))
    len <- as.numeric(GenomicRanges::width(idx))[chr]
    ok  <- len > width + 1L
    if (!any(ok)) next
    chr <- chr[ok]; len <- len[ok]
    pos <- floor(stats::runif(length(chr), opt$flank + 1L, len - opt$flank))
    gr  <- GenomicRanges::GRanges(
      as.character(GenomicRanges::seqnames(idx))[chr],
      IRanges::IRanges(pos - opt$flank, pos + opt$flank))
    s   <- Rsamtools::scanFa(fa, gr)
    txt <- as.character(s)
    good <- !grepl("N", txt, fixed = TRUE) & max_homopolymer(txt) <= opt$max_homopolymer
    s <- s[good]
    if (length(s)) {
      names(s) <- sprintf("%s%03d|+|%d", opt$decoy_prefix,
                          length(out) + seq_along(s), opt$flank + 1L)
      out <- c(out, s)
    }
  }
  if (length(out) > want) out <- out[seq_len(want)]
  out
}

#' Split assignment hits into real and decoy, and turn the decoy rate into epsilon.
#'
#' Under the null a decoy flank is exchangeable with a real one, so the per-flank decoy hit rate is
#' the per-flank chance-matching rate. Scaled to the number of real flanks it is the expected number
#' of spurious assignments in the reported table. Reported per clip-length band because that is the
#' parameter the floor controls: if epsilon is flat across bands the floor is not the binding
#' constraint, and if it climbs at the short end the band where it climbs is where to stop.
decoy_error_rate <- function(mini, n_real, opt = IS_DEFAULTS) {
  if (nrow(mini) == 0L) return(dplyr::tibble())
  is_decoy <- startsWith(mini$site_id, opt$decoy_prefix)
  band <- function(x) cut(x, breaks = c(0, 9, 12, 15, 20, 30, Inf),
                          labels = c("<=9", "10-12", "13-15", "16-20", "21-30", ">30"),
                          right = TRUE)
  # Built by summarise() rather than pivot_wider(): when there are no decoy hits at all - the
  # outcome this measurement hopes for - the wide form has no TRUE column and every downstream
  # reference to it errors. The table must be able to report zero.
  dplyr::tibble(clip_band = band(mini$clip_len), is_decoy = is_decoy) |>
    dplyr::group_by(clip_band) |>
    dplyr::summarise(real = sum(!is_decoy), decoy = sum(is_decoy), .groups = "drop") |>
    dplyr::mutate(
      per_decoy_flank = decoy / opt$n_decoy,
      expected_spurious = round(per_decoy_flank * n_real, 2),
      epsilon = ifelse(real > 0L, round(expected_spurious / real, 4), NA_real_))
}

#' Every clip against every called flank, at every mismatch count. No floor, no tolerance.
#'
#' The one matcher. It replaced the bwa pass against the mini-reference, and it unified the two
#' assignment paths that had grown apart: the genomic channel already tested the clip boundary by
#' exact comparison with a +/-2 absorption tolerance, while this tier asked bwa for a best hit, threw
#' the ranking away, and tested position anyway. Only one of those was ever measured.
#'
#' The aligner was never doing useful work here. MAPQ against a handful of short flanks is
#' near-meaningless, so it was already ignored; `-a` had to be passed to stop the ranking discarding
#' true hits; and `-T 12` silently voided the 8 bp floor the thresholds advertised - clips of 8-11 bp
#' had a placement rate of exactly 0.00, 169 reads written to assign.fa and dropped with no error,
#' the same failure as `-T 30` against the 25 bp discovery floor one tier up.
#'
#' What is left is the test that was doing the work all along: does the clip match the flank at the
#' junction base, forward or reverse-complemented. This function answers that for EVERY clip against
#' EVERY flank and records the mismatch count, rather than returning a verdict.
#'
#' Recording `nm` rather than a 0/1 match is deliberate and it is what makes the clip table a
#' superset of the assignment table instead of a parallel one. A binary "does it match" bakes a
#' tolerance into the record, and then no view can ask a different question of it or reproduce the
#' old answer. With `nm` and `budget` carried per row, assign_clips_exact() is a two-line filter,
#' an exact-match analysis is `nm == 0`, and a lower floor is `min_len`. One matcher, several views,
#' and the equality gate in the assign stage can prove they agree.
clip_site_matches <- function(clip_fa, flank_fa, opt = IS_DEFAULTS, min_len = 1L) {
  if (!file.exists(clip_fa) || !file.exists(flank_fa)) return(empty_match_table())
  clips  <- Biostrings::readDNAStringSet(clip_fa)
  flanks <- Biostrings::readDNAStringSet(flank_fa)
  if (!length(clips) || !length(flanks)) return(empty_match_table())

  prov <- do.call(rbind, strsplit(names(clips), "|", fixed = TRUE))
  len  <- Biostrings::width(clips)
  txt  <- as.character(clips)
  use  <- len >= min_len
  if (!any(use)) return(empty_match_table())
  prov <- prov[use, , drop = FALSE]; len <- len[use]; txt <- txt[use]

  fref <- do.call(rbind, strsplit(names(flanks), "|", fixed = TRUE))
  fseq <- as.character(flanks)
  fjun <- as.integer(fref[, 3])

  nmis <- function(a, b) {
    am <- do.call(rbind, strsplit(a, "", fixed = TRUE))
    bm <- do.call(rbind, strsplit(b, "", fixed = TRUE))
    rowSums(am != bm)
  }
  revcomp <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))

  # Decoys are only compared against SHORT clips, and that is where most of the cost went. Chance
  # matching at 30 bases is 4^-30: there is nothing to measure there, and comparing 1000 decoy flanks
  # against the 748 molecules of 7 bases and up was paying for it anyway. Capping decoy comparisons
  # at exact_below keeps every comparison that carries information and drops the rest.
  is_decoy_flank <- startsWith(fref[, 1], opt$decoy_prefix)
  decoy_max_len <- opt$assign_floor$exact_below

  out <- list()
  for (s in seq_along(fseq)) {
    j <- fjun[s]; fl <- fseq[s]; n <- nchar(fl)
    lens <- sort(unique(len))
    if (is_decoy_flank[s]) lens <- lens[lens <= decoy_max_len]
    for (L in lens) {
      k <- which(len == L)
      # Forward: the clip's first base IS the junction base, so it starts at j and runs right.
      if (j + L - 1L <= n) {
        win <- substr(fl, j, j + L - 1L)
        out[[length(out) + 1L]] <- dplyr::tibble(
          row = k, site = s, strand_hit = "+", nm = nmis(txt[k], rep(win, length(k))),
          flank_pos = j, span = L)
      }
      # Reverse: the clip is the reverse complement, so it ENDS at j and runs left.
      if (j - L + 1L >= 1L) {
        win <- revcomp(substr(fl, j - L + 1L, j))
        out[[length(out) + 1L]] <- dplyr::tibble(
          row = k, site = s, strand_hit = "-", nm = nmis(txt[k], rep(win, length(k))),
          flank_pos = j - L + 1L, span = L)
      }
    }
  }
  if (!length(out)) return(empty_match_table())
  h <- dplyr::bind_rows(out)

  dplyr::tibble(
    site_id   = fref[h$site, 1],
    pos       = fjun[h$site],
    strand    = fref[h$site, 2],
    mapq      = NA_integer_,          # no aligner, so no aligner's opinion to record
    tube      = prov[h$row, 1],
    barcode   = prov[h$row, 2],
    umi       = prov[h$row, 3],
    clip_len  = len[h$row],
    clip      = txt[h$row],
    nm        = h$nm,
    budget    = pmax(1L, len[h$row] %/% 10L),
    strand_hit = h$strand_hit,
    flank_pos = h$flank_pos,
    span      = h$span,
    junction  = fjun[h$site]
  )
}

#' The shape clip_site_matches() returns when there is nothing to return.
empty_match_table <- function() {
  dplyr::tibble(site_id = character(), pos = integer(), strand = character(),
                mapq = integer(), tube = character(), barcode = character(), umi = character(),
                clip_len = integer(), clip = character(), nm = integer(), budget = integer(),
                strand_hit = character(), flank_pos = integer(), span = integer(),
                junction = integer())
}

#' The rule that turns a match into an assignment. Defined once so nothing can test a different one.
#'
#' Two thresholds held as one policy, `assign_floor`: a clip must reach `min_len`, and beneath
#' `exact_below` it must match exactly rather than within the L %/% 10 budget. Keeping them in one
#' predicate is not tidiness - the equality gate exists to prove the clip table and the pipeline
#' agree, and it can only do that if both ask the same question. When this was two copies of a
#' filter, changing the floor in one place would have left the gate passing against the old rule.
assign_keep <- function(clip_len, nm, budget, opt = IS_DEFAULTS) {
  fl <- opt$assign_floor
  clip_len >= fl$min_len & nm <= ifelse(clip_len < fl$exact_below, 0L, budget)
}

#' Assign clips to sites: the floor-and-tolerance VIEW over clip_site_matches().
#'
#' This used to be the matcher. It is now a filter over one, which is the point: the clip table and
#' the assignment table cannot disagree about whether a clip matches a flank, because only one piece
#' of code decides. The two thresholds that make it an assignment rather than a record are applied
#' here and nowhere else - `min_assign` on the clip length, and the same `L %/% 10` mismatch budget
#' the bwa path used to apply through NM.
assign_clips_exact <- function(clip_fa, flank_fa, opt = IS_DEFAULTS) {
  m <- clip_site_matches(clip_fa, flank_fa, opt, min_len = opt$assign_floor$min_len)
  if (nrow(m) == 0L) return(empty_hit_table())
  m |>
    dplyr::filter(assign_keep(clip_len, nm, budget, opt)) |>
    dplyr::select(site_id, pos, strand, mapq, tube, barcode, umi, clip_len, flank_pos, span,
                  junction)
}

#' isTRUE() elementwise, treating NA as FALSE: a site with no 5' reads has not been confirmed.
isTRUE_vec <- function(x) !is.na(x) & x

#' Longest run of a single base in each sequence.
max_homopolymer <- function(x) {
  vapply(regmatches(x, gregexpr("(.)\\1*", x)),
         function(runs) max(nchar(runs)), numeric(1))
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

#' Replace the clustering's internal ids with IS01.. labels, pinned to loci by a registry.
#'
#' Done once, before either summary, so the site table and the per-cell table cannot disagree about
#' which site is which.
#'
#' The ordering is by size, which makes the tables readable and the labels UNSTABLE, and an earlier
#' version of this function learned that the hard way twice. First three sites tied at 3 cells and
#' reshuffled between runs; ordering by locus within a tie fixed that. Then the genome-mapped channel
#' joined discovery, chr9:53,079,916 overtook chr13:83,681,311 on cell count, and IS01 through IS04
#' all moved - silently repointing every ID in the report at a different locus. A rank is not an
#' identity. Nothing that depends on which rows exist, or on how many cells they have, can be one.
#'
#' So the labels are pinned instead. `registry` is a tracked CSV of site_id, chrom, pos: a locus
#' already in it keeps its number for good, a new locus takes the lowest number not yet spent, and
#' the file is rewritten with the union. It lives under metadata/ rather than output/ because
#' output/ is ignored by git - a registry that a clean checkout does not have is not a registry.
#' Matching is by position within one cluster window, so a junction call that moves by a base or two
#' between runs is still the same site.
#'
#' With no registry the old behaviour returns: numbered by rank, stable only as long as the ranks
#' are. That path exists for one-off probes, not for anything whose output gets written down.
label_sites <- function(hits, registry = NULL, opt = IS_DEFAULTS) {
  if (nrow(hits) == 0L) return(hits)
  # Ties are broken by locus, not left to grouping order.
  order <- hits |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(cells = dplyr::n_distinct(paste(tube, barcode)),
                     chrom = dplyr::first(chrom), pos = min(pos), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(cells), chrom, pos)

  if (is.null(registry)) {
    order$label <- sprintf("IS%02d", seq_len(nrow(order)))
  } else {
    reg <- if (file.exists(registry)) {
      readr::read_csv(registry, show_col_types = FALSE, progress = FALSE)
    } else {
      dplyr::tibble(site_id = character(), chrom = character(), pos = integer())
    }
    lab <- rep(NA_character_, nrow(order))
    for (i in seq_len(nrow(order))) {
      near <- which(reg$chrom == order$chrom[i] &
                      abs(as.integer(reg$pos) - order$pos[i]) <= opt$cluster_window)
      if (length(near)) lab[i] <- reg$site_id[near[which.min(abs(reg$pos[near] - order$pos[i]))]]
    }
    # Lowest number not already spent, so a retired site's number is never reused.
    spent <- unique(c(reg$site_id, lab[!is.na(lab)]))
    free  <- setdiff(sprintf("IS%02d", seq_len(nrow(order) + length(spent))), spent)
    lab[is.na(lab)] <- utils::head(free, sum(is.na(lab)))
    order$label <- lab

    new_rows <- dplyr::tibble(site_id = order$label, chrom = order$chrom, pos = order$pos) |>
      dplyr::filter(!site_id %in% reg$site_id)
    if (nrow(new_rows)) {
      dir.create(dirname(registry), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(dplyr::arrange(dplyr::bind_rows(reg, new_rows), site_id), registry)
      message("    registry: ", nrow(new_rows), " new site(s) pinned -> ",
              paste(new_rows$site_id, collapse = ", "))
    }
  }

  hits |>
    dplyr::left_join(dplyr::select(order, site_id, label), by = "site_id") |>
    dplyr::mutate(site_id = label) |>
    dplyr::select(-label)
}

#' Promote 5' junctions that have no called 3' site, by looking for the 3' partner at q -/+ 4.
#'
#' The mini-reference logic run backwards. A 5' junction at q predicts its 3' partner at exactly one
#' position, so searching there is a single-position test rather than a scan, and that is what makes
#' a 2-read cluster admissible evidence. chrX:100,691,117 was found this way: three 3' clips at
#' MAPQ 60 in two cells at 100,691,113, signed delta +4 matching the strand - a real site the pooled
#' `>= min_cells` floor had dropped at 2 cells in the 3' channel.
#'
#' Reports its own denominator, without which one match is unfalsifiable: `candidates` is how many
#' orphan 5' junctions were tested and `rescued` how many found a partner. Forty tested and one
#' rescued is compelling; forty and eight is a coincidence rate.
rescue_orphan_junctions <- function(five_hits, discovery_hits, sites, opt = IS_DEFAULTS,
                                    expect = 4L, tol = 2L) {
  if (nrow(five_hits) == 0L) return(dplyr::tibble())
  five <- cluster_hits(five_hits, opt) |>
    dplyr::group_by(site_id) |>
    dplyr::summarise(chrom = dplyr::first(chrom), pos_5p = min(pos),
                     strand_5p = dplyr::first(strand),
                     reads_5p = dplyr::n(),
                     cells_5p = dplyr::n_distinct(paste(tube, barcode)), .groups = "drop") |>
    dplyr::select(-site_id)

  # Orphans only: a 5' junction already within a TSD's reach of a called site is not a candidate.
  orphan <- vapply(seq_len(nrow(five)), function(i) {
    if (nrow(sites) == 0L) return(TRUE)
    !any(sites$chrom == five$chrom[i] & abs(sites$pos - five$pos_5p[i]) <= 50L)
  }, logical(1))
  cand <- five[orphan, , drop = FALSE]
  if (!nrow(cand)) return(dplyr::tibble())

  rows <- lapply(seq_len(nrow(cand)), function(i) {
    # The partner sits at q - 4 or q + 4 depending on orientation; both are tested, and which one
    # hits also tells you the strand. Tolerance is only for the absorbed-base case.
    want <- c(cand$pos_5p[i] - expect, cand$pos_5p[i] + expect)
    hit  <- discovery_hits[discovery_hits$chrom == cand$chrom[i] &
                             vapply(discovery_hits$pos, function(p)
                               any(abs(p - want) <= tol), logical(1)), , drop = FALSE]
    dplyr::tibble(
      chrom = cand$chrom[i], pos_5p = cand$pos_5p[i],
      reads_5p = cand$reads_5p[i], cells_5p = cand$cells_5p[i],
      pos_3p = if (nrow(hit)) as.integer(stats::median(hit$pos)) else NA_integer_,
      reads_3p = nrow(hit),
      cells_3p = if (nrow(hit)) dplyr::n_distinct(paste(hit$tube, hit$barcode)) else 0L,
      umis_3p  = if (nrow(hit)) dplyr::n_distinct(paste(hit$tube, hit$barcode, hit$umi)) else 0L,
      max_mapq = if (nrow(hit)) max(hit$mapq) else NA_integer_,
      delta    = if (nrow(hit)) as.integer(cand$pos_5p[i] - stats::median(hit$pos)) else NA_integer_)
  })
  out <- dplyr::bind_rows(rows)
  out$rescued <- !is.na(out$pos_3p) & abs(out$delta) == expect
  # The denominator travels with the result, as an attribute so it survives a write_csv of the rows.
  attr(out, "candidates") <- nrow(cand)
  attr(out, "rescued")    <- sum(out$rescued)
  out |> dplyr::arrange(dplyr::desc(rescued), dplyr::desc(cells_5p + cells_3p))
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
    # Does this site's read-through signal get counted as its HOST GENE's expression? Cell Ranger
    # tags these reads with the host gene and sets the xf bit that makes them a UMI count, so for an
    # affected site the carrying cells gain host-gene UMIs that are vector transcript - in exactly
    # the cells any cis-effect test would compare. Measured: 20.9% of carriers' Abhd17c UMIs at IS04,
    # 23.5% of Gm35330 at IS10, 12.1% of Stau2 at IS11, 2.3% of Hells at IS03.
    #
    # It needs BOTH conditions, and the second is what keeps it off most sites. The read has to land
    # inside the gene, and it has to be in the gene's own orientation: an antisense read is tagged
    # but not counted. Checked against the data at six sites and the rule predicted all six - the
    # three same-sense intragenic ones carry 2-24% junction UMIs, the three antisense ones 0.0%,
    # 0.4% and 0.0%.
    dplyr::mutate(cis_confound = grepl("^intragenic", context) & strand == gene_strand) |>
    dplyr::select(site_id, chrom, pos, strand, gene, gene_type, gene_strand, context, cis_confound,
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

#' Scrambled-terminus decoy for the GENOME-MAPPED channel.
#'
#' The 100 decoy flanks bound the assignment tier. Nothing bounded this one, and it needed bounding
#' more than the other, because four of the fifteen sites rest on this channel alone and IS10 - 99
#' cells, the fifth largest clone - is invisible without it. The channel's claim is that a read the
#' aligner placed at MAPQ 255 whose soft-clipped end matches the last `min_genomic` bases of the
#' provirus terminus, within the absorption window, is a junction read. This measures how often that
#' is true of a terminus that cannot be there.
#'
#' The decoy is a PERMUTATION OF THE TESTED K-MER, not of the whole terminus. Those differ and the
#' difference flatters the channel: the real 13-mer AAAATCTCTAGCA is A-rich (6 A), a permutation of
#' the full 18-mer leaves the tested 13 with 3, and A-rich k-mers are commoner in a mammalian genome,
#' so a composition-mismatched decoy understates the false-positive rate. Composition is matched
#' here and the shuffle is reseeded until it is neither the real k-mer nor a substring of the real
#' terminus.
#'
#' `region` exists because a whole-genome pass is ~105 minutes a tube. One chromosome is a real
#' measurement on a real slice; the argument is there so the bound can be tightened without
#' rewriting anything.
#'
#' Returns one row a tube: decoy hits, and the real hits on the same slice for the contrast. A
#' non-zero decoy count is not a rate to report but a threshold to raise.
#' `real_cache` skips the real pass and takes its counts from the genomic stage's cache instead.
#' Rescanning for them doubles the cost of the most expensive thing the pipeline does, to recompute
#' a number already on disk. Supply it for a whole-genome run; leave it NULL only when the cache does
#' not cover the slice being tested.
genomic_decoy_rate <- function(bam, tube, region = NULL, opt = IS_DEFAULTS,
                               seed = 20261118L, real_cache = NULL) {
  k      <- opt$min_genomic
  real_k <- substr(opt$terminus, nchar(opt$terminus) - k + 1L, nchar(opt$terminus))
  set.seed(seed)
  repeat {
    dec_k <- paste(sample(strsplit(real_k, "")[[1]]), collapse = "")
    if (dec_k != real_k && !grepl(dec_k, opt$terminus, fixed = TRUE)) break
  }
  opt_d <- opt
  # Only the last k bases are ever tested, so the prefix is padding and is deliberately not DNA.
  opt_d$terminus <- paste0(strrep("X", max(0L, nchar(opt$terminus) - k)), dec_k)

  real <- if (is.null(real_cache)) {
    read_genomic_clips(bam, tube, region = region, opt = opt)
  } else {
    r <- read_genomic_cache(real_cache)
    if (!is.null(region)) r <- r[r$chrom %in% as.character(GenomeInfoDb::seqnames(region)), ]
    r
  }
  decoy <- read_genomic_clips(bam, tube, region = region, opt = opt_d)
  dplyr::tibble(
    tube = tube, region = if (is.null(region)) "whole genome" else as.character(region),
    real_kmer = real_k, decoy_kmer = dec_k,
    real_reads = nrow(real), real_pos = dplyr::n_distinct(real$pos),
    decoy_reads = nrow(decoy), decoy_pos = dplyr::n_distinct(decoy$pos))
}

#' write_csv, and then check that what landed is what was sent.
#'
#' On 2026-11-18 the project's GPFS fileset reached its byte limit - 391T of 391T, ~23 GB free - and
#' writes into output/ began truncating. Not failing: TRUNCATING, keeping only the last buffer, with
#' no header, no warning and no error. A deliberate 5000-row write returned 88 lines starting
#' mid-table, and the same write to /tmp was complete. Four outputs were corrupted that way, one of
#' them cell_site_assignments.csv, which the report reads; nothing anywhere said so.
#'
#' A silent short write is worse than a failed one, because every consumer downstream - the report,
#' the notebook chunks, the next stage - reads whatever is there and produces plausible numbers from
#' it. So every table this pipeline writes goes through here, and a row count that does not match is
#' an error that stops the run.
#'
#' The check is a line count, which assumes no embedded newlines in any field. True here: the fields
#' are identifiers, integers and DNA. If that ever stops being true, compare on nrow() after a
#' read-back instead.
write_csv_checked <- function(x, path, ...) {
  readr::write_csv(x, path, ...)
  got <- length(readLines(path, warn = FALSE))
  want <- nrow(x) + 1L
  if (got != want) {
    stop("short write: ", path, " has ", got, " lines, expected ", want,
         ". The filesystem accepted a truncated file without erroring - check free space on ",
         "the fileset (df -h) before trusting anything else this run wrote.", call. = FALSE)
  }
  record_in_manifest(path, nrow(x))
  invisible(path)
}

#' Record a verified write in the output directory's manifest, so a report can prove the files behind
#' its numbers are the pipeline's own output.
#'
#' Why this and not targets' own file hashing. The report declares these CSVs as inputs, so targets
#' hashes them - but it hashes them when the report renders, which cannot tell a pipeline output from
#' one edited by hand, half-regenerated, or left over from a run with different parameters. The
#' manifest can, because the entry is written by the same call that wrote the data, from the frame
#' that was just line-counted, and it carries the commit the pipeline ran from. The provenance banner
#' compares each declared input against it: an md5 that is not the recorded one means the report is
#' reading something this pipeline did not produce, which is the one failure that otherwise looks
#' exactly like a normal render.
#'
#' A re-run stage replaces its own entries and leaves the others, since the five stages write
#' different files and are deliberately runnable one at a time.
#'
#' The commit is the full 40-character SHA, not the short one, because a short hash is often a valid
#' number: `8e55518` is 8 x 10^55518, so anything that guesses column types reads it as `Inf`. Read
#' this file with every column as character.
record_in_manifest <- function(path, rows) {
  man <- file.path(dirname(path), "manifest.tsv")
  entry <- data.frame(
    file    = basename(path),
    rows    = as.character(rows),
    bytes   = as.character(file.info(path)$size),
    md5     = unname(tools::md5sum(path)),
    written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    commit  = tryCatch(system2("git", c("rev-parse", "HEAD"),
                               stdout = TRUE, stderr = FALSE)[1],
                       error = function(e) NA_character_),
    stringsAsFactors = FALSE
  )
  old <- if (file.exists(man)) {
    o <- utils::read.delim(man, colClasses = "character", na.strings = character())
    o[o$file != entry$file, , drop = FALSE]
  } else entry[0, , drop = FALSE]
  out <- rbind(old, entry)
  utils::write.table(out[order(out$file), ], man, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(man)
}

#' Read a cached genomic-channel table with its schema fixed, never guessed.
#'
#' The genomic stage caches one CSV a tube so the expensive scan can be restarted without redoing
#' what is already done. Read back with guessed column types that cache is a trap: the two wild-type
#' tubes legitimately return ZERO rows, readr types every column in a header-only file as character,
#' and binding those to a tube with real rows fails with "Can't combine `pos` <character> and `pos`
#' <double>". The failure lands at the START of the rerun, after the first run has already been paid
#' for, and it is caused by the two tubes that are working correctly.
#'
#' So the schema is declared. It matches empty_genomic_table(), and a cache that does not fit it is
#' an error worth seeing rather than something to coerce around.
read_genomic_cache <- function(path) {
  readr::read_csv(
    path, progress = FALSE, show_col_types = FALSE,
    col_types = readr::cols(
      chrom = readr::col_character(), pos = readr::col_integer(),
      strand = readr::col_character(), mapq = readr::col_integer(),
      tube = readr::col_character(), barcode = readr::col_character(),
      umi = readr::col_character(), clip_len = readr::col_integer(),
      host_len = readr::col_integer(), absorbed = readr::col_integer()))
}

#' Which channel calls each site on its own, and which sites need both.
#'
#' The headline counts - 9 sites from the 3' channel, 13 from the genome-mapped one, 15 pooled - look
#' like they do not add up, and the reason they do is worth carrying in the table rather than in a
#' sentence. The 3' channel's sites are a subset of the genomic channel's, so pooling should give 13.
#' It gives 15 because two sites clear the cell floor only when the two channels' reads are SUMMED:
#' neither channel reaches `min_cells` there alone.
#'
#' That is legitimate evidence and it is also the weakest evidence in the set, so those two sites are
#' flagged `pooled_only`. They should be excluded from any repertoire count, diversity statistic or
#' denominator until a second line of evidence arrives - a TSD, a 5' read, or the spliced channel.
#' Reporting "15 sites" without that distinction overstates what was found.
#'
#' Computed by rerunning the same caller on each channel's reads alone, so it cannot drift from the
#' pooled call: same function, same thresholds, one subset at a time.
site_provenance <- function(discovery_hits, sites, opt = IS_DEFAULTS) {
  if (nrow(sites) == 0L) return(dplyr::tibble())
  loci <- function(h) {
    if (!nrow(h)) return(dplyr::tibble(chrom = character(), pos = integer()))
    call_sites(h, opt) |>
      dplyr::group_by(site_id) |>
      dplyr::summarise(chrom = dplyr::first(chrom), pos = min(pos), .groups = "drop") |>
      dplyr::select(chrom, pos)
  }
  near <- function(tbl, chrom, pos) {
    if (!nrow(tbl)) return(FALSE)
    any(tbl$chrom == chrom & abs(tbl$pos - pos) <= opt$cluster_window)
  }
  ch <- if ("channel" %in% names(discovery_hits)) discovery_hits$channel else "three_prime"
  three <- loci(discovery_hits[ch == "three_prime", , drop = FALSE])
  gen   <- loci(discovery_hits[ch == "genomic", , drop = FALSE])

  out <- vapply(seq_len(nrow(sites)), function(i) {
    a <- near(three, sites$chrom[i], sites$pos[i])
    b <- near(gen,   sites$chrom[i], sites$pos[i])
    if (a && b) "both channels" else if (a) "3' channel only" else if (b) "genomic channel only"
    else "pooled only"
  }, character(1))
  dplyr::tibble(site_id = sites$site_id, called_by = out,
                pooled_only = out == "pooled only")
}

#' Attach the genome-mapped channel's reads to the called sites, as cell evidence.
#'
#' The genomic channel already earned its place in DISCOVERY - it calls sites on its own. Until this
#' function existed it then dropped out: the cell table and the site counts were rebuilt from the
#' mini-reference assignments alone, so a site's clone size counted only the cells that happened to
#' carry a soft-clipped tail read against the vector contig and discarded every cell whose evidence
#' was a host-mapped read. The two are the same observation seen from opposite ends, and for IS01
#' the omission halved the count (308 cells at discovery, 154 after assignment).
#'
#' Matching is by position, to the nearest called site within one cluster window - the same window
#' call_sites() used to merge these reads into those sites in the first place, so this reproduces the
#' calling assignment rather than inventing a second one. Reads further than that from every site are
#' dropped: they are the positions that did not reach the cell floor, and they are not evidence for
#' any site that did.
#'
#' clip_len carries the HOST alignment length for these rows, not the vector clip. In this channel
#' the host is the anchor - 120-140 bases of it at MAPQ 255 - so host_len is the quantity that plays
#' the part clip_len plays in the 3' channel, and max_clip_len in the cell table stays interpretable
#' as "the longest host anchor behind this cell's call".
genomic_site_hits <- function(gen, site_table, opt = IS_DEFAULTS) {
  if (nrow(gen) == 0L || nrow(site_table) == 0L) {
    return(dplyr::tibble(site_id = character(), tube = character(), barcode = character(),
                         umi = character(), clip_len = integer(), channel = character(),
                         tier = character()))
  }
  n   <- nrow(gen)
  sid <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    near <- which(site_table$chrom == gen$chrom[i] &
                    abs(site_table$pos - gen$pos[i]) <= opt$cluster_window)
    if (!length(near)) next
    sid[i] <- site_table$site_id[near[which.min(abs(site_table$pos[near] - gen$pos[i]))]]
  }
  dplyr::tibble(
    site_id  = sid,
    tube     = gen$tube,
    barcode  = gen$barcode,
    umi      = gen$umi,
    clip_len = as.integer(gen$host_len),
    channel  = "genomic",
    tier     = "genomic"
  ) |>
    dplyr::filter(!is.na(site_id), !is.na(barcode), !is.na(umi))
}

#' Collapse the matcher's output to ONE ROW PER MOLECULE, with no floor applied.
#'
#' The unit is the molecule - tube + barcode + UMI - because that is the unit of capture. Reads of
#' one UMI are duplicates of one fragment, and counting them as separate evidence is the mistake that
#' put 150 cells in both channels when 83 belong there.
#'
#' One row per molecule, with channel membership as FLAGS. This shape is deliberate and it is a fix
#' rather than a preference. The table used to carry one row per molecule PER CHANNEL, and the cell
#' table carried `umis_three` and `umis_genomic` as separate counts - and as long as two counts sit
#' side by side, something will add them. Something did: the depth variable in the sharing analysis
#' was `umis_three + umis_genomic`, which double-counts the 83 molecules present in both channels and
#' pushes single-molecule cells into the two-molecule stratum, where a cell cannot show two sites by
#' construction. The same error was fixed for lambda and not propagated. With `in_three`/`in_genomic`
#' as logicals and one row per molecule, `dplyr::n()` is the molecule count and the wrong sum is not
#' expressible.
#'
#' `clip_len` and `clip` come from the molecule's LONGEST three-prime read, which is its most
#' informative one; `host_len` carries the genomic channel's host anchor separately, because the
#' two are different quantities and were previously merged into one column. `match_sites` comes
#' from the representative read alone - a union across a molecule's reads let a 5-base read's
#' matches through a floor applied to its 30-base sibling.
#'
#' Decoy flanks are matched too and are kept, labelled, because a decoy match at a given clip length
#' is the error rate and dropping it here would throw away the measurement.
molecule_clip_table <- function(raw, gen_site_hits = NULL, sites = NULL, opt = IS_DEFAULTS) {
  is_decoy <- function(x) startsWith(x, opt$decoy_prefix)
  three <- if (nrow(raw)) {
    rep <- raw |>
      dplyr::group_by(tube, barcode, umi) |>
      dplyr::mutate(reads = dplyr::n_distinct(clip_len), longest = max(clip_len)) |>
      dplyr::filter(clip_len == longest) |>
      dplyr::ungroup()
    rep |>
      dplyr::group_by(tube, barcode, umi) |>
      dplyr::summarise(
        n_reads_three = dplyr::first(reads),
        clip_len  = dplyr::first(clip_len),
        clip      = dplyr::first(clip),
        best_nm   = min(nm),
        n_match_three = dplyr::n_distinct(site_id[nm <= budget & !is_decoy(site_id)]),
        n_decoy_match = dplyr::n_distinct(site_id[nm <= budget & is_decoy(site_id)]),
        n_match_exact = dplyr::n_distinct(site_id[nm == 0L & !is_decoy(site_id)]),
        n_decoy_exact = dplyr::n_distinct(site_id[nm == 0L & is_decoy(site_id)]),
        sites_three = paste(sort(unique(site_id[nm <= budget & !is_decoy(site_id)])),
                            collapse = ";"),
        .groups = "drop")
  } else {
    dplyr::tibble(tube = character(), barcode = character(), umi = character())
  }
  gen <- if (!is.null(gen_site_hits) && nrow(gen_site_hits)) {
    gen_site_hits |>
      dplyr::group_by(tube, barcode, umi) |>
      dplyr::summarise(n_reads_genomic = dplyr::n(),
                       host_len = max(clip_len),
                       sites_genomic = paste(sort(unique(site_id)), collapse = ";"),
                       n_match_genomic = dplyr::n_distinct(site_id),
                       .groups = "drop")
  } else {
    dplyr::tibble(tube = character(), barcode = character(), umi = character())
  }

  dplyr::full_join(three, gen, by = c("tube", "barcode", "umi")) |>
    dplyr::mutate(
      in_three   = !is.na(n_reads_three),
      in_genomic = !is.na(n_reads_genomic),
      channels   = dplyr::case_when(in_three & in_genomic ~ "both",
                                    in_three ~ "three_prime", TRUE ~ "genomic"),
      dplyr::across(c(n_reads_three, n_reads_genomic, n_match_three, n_match_genomic,
                      n_decoy_match, n_match_exact, n_decoy_exact),
                    ~ tidyr::replace_na(.x, 0L)),
      # The union across channels, which is what "does this molecule name a site" means.
      match_sites = vapply(seq_along(in_three), function(i) {
        s <- c(strsplit(tidyr::replace_na(sites_three[i], ""), ";")[[1]],
               strsplit(tidyr::replace_na(sites_genomic[i], ""), ";")[[1]])
        paste(sort(unique(s[nzchar(s)])), collapse = ";")
      }, character(1)),
      n_match = pmax(n_match_three, n_match_genomic)) |>
    dplyr::select(tube, barcode, umi, channels, in_three, in_genomic,
                  n_reads_three, n_reads_genomic, clip_len, host_len, clip, best_nm,
                  n_match, n_match_three, n_match_genomic, n_match_exact,
                  n_decoy_match, n_decoy_exact, match_sites) |>
    dplyr::arrange(tube, barcode, umi)
}

#' Does the matcher, filtered at the current floor, reproduce the pipeline's own output?
#'
#' The gate the whole refactor exists to pass. It runs on the READ-level match table, not the
#' molecule view, and that distinction is the result of getting it wrong once: the pipeline's
#' assignment tier is a per-read filter, so only a per-read view can be expected to reproduce it to
#' the row. The molecule table is a convenience built on top for steps that need molecules; it is
#' not the thing under test.
#'
#' Filtering at `min_assign` with the same mismatch budget, restricted to called sites and unioned
#' with the genomic channel's cells, must give back the site count, the assignment count and the cell
#' count the pipeline produced by its usual route. Anything else means the table and the pipeline
#' disagree about the same quantity, and then nothing taken from the table - the floor curve, lambda,
#' the posteriors - means anything.
#'
#' It is deliberately an equality and not a tolerance. A clip table that is 1% off is not 99% right,
#' it is wrong in a way nobody will track down later.
clip_table_gate <- function(matches, site_table, cell_table, gen_site_hits = NULL,
                            called = NULL, opt = IS_DEFAULTS) {
  view <- matches |>
    dplyr::filter(assign_keep(clip_len, nm, budget, opt)) |>
    dplyr::semi_join(dplyr::select(site_table, site_id), by = "site_id") |>
    dplyr::select(site_id, tube, barcode)
  if (!is.null(gen_site_hits) && nrow(gen_site_hits)) {
    view <- dplyr::bind_rows(view, dplyr::select(gen_site_hits, site_id, tube, barcode))
  }
  # The cell table is restricted to called cells, so the view must be too, or the gate compares a
  # filtered count against an unfiltered one and fails for a reason that has nothing to do with the
  # matcher.
  if (!is.null(called) && nrow(called)) {
    view <- dplyr::semi_join(view, called, by = c("tube", "barcode"))
  }

  dplyr::tibble(
    quantity = c("sites", "cell-site assignments", "distinct cells"),
    pipeline = c(nrow(site_table), nrow(cell_table),
                 dplyr::n_distinct(paste(cell_table$tube, cell_table$barcode))),
    clip_table = c(nrow(site_table),
                   dplyr::n_distinct(paste(view$site_id, view$tube, view$barcode)),
                   dplyr::n_distinct(paste(view$tube, view$barcode)))
  ) |>
    dplyr::mutate(agrees = pipeline == clip_table)
}

#' Per-cell vector-transcript UMIs, for every called cell that carries any.
#'
#' The per-tube share answers whether the sort gate is pure. This answers the per-cell questions
#' that depend on it: whether a cell with no junction molecule is nonetheless transduced, and
#' whether the per-clone expression differences that this report reads off vector UMIs are
#' confounded by anything. Same three reporter genes, same matrix, one row per cell rather than one
#' per tube.
#'
#' Only cells with at least one vector UMI are returned - about 26,000 of 161,000 - because a table
#' of zeros is the same information at six times the size. A cell absent from this table has none.
vector_expression_cells <- function(bam_dir, tubes, genes = c("mCherry", "DRE", "WPRE-LTR")) {
  dplyr::bind_rows(lapply(tubes, function(tube) {
    d <- file.path(bam_dir, paste0("cellranger_", tube), "sample_filtered_feature_bc_matrix")
    f_feat <- file.path(d, "features.tsv.gz"); f_mtx <- file.path(d, "matrix.mtx.gz")
    f_bc   <- file.path(d, "barcodes.tsv.gz")
    if (!all(file.exists(f_feat, f_mtx, f_bc))) return(dplyr::tibble())
    feat <- readLines(gzfile(f_feat))
    want <- which(sub("\t.*", "", feat) %in% genes)
    bcs  <- readLines(gzfile(f_bc))
    hdr_raw <- readLines(gzfile(f_mtx), n = 10L)
    n_comment <- sum(startsWith(hdr_raw, "%"))
    trip <- data.table::fread(f_mtx, skip = n_comment + 1L, showProgress = FALSE,
                              col.names = c("row", "col", "val"))
    hit <- trip[trip$row %in% want, ]
    if (!nrow(hit)) return(dplyr::tibble())
    dplyr::tibble(tube = tube, barcode = bcs[hit$col],
                  gene = sub("\t.*", "", feat[hit$row]), umis = hit$val) |>
      dplyr::group_by(tube, barcode) |>
      dplyr::summarise(vector_umis = sum(umis),
                       genes_detected = dplyr::n_distinct(gene), .groups = "drop")
  }))
}

#' One row per site, every flag that bears on whether to trust it, in one place.
#'
#' The flags exist already and are scattered across five outputs, so deciding which sites are
#' defensible means re-deriving a judgement that has already been made three times in this analysis.
#' Collected here, downstream work can filter to the defensible set without reconstructing the
#' reasoning: `defensible` is the single column to filter on, and the columns beside it say why.
#'
#' `defensible` is FALSE for a site that clears the cell floor only on the two channels' reads pooled
#' (`pooled_only`) or whose flank cannot support quantification (`flank_gated`). It is deliberately
#' not a quality score - there is no threshold being applied to a continuum here, just two named
#' disqualifications and the evidence for each.
site_evidence <- function(site_table, tsd = NULL, detection = NULL, five = NULL,
                          flank_fa = NULL, opt = IS_DEFAULTS) {
  if (nrow(site_table) == 0L) return(dplyr::tibble())
  in_mini <- if (!is.null(flank_fa) && file.exists(flank_fa)) {
    nm <- grep("^>", readLines(flank_fa), value = TRUE)
    sub("\\|.*", "", sub("^>", "", nm))
  } else NA_character_

  ev <- dplyr::tibble(
    site_id      = site_table$site_id,
    chrom        = site_table$chrom,
    pos          = site_table$pos,
    strand       = site_table$strand,
    gene         = site_table$gene,
    cells        = site_table$cells,
    umis         = site_table$umis,
    called_by    = site_table$called_by,
    pooled_only  = site_table$pooled_only,
    cis_confound = site_table$cis_confound,
    flank_gated  = if (length(in_mini) == 1L && is.na(in_mini[1])) NA
                   else !site_table$site_id %in% in_mini)

  if (!is.null(tsd) && nrow(tsd)) {
    st <- tsd |>
      dplyr::filter(!is.na(site_id)) |>
      dplyr::group_by(site_id) |>
      dplyr::summarise(tsd_status = dplyr::first(tsd),
                       tsd_seq = dplyr::first(tsd_seq),
                       n_5prime_reads = sum(reads_5p), .groups = "drop")
    ev <- dplyr::left_join(ev, st, by = "site_id")
  }
  if (!is.null(detection) && nrow(detection)) {
    ev <- dplyr::left_join(
      ev, dplyr::transmute(detection, site_id, molecules_per_cell = mean_molecules,
                           lambda_estimable = !is.na(lambda), detection, clone_size),
      by = "site_id")
  }
  ev |>
    dplyr::mutate(
      tsd_status = tidyr::replace_na(tsd_status, "no 5' reads"),
      n_5prime_reads = tidyr::replace_na(n_5prime_reads, 0L),
      defensible = !pooled_only & !tidyr::replace_na(flank_gated, FALSE)) |>
    dplyr::arrange(dplyr::desc(cells))
}

#' Cells carrying vector transcript, per tube, from the filtered feature-barcode matrix.
#'
#' The discriminator the completeness figure needed, and it is independent of junctions entirely:
#' vector expression is measured from reads on the reporter contig, which the reference carries as
#' three genes - mCherry, DRE and WPRE-LTR. A sorted mCherry-positive cell with no vector transcript
#' at all is not a cell this assay can ever place, whether because the sort admitted it or because
#' its provirus is silent.
#'
#' It settles a question the junction data cannot. 93% of sorted mCherry+ cells belong to no called
#' site, and that is either a repertoire far larger than fifteen clones or an impure gate. The answer
#' is the gate: 65.2% of Rx_4 and 43.2% of Rx_6 called cells carry any vector UMI, 53.6% pooled - not
#' the ~93% the large-repertoire reading requires. The wild-type tube is the control and it works: 4
#' cells of 26,822, so the measurement has essentially no false-positive rate.
#'
#' Read the matrix directly rather than through Seurat: only three rows of it are needed, and the
#' triplet file can be streamed.
vector_expression <- function(bam_dir, tubes, genes = c("mCherry", "DRE", "WPRE-LTR")) {
  dplyr::bind_rows(lapply(tubes, function(tube) {
    d <- file.path(bam_dir, paste0("cellranger_", tube), "sample_filtered_feature_bc_matrix")
    f_feat <- file.path(d, "features.tsv.gz"); f_mtx <- file.path(d, "matrix.mtx.gz")
    if (!file.exists(f_feat) || !file.exists(f_mtx)) {
      return(dplyr::tibble(tube = tube, cells = NA_integer_, vector_positive = NA_integer_,
                           vector_umis = NA_integer_))
    }
    feat <- readLines(gzfile(f_feat))
    want <- which(sub("\t.*", "", feat) %in% genes)
    # fread, not readLines + strsplit and not scan(): the triplet body of a filtered matrix runs to
    # tens of millions of lines. Splitting it in R takes minutes and scan() in chunks was no better;
    # fread reads it in seconds. Only three of ~33,000 feature rows are wanted, but the triplets are
    # ordered by cell, so there is nothing to skip - the file has to be read and filtered.
    head_lines <- readLines(gzfile(f_feat), n = 1L)  # force the connection, then read the matrix
    hdr_raw <- readLines(gzfile(f_mtx), n = 10L)
    n_comment <- sum(startsWith(hdr_raw, "%"))
    hdr <- scan(text = hdr_raw[n_comment + 1L], what = integer(), quiet = TRUE)
    ncell <- if (length(hdr) >= 2L) hdr[2] else NA_integer_
    trip <- data.table::fread(f_mtx, skip = n_comment + 1L, showProgress = FALSE,
                              col.names = c("row", "col", "val"))
    hit <- trip[trip$row %in% want, ]
    pos <- hit$col
    umis <- sum(hit$val)

    dplyr::tibble(tube = tube, cells = ncell,
                  vector_positive = dplyr::n_distinct(pos), vector_umis = umis)
  })) |>
    dplyr::mutate(vector_share = vector_positive / cells)
}

#' Repertoire completeness: what share of the transduced population the called sites account for.
#'
#' The denominator this needed turned out to be in hand all along - the mCherry-positive cells
#' SEQUENCED, not the FACS record. Two numbers come out, and they should always be reported together
#' because they bracket the answer from opposite sides.
#'
#' TWO DENOMINATORS, because the obvious one is wrong. Sorted mCherry+ cells were used at first, and
#' vector_expression() shows that gate is impure: only 65.2% of Rx_4 and 43.2% of Rx_6 called cells
#' carry any vector transcript, 53.6% pooled. A cell with no vector transcript cannot produce a
#' junction read - whether it was never transduced or its provirus is wholly silent - so it does not
#' belong in the denominator of a junction-based completeness measure. Both are reported: the sorted
#' count is the upper bound on the population, the vector-positive count the defensible one.
#'
#' `floor_share` is labelled cells over the denominator. It assumes nothing: even at perfect
#' detection the called sites cannot account for less than this. `corrected_share` divides each
#' site's detected count by its own per-cell detection rate first, and so inherits lambda's
#' uncertainty.
#'
#' THE IMPUTATION HAS A DIRECTION, and it is the one that flatters the result. lambda is estimable
#' only at the six sites with the highest molecules per cell - a sample selected for high detection -
#' so imputing their median to the other nine OVERSTATES detection, understates those clone sizes and
#' overstates completeness. `corrected_share` is therefore an UPPER bound and `floor_share` the only
#' model-free figure. Never quote the corrected number without the floor beside it.
#'
#' What it cannot do, stated because this was got wrong once: the unassigned-molecule census does NOT
#' bound the hidden population. A clone that yields no junction molecule leaves none to count, so a
#' molecule census has no power over exactly the population it was being used to reassure about.
repertoire_completeness <- function(detection_tbl, mcherry_cells,
                                   vector_positive_cells = NA_integer_) {
  if (nrow(detection_tbl) == 0L) return(dplyr::tibble())
  med <- stats::median(detection_tbl$detection, na.rm = TRUE)
  d <- dplyr::mutate(detection_tbl,
                     detection_used = ifelse(is.na(detection), med, detection),
                     imputed = is.na(detection),
                     size_used = cells / detection_used)
  corrected <- sum(d$size_used)
  dplyr::tibble(
    denominator = c("sorted mCherry+ cells", "cells carrying vector transcript"),
    cells = c(mcherry_cells, vector_positive_cells),
    sites = nrow(d),
    detected_cells = sum(d$cells),
    floor_share = sum(d$cells) / c(mcherry_cells, vector_positive_cells),
    median_detection_used = med,
    sites_with_imputed_detection = sum(d$imputed),
    corrected_cells = round(corrected),
    corrected_share = corrected / c(mcherry_cells, vector_positive_cells),
    unaccounted_share = 1 - corrected / c(mcherry_cells, vector_positive_cells))
}

#' The floor curve: what a clip of each length matches, real sites and decoy flanks alike.
#'
#' The measurement that sets `min_assign`, computed rather than quoted. For every clip length it
#' reports how many molecules matched exactly one called site, how many matched more than one, and
#' how many matched any of the 100 decoy flanks - and a decoy has no junction, so a decoy match at a
#' given length IS the chance rate at that length. Decoys can be cut at will, so this has a sample
#' size the mono-replicate control cannot reach.
#'
#' Reported under both rules, because the floor and the mismatch tolerance are one decision: `_exact`
#' columns require nm == 0, the others allow the L %/% 10 budget. Under the budget a 6-base clip is
#' allowed a mismatch and decoy hits persist to 8 bp; under exact matching they stop at 6.
#'
#' AND WITH THE POWER OF EACH ZERO, because "no decoy hits" is not one statement. The decoy dimension
#' has no sample-size limit - 100 flanks can be 1,000 - but the binding constraint is MOLECULES, and
#' there are only 35 at 6 bp. `expected_decoy` is how many of them chance alone should have put on a
#' decoy: a specific window of length L matches by chance at 4^-L, each flank is tested in both
#' orientations, so the per-molecule rate is 1 - (1 - 2*4^-L)^n_decoy. At 6 bp that is 4.8% and
#' expects 1.7 hits across 35 molecules, so observing zero carries p = 0.18 under chance alone -
#' unremarkable. `p_zero_under_chance` states that directly.
#'
#' `upper_bound` is the number the floor should actually be read from: with zero observed, the 95%
#' Poisson ceiling is 3 events, which bounds the per-molecule error against the real flank set at
#' 3 * n_real / (molecules * n_decoy * orientations). At 6 bp that is ~0.5-1%, sitting AT the 1%
#' tolerance rather than below it. A zero with no power behind it should not be mistaken for a clean
#' one, which is why these columns exist rather than the bare count.
clip_floor_curve <- function(clip_tbl, sites = NULL, opt = IS_DEFAULTS, orientations = 2L) {
  if (nrow(clip_tbl) == 0L) return(dplyr::tibble())
  n_real <- if (is.null(sites)) NA_integer_ else nrow(sites)
  clip_tbl |>
    dplyr::filter(in_three) |>
    dplyr::group_by(clip_len) |>
    dplyr::summarise(
      molecules       = dplyr::n(),
      unique_budget   = sum(n_match == 1L),
      ambiguous_budget = sum(n_match > 1L),
      decoy_budget    = sum(n_decoy_match > 0L),
      unique_exact    = sum(n_match_exact == 1L),
      ambiguous_exact = sum(n_match_exact > 1L),
      decoy_exact     = sum(n_decoy_exact > 0L),
      decoy_flanks_exact = mean(n_decoy_exact),
      .groups = "drop") |>
    dplyr::mutate(
      at_or_above_floor = clip_len >= opt$assign_floor$min_len,
      exact_required    = clip_len < opt$assign_floor$exact_below,
      n_decoy           = opt$n_decoy,
      # chance of a specific length-L window matching, per flank, per orientation
      per_flank_chance  = 4^-clip_len,
      per_molecule_chance = 1 - (1 - orientations * per_flank_chance)^opt$n_decoy,
      expected_decoy    = molecules * per_molecule_chance,
      p_zero_under_chance = exp(-expected_decoy),
      # The standing test is this ratio, not the raw count. A count tells you how many chance
      # matches a given decoy set produced; the ratio tells you whether the matcher still behaves as
      # theory says, which is what a regression check needs and what the 100-decoy version could not
      # measure at all. 15 observed against 13.5 predicted at 6 bp is the form to read.
      obs_over_expected = ifelse(expected_decoy > 0, decoy_exact / expected_decoy, NA_real_),
      # ...and where that ratio can be read at all. With an expectation below a few events a single
      # chance hit sends the ratio to 12 or 30 - which is what 9 bp and 10 bp do here - and a
      # standing check that alarms on Poisson noise is the monotonicity test's mistake again. Compare
      # obs/exp only where this is TRUE; everywhere else read expected_decoy and decoy_exact directly.
      ratio_interpretable = expected_decoy >= 5,
      # 95% Poisson ceiling on zero observed, expressed against the REAL flank set
      upper_bound = ifelse(decoy_exact == 0L & molecules > 0L & !is.na(n_real),
                           3 * n_real / (molecules * opt$n_decoy * orientations), NA_real_),
      # THE NUMBER THAT LICENSES THE FLOOR. The decoy zero does not: observed hits track
      # expected-under-chance at every length (10 against 9.24 at 4 bp, 1 against 1.60 at 5 bp,
      # 0 against 1.67 at 6 bp), so with 100 flanks the decoys measure 4^-L and resolve no
      # transition. What they do establish, where they have power, is that real genomic sequence
      # matches at about the analytic rate - and that rate can then be applied to the flank set that
      # actually matters, the real one.
      chance_vs_real = 1 - (1 - orientations * per_flank_chance)^n_real,
      expected_spurious = molecules * chance_vs_real) |>
    dplyr::arrange(clip_len)
}

#' Where the junction molecules that match no called site actually go.
#'
#' The honest answer to "if most of the transduced population belongs to no called site, why are
#' there no junction reads from it". There are, and this says what becomes of them: too short to
#' place, long enough but matching nothing, and of those, how many bwa could not place uniquely
#' versus how many sit at a real position that failed the cell floor. The second group is the
#' repertoire genuinely missed; the first is the repetitive-flank blind spot.
#'
#' `discovery` is read_bwa_hits() output, used to separate the two.
unassigned_molecules <- function(clip_tbl, discovery, cell_table, opt = IS_DEFAULTS) {
  if (nrow(clip_tbl) == 0L) return(dplyr::tibble())
  assigned <- unique(paste(cell_table$tube, cell_table$barcode))
  un <- clip_tbl |> dplyr::filter(n_match == 0L)
  placed <- if (nrow(discovery)) {
    dplyr::distinct(discovery, tube, barcode, umi) |> dplyr::mutate(placed_uniquely = TRUE)
  } else {
    dplyr::tibble(tube = character(), barcode = character(), umi = character(),
                  placed_uniquely = logical())
  }
  un |>
    dplyr::left_join(placed, by = c("tube", "barcode", "umi")) |>
    dplyr::mutate(
      placed_uniquely = tidyr::replace_na(placed_uniquely, FALSE),
      class = dplyr::case_when(
        !in_three                       ~ "genome-mapped, no called site within the window",
        clip_len < opt$assign_floor$min_len       ~ "below the assignment floor",
        clip_len < opt$min_discovery    ~ "above the floor, below the discovery length",
        placed_uniquely                 ~ "long and placed uniquely: a locus below the cell floor",
        TRUE                            ~ "long but unplaceable: MAPQ < 30 or unmapped")) |>
    dplyr::group_by(class) |>
    dplyr::summarise(molecules = dplyr::n(),
                     cells = dplyr::n_distinct(paste(tube, barcode)),
                     cells_with_no_assignment =
                       dplyr::n_distinct(paste(tube, barcode)[!paste(tube, barcode) %in% assigned]),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(molecules))
}

#' Ambient: junction-bearing barcodes that Cell Ranger did not call as cells.
#'
#' Per-cell against per-empty-droplet, which is the comparison that makes the number mean something.
#' It needs the barcodes the called-cell filter removes, which is why that filter is applied to the
#' cell table and not to the clip and genomic tables.
ambient_rate <- function(clip_tbl, called, raw_counts) {
  if (nrow(clip_tbl) == 0L) return(dplyr::tibble())
  clip_tbl |>
    dplyr::distinct(tube, barcode) |>
    dplyr::left_join(dplyr::mutate(called, is_cell = TRUE), by = c("tube", "barcode")) |>
    dplyr::mutate(is_cell = tidyr::replace_na(is_cell, FALSE)) |>
    dplyr::group_by(tube) |>
    dplyr::summarise(junction_barcodes = dplyr::n(),
                     in_called_cells = sum(is_cell),
                     not_cells = sum(!is_cell),
                     ambient_share = sum(!is_cell) / dplyr::n(), .groups = "drop") |>
    dplyr::left_join(raw_counts, by = "tube") |>
    dplyr::mutate(rate_per_cell = in_called_cells / cells,
                  rate_per_empty = not_cells / pmax(raw_barcodes - cells, 1L),
                  # NA, not a guard constant. Dividing by pmax(rate, 1e-12) turned Rx_7's zero
                  # non-cell barcodes into a printed ratio of 816,888,086, which is not a large
                  # number - it is the absence of a denominator. A tube with no ambient barcodes has
                  # no ratio to report, and the same goes for one with no labelled cells.
                  ratio = ifelse(rate_per_empty > 0 & rate_per_cell > 0,
                                 rate_per_cell / rate_per_empty, NA_real_)) |>
    dplyr::arrange(dplyr::desc(junction_barcodes))
}

#' Cell Ranger's called barcodes, per tube, read from the filtered matrix.
#'
#' The BAM carries reads from barcodes Cell Ranger error-corrected to a CB but did NOT call as cells,
#' and this pipeline takes CB tags straight from the BAM. So 33 of what the cell table called cells
#' were not cells: 18 in Rx_4, 12 in Rx_6, one each in Rx_2, Rx_3 and Rx_5, none of them in any
#' filtered barcode list and therefore in none of the Seurat objects either. 2.8% of 1186.
#'
#' They are not noise to be discarded quietly - they are the ambient population, the same droplets
#' that give the 2.6% ambient fraction in Rx_4. So the list is read here and the FILTER is applied
#' where cells are counted, leaving the raw clip and genomic tables untouched: the ambient
#' measurement needs the non-cells, and the clone sizes must not include them.
#'
#' Written to work/called_barcodes.csv.gz by the extract stage so every count downstream is
#' reproducible from the run's own outputs, without reaching back into data/.
read_called_barcodes <- function(bam_dir, tubes) {
  dplyr::bind_rows(lapply(tubes, function(tube) {
    f <- file.path(bam_dir, paste0("cellranger_", tube),
                   "sample_filtered_feature_bc_matrix", "barcodes.tsv.gz")
    if (!file.exists(f)) {
      warning("no filtered barcode list for ", tube, "; its cells cannot be checked", call. = FALSE)
      return(dplyr::tibble(tube = character(), barcode = character()))
    }
    dplyr::tibble(tube = tube, barcode = readLines(gzfile(f)))
  }))
}

#' Called cells and total droplets per tube, for the ambient denominator.
#'
#' The ambient rate is only interpretable per droplet: a junction in one of 2.5 million empty
#' droplets and a junction in one of 20,000 cells are different events. The raw barcode count is the
#' denominator for the first, so it is captured with the run rather than recomputed from data/ later.
count_barcodes <- function(bam_dir, tubes) {
  dplyr::bind_rows(lapply(tubes, function(tube) {
    base <- file.path(bam_dir, paste0("cellranger_", tube))
    f_called <- file.path(base, "sample_filtered_feature_bc_matrix", "barcodes.tsv.gz")
    f_raw    <- file.path(base, "sample_raw_feature_bc_matrix", "barcodes.tsv.gz")
    dplyr::tibble(
      tube = tube,
      cells = if (file.exists(f_called)) length(readLines(gzfile(f_called))) else NA_integer_,
      raw_barcodes = if (file.exists(f_raw)) length(readLines(gzfile(f_raw))) else NA_integer_)
  }))
}

#' Junction molecules per cell, counted once each, with no floor.
#'
#' The quantity the detection model needs, and neither of the obvious ways of getting it is right.
#' Summing umis_three + umis_genomic over a cell's rows in the cell table double-counts the 83
#' molecules that appear in both channels, and it also drops every molecule too short to assign - but
#' a clip too short to name a site still proves a junction molecule was captured in that cell, which
#' is exactly what lambda is estimated from. Counting distinct UMIs in the clip table fixes both at
#' once, and the two errors happen to point in opposite directions: the mean goes from 1.271 to
#' 1.236, so the inflated one was closer to the truth by accident.
#'
#' `unattributed` is the molecules that match no called site. They are real captures at an uncalled
#' or unplaceable position, and they belong in the capture depth even though they name no clone.
cell_junction_molecules <- function(clip_tbl) {
  # One row per molecule in the clip table, so n() IS the molecule count. It used to be
  # n_distinct(umi) over a table with one row per molecule per channel, which was correct but only
  # by vigilance; the table's shape now carries the invariant instead.
  clip_tbl |>
    dplyr::group_by(tube, barcode) |>
    dplyr::summarise(molecules = dplyr::n(),
                     unattributed = sum(n_match == 0L),
                     in_both_channels = sum(in_three & in_genomic),
                     .groups = "drop")
}

#' Per-cell detection rate and detection-corrected clone size, by zero-truncated Poisson.
#'
#' The estimate the clone-size section needs. If the cells of a clone each yield Poisson(lambda)
#' DETECTED junction molecules, the cells with zero are exactly the ones missing from the table, and
#' the mean of the zero-truncated distribution, lambda / (1 - exp(-lambda)), identifies lambda from
#' the cells in hand. Detection is then 1 - exp(-lambda) and the clone size detected / detection. No
#' second sample is needed, which is what makes this usable where the channel overlap is not: the two
#' channels are a partition of one capture rather than two draws, so n1*n2/m has no referent.
#'
#' Fitted PER SITE, because a global lambda would push the measured seven-fold span in per-clone
#' vector expression into the clone sizes as bias.
#'
#' It is sharp only where the truncated mean sits clearly above 1. At exactly 1.00 the fit is
#' degenerate - lambda -> 0, clone size -> infinity - so sites at or below `min_mean` return NA
#' rather than a number. That is the estimator declining, not a gap to be filled: for a 2-cell site
#' the correction would be meaningless against its own uncertainty anyway. The bootstrap interval is
#' over cells and is wide by construction at these depths; report it, because a point estimate here
#' invites more confidence than the data supports.
#'
#' `obs_1`/`exp_1` and friends are the fit check: if the observed molecules-per-cell distribution is
#' overdispersed relative to the fitted Poisson, the Poisson understates the zero class and the clone
#' size is an underestimate. It fits well here - IS02 gives 240 cells at one molecule against 232
#' expected and 76 at two against 85, IS01 273 against 269 and 37 against 43 - so there is no
#' negative-binomial refit to report. That was not expected: the overdispersion the plan anticipated
#' was inferred from 150 both-channel cells against a predicted 110, and that excess turned out to be
#' the same molecule counted in two channels.
#'
#' Read the fit from obs against exp and NOT from var_mean_ratio. A zero-truncated Poisson with small
#' lambda is underdispersed by construction - the ratio is 0.41 at IS02 and 0.17 at IS01 - so
#' comparing it to 1 as if these were untruncated counts suggests a problem that is not there. The
#' column is kept because it is cheap and because its value is worth knowing; it is not the test.
detection_model <- function(cell_table, clip_tbl, opt = IS_DEFAULTS, min_mean = 1.10,
                            boots = 400L, seed = 20261118L) {
  if (nrow(cell_table) == 0L || nrow(clip_tbl) == 0L) return(dplyr::tibble())
  per <- cell_junction_molecules(clip_tbl)
  d   <- cell_table |>
    dplyr::distinct(site_id, tube, barcode) |>
    dplyr::left_join(per, by = c("tube", "barcode")) |>
    dplyr::filter(!is.na(molecules))

  ztp <- function(k) {
    if (length(k) < 5L || mean(k) <= min_mean) return(NA_real_)
    tryCatch(stats::uniroot(function(l) l / (1 - exp(-l)) - mean(k), c(1e-6, 50))$root,
             error = function(e) NA_real_)
  }
  set.seed(seed)
  d |>
    dplyr::group_by(site_id) |>
    dplyr::group_modify(function(x, key) {
      k <- x$molecules
      lam <- ztp(k)
      bs <- if (is.na(lam)) rep(NA_real_, boots) else
        vapply(seq_len(boots), function(i) ztp(sample(k, replace = TRUE)), numeric(1))
      det <- 1 - exp(-lam)
      # observed against the fitted zero-truncated Poisson, for k = 1, 2, 3+
      ex <- if (is.na(lam)) rep(NA_real_, 3) else {
        p <- stats::dpois(1:2, lam) / (1 - exp(-lam))
        c(p, 1 - sum(p)) * length(k)
      }
      dplyr::tibble(
        cells = length(k), mean_molecules = mean(k), unattributed = sum(x$unattributed),
        lambda = lam, detection = det, clone_size = round(length(k) / det),
        lambda_lo = stats::quantile(bs, 0.025, na.rm = TRUE, names = FALSE),
        lambda_hi = stats::quantile(bs, 0.975, na.rm = TRUE, names = FALSE),
        clone_lo = round(length(k) / (1 - exp(-stats::quantile(bs, 0.975, na.rm = TRUE,
                                                               names = FALSE)))),
        clone_hi = round(length(k) / (1 - exp(-stats::quantile(bs, 0.025, na.rm = TRUE,
                                                               names = FALSE)))),
        obs_1 = sum(k == 1L), exp_1 = ex[1], obs_2 = sum(k == 2L), exp_2 = ex[2],
        obs_3plus = sum(k >= 3L), exp_3plus = ex[3],
        var_mean_ratio = if (length(k) > 1L) stats::var(k) / mean(k) else NA_real_)
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(cells))
}

#' Per-site channel overlap: the detection rate it does support, and the clone size it does not.
#'
#' The tempting argument: the two channels read the same junctions through different reads, so a cell
#' carrying a clone can be captured by either, the overlap between them estimates the capture
#' probability, and Chapman's (n1+1)(n2+1)/(m+1) - 1 turns that into the number of cells carrying the
#' clone INCLUDING those neither channel saw. Run on this data it returns 650-950 cells for the large
#' sites against 130-342 detected, i.e. it claims the assay sees a third to a half of each clone.
#'
#' The columns after it are there to stop that being reported as a fact. Lincoln-Petersen needs the
#' two captures to be independent given the cell. Here they are not even close: which channel a
#' junction read falls in is decided by where its fragment straddles the junction - more vector than
#' host and it is read against the contig, more host than vector and it is read against the genome -
#' and every read of one UMI is a duplicate of one fragment. So the draw happens once per MOLECULE,
#' and a cell can appear in both channels only if it carries at least two junction molecules. Between
#' 73% and 100% of cells carry exactly one. m is therefore small for a reason that has nothing to do
#' with undetected cells, and a small m is exactly what inflates N_hat.
#'
#' m_exp_molecule is that model's prediction with NO hidden population at all: for a cell with k
#' molecules and a per-molecule probability p of landing in the 3' channel, P(both) = 1 - p^k -
#' (1-p)^k, summed over the cells actually observed. Across the eleven sites with at least ten cells
#' it predicts 110 overlapping cells against 150 observed - the same calculation done on reads
#' instead of molecules predicts 430, which is the whole reason this is counted in molecules.
#'
#' An earlier version of this comment reported a 40-cell excess of observed overlap over that
#' prediction, z = 6, and reached for residual correlation in capture to explain it. The excess was
#' an artefact of counting. umis_three + umis_genomic double-counts a molecule seen in both channels,
#' and 73 of 1434 molecules are: of the 150 cells that looked like two captures, 72 share a UMI
#' between the channels and 52 hold one molecule in total. m_distinct requires a three-prime-only
#' molecule AND a genomic-only one, and it reads 83 against 150. The prediction, 110, sits between
#' them - so the overlap is consistent with the per-molecule model and the excess is withdrawn.
#'
#' The estimator suffers more than the model does. Any lambda fitted from the channel overlap
#' inherits the inflated m and comes out 1.17x to 2.22x above the molecule-count lambda at every
#' site where both are defined, in exactly the direction double counting predicts. That is why this
#' function reports lambda from the molecule counts and not from the overlap: the molecule counts
#' do not require the partition to hold.
#'
#' The underlying objection is a partition, not a bias: Lincoln-Petersen needs two SAMPLES of one
#' population, and these two channels are a partition of one sample. Every junction molecule is in
#' exactly one of them, assigned by where its fragment fell. There is no second draw to recapture
#' anything, so n1*n2/m is not a mis-estimate of the clone - it has no referent in this design at
#' all, and no amount of bias correction gives it one. The molecule model above is what that
#' partition predicts quantitatively; the partition is why the quantity is meaningless.
#'
#' Hence: chapman is written to the table as an audit trail, because an argument this tempting should
#' be inspectable rather than deleted, and `lambda`/`detection`/`n_ztp` are what the section should
#' quote instead. Report detected as the floor and n_ztp where the fit is conditioned well enough to
#' mean anything.
#'
#' n1, n2 and m are counted in MOLECULES for the reason above. 73 of 1434 molecules do appear in both
#' channels, so a molecule's channel is not strictly fixed - fragmentation happens after
#' amplification - which makes umis_three + umis_genomic a slight over-count of k, and so pushes
#' m_exp_molecule up rather than down: the shortfall against the observed overlap is not an artefact
#' of that.
channel_capture_recapture <- function(cell_table) {
  if (nrow(cell_table) == 0L) return(dplyr::tibble())
  cell_table |>
    dplyr::group_by(site_id) |>
    dplyr::group_modify(function(d, key) {
      k <- d$umis_three + d$umis_genomic
      p <- if (sum(k) > 0L) sum(d$umis_three) / sum(k) else NA_real_
      n1 <- sum(d$umis_three > 0L)
      n2 <- sum(d$umis_genomic > 0L)
      m  <- sum(d$umis_three > 0L & d$umis_genomic > 0L)
      # m counted on DISTINCT molecules: the cell must hold a three-prime-only molecule AND a
      # genomic-only one. m_both counts a single shared molecule as two captures, which it is not.
      sh <- if ("umis_shared" %in% names(d)) d$umis_shared else rep(0L, nrow(d))
      m_d <- sum(d$umis_three > sh & d$umis_genomic > sh)
      # Zero-truncated Poisson on the per-cell molecule counts. This is the estimate the section
      # actually needs, and it comes from a different place than the overlap: if a clone's cells each
      # carry Poisson(lambda) DETECTED junction molecules, the cells with zero are exactly the ones
      # missing from the table, and the truncated mean k/(1 - exp(-k)) identifies lambda without
      # them. detection = 1 - exp(-lambda) is then a per-cell detection rate and detected/detection
      # a clone size.
      #
      # It is sharp only where the mean is comfortably above 1. At a mean of 1.00 the fit is
      # degenerate - lambda -> 0, detection -> 0, clone size -> infinity - so it returns NA there
      # rather than a number, and at 1.04 (IS03) it is formally defined and still worthless. IS05 is
      # the one well-conditioned case in this data at 2.65 molecules a cell, and it returns 63
      # against 57 detected, which is the kind of answer the estimator gives when it can see.
      lam <- if (length(k) >= 5L && mean(k) > 1.10) {
        tryCatch(stats::uniroot(function(l) l / (1 - exp(-l)) - mean(k), c(1e-6, 50))$root,
                 error = function(e) NA_real_)
      } else NA_real_
      dplyr::tibble(
        detected = nrow(d), umis_per_cell = mean(k), frac_one_umi = mean(k == 1L),
        lambda = lam, detection = 1 - exp(-lam),
        n_ztp = round(nrow(d) / (1 - exp(-lam))),
        n_three = n1, n_genomic = n2, m_both = m, m_distinct = m_d,
        cells_sharing_a_umi = sum(sh > 0L),
        lincoln_petersen = if (m > 0L) n1 * n2 / m else NA_real_,
        chapman = (n1 + 1) * (n2 + 1) / (m + 1) - 1,
        # Chapman's variance, Seber 1970.
        chapman_se = sqrt((n1 + 1) * (n2 + 1) * (n1 - m) * (n2 - m) / ((m + 1)^2 * (m + 2))),
        m_exp_molecule = if (is.na(p)) NA_real_ else sum(1 - p^k - (1 - p)^k),
        m_ratio = ifelse(m_exp_molecule > 0, m / m_exp_molecule, NA_real_))
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(detected))
}

#' Per-base coverage of the host flanks either side of each junction, by channel.
#'
#' The two channels approach a provirus from opposite ends: 3' clips are the host sequence that
#' FOLLOWS it and 5' clips the host that PRECEDES it. Plotted on a common axis anchored at the 3'
#' junction, they should therefore pile up on opposite sides of it and meet - overlapping by the
#' 5 bp target-site duplication. The TSD stops being a number in a table and becomes visible.
#'
#' Offsets are signed distances from the 3' junction in genomic coordinates, so 0 is the junction
#' base itself. Mini-reference hits are converted out of flank coordinates using the junction offset
#' carried in the reference name; 5' hits are already genomic.
junction_coverage <- function(mini_hits, five_sam, sites, opt = IS_DEFAULTS, window = 120L) {
  span_of <- function(cigar) ref_span(cigar)

  read_sam <- function(path) {
    if (!file.exists(path) || file.size(path) == 0L) return(NULL)
    lines <- readLines(path, warn = FALSE)
    lines <- lines[!startsWith(lines, "@")]
    if (!length(lines)) return(NULL)
    f <- strsplit(lines, "\t", fixed = TRUE)
    fl <- as.integer(vapply(f, `[`, "", 2L))
    keep <- bitwAnd(fl, 0x100) == 0L & bitwAnd(fl, 0x800) == 0L &
      vapply(f, `[`, "", 3L) != "*"
    if (!any(keep)) return(NULL)
    f <- f[keep]
    dplyr::tibble(rname = vapply(f, `[`, "", 3L),
                  pos   = as.integer(vapply(f, `[`, "", 4L)),
                  mapq  = as.integer(vapply(f, `[`, "", 5L)),
                  cigar = vapply(f, `[`, "", 6L))
  }

  expand <- function(site_id, start, len) {
    n <- rep(site_id, len)
    dplyr::tibble(site_id = n, offset = unlist(Map(function(a, l) a:(a + l - 1L), start, len)))
  }

  out <- list()

  # Takes the already-filtered hit table rather than a SAM path: only matches that passed the
  # positional test count, and drawing the unfiltered set spread IS01 over -49..+119, wider than a
  # 125 bp read can reach. Decoy flanks are not sites, so they are dropped here too.
  m <- mini_hits[mini_hits$site_id %in% sites$site_id, , drop = FALSE]
  if (nrow(m)) {
    out$three <- expand(m$site_id, m$flank_pos - m$junction, m$span) |>
      dplyr::mutate(channel = "3' junction (host after the provirus)")
  }

  f5 <- read_sam(five_sam)
  if (!is.null(f5)) {
    f5 <- dplyr::filter(f5, mapq >= opt$min_mapq)
    hit <- lapply(seq_len(nrow(f5)), function(i) {
      k <- which(sites$chrom == f5$rname[i] & abs(sites$pos - f5$pos[i]) <= 2000L)
      if (!length(k)) return(NULL)
      k <- k[which.min(abs(sites$pos[k] - f5$pos[i]))]
      len <- as.integer(span_of(f5$cigar[i]))
      if (len <= 0L) return(NULL)
      dplyr::tibble(site_id = sites$site_id[k],
                    offset  = (f5$pos[i] - sites$pos[k]):(f5$pos[i] - sites$pos[k] + len - 1L))
    })
    hit <- dplyr::bind_rows(hit)
    if (nrow(hit)) out$five <- dplyr::mutate(hit, channel = "5' junction (host before the provirus)")
  }

  if (!length(out)) return(dplyr::tibble())
  dplyr::bind_rows(out) |>
    dplyr::filter(abs(offset) <= window) |>
    dplyr::count(site_id, channel, offset, name = "depth")
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
#' `genome_fa`, when given, turns the coordinate test into a sequence test. Signed delta is reported
#' rather than absolute, because its SIGN is a free check on the strand call: a plus-strand provirus
#' puts its 5' junction above its 3' junction and a minus-strand one below, so sign(delta) must equal
#' the strand. All four sites with both channels already agree, and a disagreement would mean either
#' the strand call or one junction is wrong without saying which - which is why the sequence test
#' exists alongside it.
tsd_check <- function(sites, five_hits, opt = IS_DEFAULTS, search = 2000L, expect = 4L,
                      genome_fa = NULL) {
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
  strand  <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    near <- which(sites$chrom == five$chrom[i] & abs(sites$pos - five$pos_5p[i]) <= search)
    if (!length(near)) next
    k <- near[which.min(abs(sites$pos[near] - five$pos_5p[i]))]
    site_id[i] <- sites$site_id[k]
    pos_3p[i]  <- as.integer(sites$pos[k])
    strand[i]  <- as.character(sites$strand[k])
    # SIGNED, not absolute. The sign is the strand check; taking the modulus throws it away.
    delta[i]   <- as.integer(five$pos_5p[i] - sites$pos[k])
  }

  out <- five |>
    dplyr::mutate(site_id = site_id, pos_3p = pos_3p, delta = delta, strand = strand) |>
    # A TSD is 4 bp, so only a junction within a few bases of a site is a claim about that site.
    # Anything further is a distinct feature that merely happens to be nearby, and must not be
    # labelled FAIL against a site that passed on its own junction - chr19:38,932,385 sits 1.4 kb
    # from IS03 and would otherwise read as IS03 failing. That one is now explained: it is an
    # annotated Hells exon 3' end with a canonical GT after it, so those reads are a spliced
    # host-exon-to-provirus fusion rather than a second integration.
    dplyr::mutate(tsd = dplyr::case_when(
      is.na(delta)              ~ "no 3' site nearby",
      abs(delta) == expect      ~ "PASS",
      abs(delta) <= expect + 2L ~ "near",
      abs(delta) <= 50L         ~ "FAIL",
      TRUE                      ~ "distinct locus"
    )) |>
    # Free consistency check, reported rather than merely assumed: a plus-strand provirus has its 5'
    # junction above its 3' junction and a minus-strand one below.
    dplyr::mutate(strand_ok = dplyr::case_when(
      is.na(delta) | tsd == "distinct locus" ~ NA,
      TRUE ~ (delta > 0L) == (strand == "+")
    )) |>
    dplyr::mutate(site_id = ifelse(tsd == "distinct locus", NA_character_, site_id))

  # THE SEQUENCE TEST. A coordinate difference of 4 can be produced by an off-by-one in either
  # channel, by a 0- versus 1-based slip, or by taking the wrong end of an alignment - this analysis
  # has had all three. The duplicated bases themselves cannot: the two junctions must share one
  # 5-base span of reference, so deriving that span from pos_3p and from pos_5p independently and
  # requiring the two to be the same string is immune to the whole class. abs(delta) == 4 is kept as
  # the cross-check, and a disagreement between the two localises which convention is wrong.
  out$tsd_seq  <- NA_character_
  out$tsd_obs  <- NA_real_
  if (!is.null(genome_fa) && file.exists(genome_fa) && any(out$tsd %in% c("PASS", "near"))) {
    fa <- Rsamtools::FaFile(genome_fa)
    i  <- which(out$tsd %in% c("PASS", "near"))
    lo <- pmin(out$pos_3p[i], out$pos_5p[i])
    hi <- pmax(out$pos_3p[i], out$pos_5p[i])
    span <- as.character(Rsamtools::scanFa(
      fa, GenomicRanges::GRanges(out$chrom[i], IRanges::IRanges(lo, hi))))
    out$tsd_seq[i] <- span

    # The non-circular half. Deriving the span twice from the reference proves nothing - addressed
    # from either junction it is the same interval, so the comparison is a tautology (an earlier
    # version of this function made exactly that mistake and reported a pass). The reference is one
    # witness; the OBSERVED clip sequence is a second, independent one. A real TSD means the
    # duplicated bases are physically present in the 5' reads, so the span must appear in their
    # sequence - in either orientation, since SEQ is stored reference-forward.
    if (!is.null(five_hits$seq)) {
      out$tsd_obs[i] <- vapply(seq_along(i), function(j) {
        rs <- five_hits$seq[five_hits$chrom == out$chrom[i][j] &
                              abs(five_hits$pos - out$pos_5p[i][j]) <= opt$cluster_window]
        if (!length(rs)) return(NA_real_)
        rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(span[j])))
        mean(grepl(span[j], rs, fixed = TRUE) | grepl(rc, rs, fixed = TRUE))
      }, numeric(1))
    }
  }

  out |>
    dplyr::arrange(dplyr::desc(reads_5p)) |>
    dplyr::select(site_id, chrom, pos_3p, pos_5p, delta, tsd, tsd_seq, tsd_obs, strand, strand_ok,
                  reads_5p, cells_5p)
}

#' One row per cell per site: the table that joins onto the Seurat metadata.
#'
#' A cell at more than one site is flagged, not resolved: a genuine multi-copy integrant and a
#' doublet are indistinguishable here, so the choice belongs downstream.
#' `channel` and `tier` are carried through per cell, not collapsed. The capture-recapture argument
#' needs, per site, the 3' cell set, the genomic cell set and their overlap, and a single label per
#' cell destroys that irrecoverably - cheap to keep now, a re-run to recover later.
summarise_cells <- function(hits) {
  if (nrow(hits) == 0L) return(dplyr::tibble())
  # names(), not hits$channel: a tibble warns on an unknown column rather than returning NULL.
  if (!"channel" %in% names(hits)) hits$channel <- "three_prime"
  hits |>
    dplyr::group_by(site_id, tube, barcode) |>
    dplyr::summarise(reads = dplyr::n(),
                     umis = dplyr::n_distinct(umi),
                     max_clip_len = max(clip_len),
                     tiers = paste(sort(unique(tier)), collapse = "+"),
                     channels = paste(sort(unique(channel)), collapse = "+"),
                     reads_three = sum(channel == "three_prime"),
                     reads_genomic = sum(channel == "genomic"),
                     # Per-channel MOLECULE counts, not just read counts. The channel a read falls
                     # in is decided by where its fragment straddles the junction, and every read of
                     # one UMI is a duplicate of one fragment - so the lottery is drawn once per
                     # molecule, not once per read. Any question about how the two channels overlap
                     # has to be asked of these, and asking it of reads gives an answer that is
                     # wrong by a factor of three.
                     umis_three = dplyr::n_distinct(umi[channel == "three_prime"]),
                     umis_genomic = dplyr::n_distinct(umi[channel == "genomic"]),
                     # Molecules seen in BOTH channels. umis_three + umis_genomic double-counts
                     # these, and they are not a modelling nicety: 73 of 1434 molecules are shared,
                     # and they account for 67 of the 150 cells that looked like two captures. A
                     # molecule has one anchor length and so belongs to one channel, so a shared
                     # UMI means that assumption fails for that molecule - which is exactly the
                     # assumption any channel-overlap estimator rests on.
                     umis_shared = dplyr::n_distinct(
                       intersect(umi[channel == "three_prime"], umi[channel == "genomic"])),
                     # The distinct-molecule count, so no consumer ever has to add umis_three and
                     # umis_genomic. That sum double-counts the shared molecules and it has already
                     # been added by mistake once, in the sharing analysis, after being fixed for
                     # lambda and not propagated. Use this column.
                     molecules = dplyr::n_distinct(umi),
                     .groups = "drop") |>
    dplyr::group_by(tube, barcode) |>
    dplyr::mutate(n_sites_in_cell = dplyr::n_distinct(site_id),
                  ambiguous = n_sites_in_cell > 1L) |>
    dplyr::ungroup() |>
    dplyr::arrange(tube, barcode, site_id)
}
