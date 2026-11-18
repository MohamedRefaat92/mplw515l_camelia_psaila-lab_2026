#!/usr/bin/env Rscript
# Driver for the integration-site analysis. The work is in src/R/integration-sites.R; this file only
# decides which tubes to read and where the tables go.
#
# Four stages, because bwa has to run between the first and the third, and because the genomic scan
# is expensive enough that it must be restartable without redoing anything else:
#
#   extract <out> <bam_dir> <tube>...            clips out of the BAMs -> discovery.fa, assign.fa,
#                                                five_prime.fa, and the tables the report plots
#     bwa  -> discovery.sam, five_prime.sam      against the genome
#   genomic <out> <bam_dir> <tube>...            the host-mapped channel, one streaming pass per BAM
#   genomic-decoy <out> <bam_dir> [chrom] [fa]   its error bound: the same scan, terminus shuffled
#   sites   <out> <genes.gtf.gz> <genome.fa>     call the sites, write the mini-reference
#   assign  <out> <genome.fa>                    cell assignments, TSD check, summary
#
# Splitting sites from assign is not bookkeeping: the mini-reference is built from the called sites,
# so the clips can only be assigned once those exist. Discovery searches the genome; assignment
# searches only the called flanks, with no aligner at all - see assign_clips_exact().

source("src/R/integration-sites.R")   # functions only; attaches nothing

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
out   <- args[2]
stopifnot(stage %in% c("extract", "genomic", "genomic-decoy", "sites", "assign"), !is.na(out))
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

} else if (stage == "genomic") {
  # The second, independent read of the same junctions: reads STAR placed on the HOST genome with
  # the vector soft-clipped, which the 3' channel cannot see because it only ever reads one window
  # of one contig. The coordinate is already resolved at MAPQ 255 behind 120-140 bases of host, so
  # no alignment step follows this one.
  #
  # Whole-BAM streaming scan, ~105 minutes a tube, which is why it is its own stage: it is the only
  # expensive part of the pipeline and it must be restartable without redoing the rest. (An earlier
  # comment here said 17 minutes. That was measured with samtools and grep, not with this R path,
  # which is about six times slower; the job was sized from the wrong number and had to go back in
  # the queue.)
  bam_dir <- args[3]
  tubes   <- args[-(1:3)]
  stopifnot(length(tubes) > 0L)

  # One core per tube. Measured on chr19: 27.6M alignments in 3m45s, so a whole BAM is ~104 minutes
  # and eight in series would be ~14 hours - over the walltime, and a timeout inside this stage would
  # leave stages 4 and 5 unrun and nothing on disk. Forked across the tubes it is one tube's wall
  # time, and each tube's result is cached so a restart redoes only what is missing.
  cores <- max(1L, min(length(tubes),
                       as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))), na.rm = TRUE)
  message(">>> scanning ", length(tubes), " BAMs on ", cores, " cores")

  one_tube <- function(tube) {
    cache <- file.path(work, paste0("genomic_", tube, ".csv"))
    if (file.exists(cache)) {
      message("    ", tube, ": cached")
      # Fixed schema, not guessed: the wild-type tubes cache zero rows. See read_genomic_cache().
      return(read_genomic_cache(cache))
    }
    bam <- file.path(bam_dir, paste0("cellranger_", tube), "sample_alignments.bam")
    g <- read_genomic_clips(bam, tube)
    readr::write_csv(g, cache)
    message("    ", tube, ": ", nrow(g), " reads at ", dplyr::n_distinct(g$pos), " positions")
    g
  }
  res <- parallel::mclapply(tubes, one_tube, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("genomic scan failed for: ", paste(tubes[bad], collapse = ", "))
  gen <- dplyr::bind_rows(res)
  readr::write_csv(gen, file.path(work, "genomic_junction_reads.csv"))
  message("=== ", nrow(gen), " genomic junction reads over ", dplyr::n_distinct(gen$pos),
          " distinct positions")
  if (nrow(gen)) {
    top <- gen |>
      dplyr::group_by(chrom, pos, strand) |>
      dplyr::summarise(reads = dplyr::n(), cells = dplyr::n_distinct(paste(tube, barcode)),
                       tubes = dplyr::n_distinct(tube), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(cells))
    print(as.data.frame(utils::head(top, 25)), row.names = FALSE)
  }

} else if (stage == "genomic-decoy") {
  # Error bound for the genome-mapped channel: the same scan with the terminus shuffled to the same
  # base composition. Run as its own stage because it costs what the genomic stage costs, and the
  # pipeline must not wait on a control. An optional 4th argument restricts it to one chromosome,
  # which is how it is run when a full pass is not worth the walltime.
  bam_dir <- args[3]
  region  <- if (length(args) >= 4L && nzchar(args[4])) args[4] else NULL
  tubes   <- readLines(file.path(work, "tubes.txt"))
  cores   <- max(1L, min(length(tubes),
                         as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))), na.rm = TRUE)
  gr <- if (is.null(region)) NULL else {
    fai <- utils::read.delim(paste0(args[5], ".fai"), header = FALSE)
    GenomicRanges::GRanges(region, IRanges::IRanges(1, fai$V2[fai$V1 == region]))
  }
  message(">>> scrambled-terminus decoy over ", length(tubes), " tubes on ", cores, " cores",
          if (is.null(region)) "" else paste0(", restricted to ", region))
  res <- parallel::mclapply(tubes, function(tb) {
    bam <- file.path(bam_dir, paste0("cellranger_", tb), "sample_alignments.bam")
    # Real counts come from the genomic stage's cache where it exists. Rescanning for them doubles
    # the job to recompute a number already on disk.
    cache <- file.path(work, paste0("genomic_", tb, ".csv"))
    tryCatch(genomic_decoy_rate(bam, tb, region = gr,
                                real_cache = if (file.exists(cache)) cache else NULL),
             error = function(e) conditionMessage(e))
  }, mc.cores = cores, mc.preschedule = FALSE)
  ok <- vapply(res, is.data.frame, logical(1))
  if (any(!ok)) message("    FAILED for ", sum(!ok), " tube(s): ",
                        paste(tubes[!ok], collapse = ", "))
  d <- dplyr::bind_rows(res[ok])
  # Region-tagged, so a one-chromosome probe cannot overwrite a whole-genome measurement with a
  # smaller one. Two runs of this stage at different scopes wrote the same filename, which is the
  # silent-wrong-number failure this analysis keeps meeting.
  readr::write_csv(d, file.path(out, if (is.null(region)) "genomic_decoy_rate.csv"
                                     else sprintf("genomic_decoy_rate_%s.csv", region)))
  if (nrow(d)) {
    print(as.data.frame(d[, c("tube", "real_reads", "decoy_reads", "decoy_pos")]), row.names = FALSE)
    message("=== real ", sum(d$real_reads), " reads vs decoy ", sum(d$decoy_reads),
            if (sum(d$decoy_reads) == 0L) "  (as required)" else "  <-- RAISE min_genomic")
  }

} else if (stage == "sites") {
  gtf       <- args[3]
  genome_fa <- args[4]
  tubes     <- readLines(file.path(work, "tubes.txt"))
  stopifnot(file.exists(gtf), file.exists(genome_fa))

  discovery <- read_bwa_hits(file.path(work, "discovery.sam"), channel = "three_prime")
  message(">>> ", nrow(discovery), " discovery clips placed uniquely (MAPQ >= ",
          IS_DEFAULTS$min_mapq, ")")

  # The genome-mapped channel joins DISCOVERY, not just assignment. These reads carry 120-140 bases
  # of host at MAPQ 255, which is far better placement than the 25-50 bp clips above, so requiring
  # them to coincide with an already-called site would add cells but no sites. Clustered together
  # with the 3' clips, they can call a site on their own.
  gen_csv <- file.path(work, "genomic_junction_reads.csv")
  if (file.exists(gen_csv)) {
    gen <- read_genomic_cache(gen_csv)
    if (nrow(gen)) {
      gen <- dplyr::mutate(gen, mapq = as.integer(mapq), seq = NA_character_)
      discovery <- dplyr::bind_rows(
        dplyr::mutate(discovery, channel = "three_prime"),
        dplyr::mutate(dplyr::select(gen, dplyr::any_of(names(discovery))), channel = "genomic"))
      message(">>> plus ", sum(discovery$channel == "genomic"), " genomic junction reads")
    }
  }
  if (!"channel" %in% names(discovery)) discovery$channel <- "three_prime"

  # metadata/, not output/: output/ is gitignored, and a label registry that does not survive a
  # clean checkout cannot pin anything. See label_sites().
  sites      <- label_sites(dplyr::mutate(call_sites(discovery), tier = "discovery"),
                            registry = "metadata/integration_site_registry.csv")
  site_table <- summarise_sites(sites, gtf, tubes)

  # Which channel called each site on its own. Carried in the table because "15 sites" is not one
  # kind of claim: two of them clear the floor only on the two channels' reads summed.
  prov <- site_provenance(discovery, site_table)
  site_table <- dplyr::left_join(site_table, prov, by = "site_id")
  message(">>> called by: ",
          paste(sprintf("%s = %d", names(table(prov$called_by)), table(prov$called_by)),
                collapse = ", "))
  if (any(prov$pooled_only)) {
    message("    POOLED ONLY (weakest evidence, exclude from repertoire counts): ",
            paste(prov$site_id[prov$pooled_only], collapse = ", "))
  }

  # TSD BEFORE the flanks are cut, so it can vouch for a site sitting at the cell floor. The flank
  # gate in write_site_flanks() needs this column; without it the gate is skipped entirely.
  five_sam <- file.path(work, "five_prime.sam")
  if (file.exists(five_sam)) {
    five_pre <- read_bwa_hits(five_sam, channel = "five_prime")
    tsd_pre  <- tsd_check(site_table, five_pre, genome_fa = genome_fa)
    pass     <- unique(tsd_pre$site_id[tsd_pre$tsd %in% c("PASS", "near") & !is.na(tsd_pre$site_id)])
    site_table$tsd_pass <- site_table$site_id %in% pass
    message(">>> 5' channel confirms: ",
            if (length(pass)) paste(sort(pass), collapse = ", ") else "no site")
  }
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

  # No aligner: exact positional matching against the called flanks. See assign_clips_exact() for
  # why bwa left this tier, and for the -T 12 floor it was silently enforcing.
  mini_all <- assign_clips_exact(file.path(work, "assign.fa"), file.path(work, "site_flanks.fa"))

  # The decoy measurement comes first, because the semi_join below drops the decoy flanks. A decoy
  # has no junction, so a non-zero count here is not a rate to report but a bug to fix: the
  # positional anchoring is not doing what it claims and nothing downstream can be trusted.
  eps <- decoy_error_rate(mini_all, n_real = nrow(site_table))
  readr::write_csv(eps, file.path(out, "decoy_error_rate.csv"))
  n_decoy_hits <- sum(startsWith(mini_all$site_id, IS_DEFAULTS$decoy_prefix))
  message("=== decoy flanks hit at the junction base: ", n_decoy_hits,
          if (n_decoy_hits == 0L) "  (as required)" else "  <-- ANCHORING IS BROKEN")
  if (nrow(eps)) print(as.data.frame(eps), row.names = FALSE)

  mini <- dplyr::semi_join(mini_all, dplyr::select(site_table, site_id), by = "site_id") |>
    dplyr::mutate(tier = "assign", channel = "three_prime")

  # The genome-mapped channel's CELLS, folded in here. These reads already called sites back in
  # stage 4; without this they then vanished, and the site table reported clone sizes built from the
  # 3' channel alone - half the evidence at the larger sites. genomic_site_hits() puts each read back
  # on the site it was clustered into, so nothing is assigned twice and nothing is invented.
  #
  # Kept as a separate object from mini: junction_coverage() needs the flank coordinates that only
  # the mini-reference hits carry, so it still gets mini, while the cell table and the site counts
  # get both channels.
  evidence <- mini
  gen_csv  <- file.path(work, "genomic_junction_reads.csv")
  if (file.exists(gen_csv)) {
    gen   <- read_genomic_cache(gen_csv)
    gmini <- genomic_site_hits(gen, site_table)
    message("=== genomic channel folded in: ", nrow(gmini), " reads, ",
            dplyr::n_distinct(paste(gmini$tube, gmini$barcode)), " cells, ",
            dplyr::n_distinct(gmini$site_id), " sites  (",
            nrow(gen) - nrow(gmini), " reads matched no called site)")
    evidence <- dplyr::bind_rows(mini, gmini)
  }

  cell_table <- summarise_cells(evidence)
  readr::write_csv(cell_table, file.path(out, "cell_site_assignments.csv"))

  # The counts in the site table must come from the same place the cell table does.
  site_table <- update_site_counts(site_table, evidence, readLines(file.path(work, "tubes.txt")))
  readr::write_csv(site_table, file.path(out, "integration_sites.csv"))

  # What the two channels say about each other. Written as a table rather than computed in the
  # report, so the estimator and its caveats live next to the code that produced them.
  cr <- channel_capture_recapture(cell_table)
  readr::write_csv(cr, file.path(out, "channel_capture_recapture.csv"))
  if (nrow(cr)) {
    message("=== per-site capture-recapture between the channels")
    print(as.data.frame(cr), row.names = FALSE, digits = 3)
  }

  five      <- read_bwa_hits(file.path(work, "five_prime.sam"), channel = "five_prime")
  genome_fa <- if (length(args) >= 3L) args[3] else NULL

  tsd <- tsd_check(site_table, five, genome_fa = genome_fa)
  readr::write_csv(tsd, file.path(out, "tsd_check.csv"))

  # Orphan 5' junctions, tested at the one position their partner must occupy. Reports the
  # denominator alongside the hits, because a single Delta = 4 match with no idea how many positions
  # were tried is not evidence of anything.
  discovery <- read_bwa_hits(file.path(work, "discovery.sam"), channel = "three_prime")
  resc <- rescue_orphan_junctions(five, discovery, site_table)
  readr::write_csv(resc, file.path(out, "orphan_rescue.csv"))
  if (nrow(resc)) {
    message("=== orphan 5' junctions tested: ", attr(resc, "candidates"),
            "   rescued at delta = 4: ", attr(resc, "rescued"))
    print(as.data.frame(resc), row.names = FALSE)
  }

  # Per-base pile-up either side of each junction, for the report's browser-style figure.
  readr::write_csv(
    junction_coverage(mini, file.path(work, "five_prime.sam"), site_table),
    file.path(out, "junction_coverage.csv"))

  readr::write_csv(
    dplyr::tibble(
      stage = c("discovery clips placed uniquely",
                "discovery clips in a site that passed the cell and UMI floor",
                "sites called",
                "clips assigned against the mini-reference",
                "genome-mapped reads folded in",
                "cell-site assignments",
                "distinct cells labelled",
                "cells at more than one site"),
      value = c(counts$n_discovery, counts$n_in_sites, nrow(site_table), nrow(mini),
                nrow(evidence) - nrow(mini),
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
