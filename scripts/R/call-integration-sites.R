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
#   assign  <out> <genome.fa> [genes.gtf.gz]     cell assignments, TSD check, summary, and - given
#                                                the gtf - what each orphan 5' junction is
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

  write_csv_checked(clip_waterfall(classified),       file.path(out, "clip_waterfall.csv"))
  write_csv_checked(anchor_end_histogram(classified), file.path(out, "anchor_ends.csv"))
  writeLines(tubes, file.path(work, "tubes.txt"))

  # discovery.fa goes to the genome and defines the sites; assign.fa is every placeable clip,
  # discovery ones included, and goes to the mini-reference built from those sites.
  n_disc <- write_clip_fasta(classified, "discovery", file.path(work, "discovery.fa"))
  n_asn  <- write_clip_fasta(classified, c("discovery", "assign"), file.path(work, "assign.fa"))
  message("=== discovery clips: ", n_disc, "   assignable clips: ", n_asn)

  # Every junction-anchored clip that has a barcode and a UMI, INCLUDING the ones below the
  # assignment floor. The floor used to be applied here, by omission: a clip classified "too short
  # to place" was counted in the waterfall and then never written anywhere, so nothing downstream
  # could ask about it and the floor could not be varied without re-extracting. 258 reads, 151
  # molecules, 151 cells of which 102 carry no other evidence - small, but it is the pool every
  # proposal to lower the floor is aimed at, and it has to be on disk to be arguable about.
  n_all <- write_clip_fasta(classified, c("discovery", "assign", "clip too short to place"),
                            file.path(work, "all_clips.fa"))
  message("=== all junction clips, no floor: ", n_all)

  # Cell Ranger's called barcodes, captured with the run so every cell count downstream is
  # reproducible from the outputs alone. See read_called_barcodes() for why this is needed.
  called <- read_called_barcodes(bam_dir, tubes)
  write_csv_checked(called, file.path(work, "called_barcodes.csv.gz"))
  write_csv_checked(count_barcodes(bam_dir, tubes), file.path(work, "barcode_counts.csv"))

  # Vector transcript per tube, independent of junctions. The discriminator for what the completeness
  # denominator should be: a sorted mCherry+ cell with no vector transcript is outside this assay's
  # reach whatever the reason. ~75 s for eight tubes.
  write_csv_checked(vector_expression(bam_dir, tubes), file.path(out, "vector_expression.csv"))
  # Per cell as well as per tube: the per-tube share says whether the sort gate is pure, the per-cell
  # table is what the per-clone expression comparisons and any "transduced but unplaced" question
  # need. Only cells with at least one vector UMI, which is about a sixth of them.
  write_csv_checked(vector_expression_cells(bam_dir, tubes),
                    file.path(out, "vector_expression_cells.csv.gz"))
  message("=== called barcodes: ", nrow(called), " over ", dplyr::n_distinct(called$tube), " tubes")

  # Second channel: the provirus's 5' junction, read off the same BAMs.
  five <- dplyr::bind_rows(lapply(tubes, function(tube) {
    bam <- file.path(bam_dir, paste0("cellranger_", tube), "sample_alignments.bam")
    classify_five_prime_clips(read_five_prime_clips(bam, tube))
  }))
  write_csv_checked(clip_waterfall(five), file.path(out, "clip_waterfall_5p.csv"))
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
    write_csv_checked(g, cache)
    message("    ", tube, ": ", nrow(g), " reads at ", dplyr::n_distinct(g$pos), " positions")
    g
  }
  res <- parallel::mclapply(tubes, one_tube, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("genomic scan failed for: ", paste(tubes[bad], collapse = ", "))
  gen <- dplyr::bind_rows(res)
  write_csv_checked(gen, file.path(work, "genomic_junction_reads.csv"))
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
  write_csv_checked(d, file.path(out, if (is.null(region)) "genomic_decoy_rate.csv"
                                     else sprintf("genomic_decoy_rate_%s.csv", region)))
  if (nrow(d)) {
    print(as.data.frame(d[, c("tube", "real_reads", "decoy_reads", "decoy_pos")]), row.names = FALSE)
    # Read counts, not site calls. A handful of scattered decoy reads is the chance-match rate and
    # is expected; what would matter is a decoy POSITION reaching the cell and UMI floor, because
    # that is what becomes a site. Judge on reads per position, not on the total being non-zero.
    worst <- if (nrow(d)) max(d$decoy_reads / pmax(d$decoy_pos, 1L)) else 0
    message("=== real ", sum(d$real_reads), " reads at ", sum(d$real_pos), " positions vs decoy ",
            sum(d$decoy_reads), " at ", sum(d$decoy_pos), " positions",
            sprintf("  (%.1f%% of reads; worst decoy position %.1f reads)",
                    100 * sum(d$decoy_reads) / max(1L, sum(d$real_reads)), worst))
    if (worst >= IS_DEFAULTS$min_cells) {
      message("    <-- a decoy position is dense enough to be callable; check min_genomic")
    }
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

  # What the floor admits and what it turns away, written BEFORE the floor is applied - otherwise the
  # outputs cannot be used to question it, since call_sites() returns only what passed.
  write_csv_checked(cluster_support(discovery), file.path(out, "cluster_support.csv"))
  write_csv_checked(site_floor_curve(discovery), file.path(out, "site_floor_curve.csv"))

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
  write_csv_checked(site_table, file.path(out, "integration_sites.csv"))
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
  write_csv_checked(eps, file.path(out, "decoy_error_rate.csv"))
  # Decoy hits are NOT required to be zero, and the old assertion that they were only looked right
  # because n_decoy was 100. Chance matching scales with the number of decoy flanks: at 1000 flanks a
  # short clip matches one by chance often, and the expected count runs to tens. The test is observed
  # against expected-under-chance per clip length, which decoy_error_rate() and clip_floor_curve()
  # both report. What would mean broken anchoring is hits at long clip lengths, where chance is
  # negligible - a 20-base clip has a 4^-20 chance per flank and should never hit one.
  n_decoy_hits <- sum(startsWith(mini_all$site_id, IS_DEFAULTS$decoy_prefix))
  long_decoy <- sum(startsWith(mini_all$site_id, IS_DEFAULTS$decoy_prefix) &
                      mini_all$clip_len >= 15L)
  message("=== decoy flanks hit at the junction base: ", n_decoy_hits, " over ",
          IS_DEFAULTS$n_decoy, " decoys; ", long_decoy, " of them on clips >= ",
          IS_DEFAULTS$assign_floor$exact_below, " bp",
          if (long_decoy == 0L) "  (as required)" else "  <-- ANCHORING IS BROKEN")
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
  gmini    <- NULL
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

  # Restrict CELL evidence to Cell Ranger's called cells. The BAM carries corrected CB tags for
  # barcodes that were never called, and they were arriving in the cell table as cells - the ambient
  # droplets, 2.8% of it. The raw clip and genomic tables keep them, because the ambient measurement
  # is computed from exactly those rows.
  called_f <- file.path(work, "called_barcodes.csv.gz")
  if (file.exists(called_f)) {
    called <- readr::read_csv(called_f, show_col_types = FALSE, progress = FALSE)
    before <- dplyr::n_distinct(paste(evidence$tube, evidence$barcode))
    evidence <- dplyr::semi_join(evidence, called, by = c("tube", "barcode"))
    after <- dplyr::n_distinct(paste(evidence$tube, evidence$barcode))
    message("=== called-cell filter: ", before - after, " non-cell barcodes removed (",
            after, " cells remain)")
  } else {
    warning("no called_barcodes.csv.gz; cell counts include non-cell barcodes", call. = FALSE)
  }

  cell_table <- summarise_cells(evidence)
  write_csv_checked(cell_table, file.path(out, "cell_site_assignments.csv"))

  # ---- the clip table: one row per molecule, no floor, every site's mismatch count ----------
  # Built from the same matcher assign_clips_exact() filters, so the two cannot disagree. The gate
  # below is the whole point of writing it this way.
  # Needed by unassigned_molecules() to separate "placed uniquely at a sub-threshold locus" from
  # "could not be placed at all", which is the difference between a missed site and a blind spot.
  discovery_for_unassigned <- read_bwa_hits(file.path(work, "discovery.sam"),
                                            channel = "three_prime")

  all_fa <- file.path(work, "all_clips.fa")
  if (file.exists(all_fa) && file.size(all_fa) > 0L) {
    raw <- clip_site_matches(all_fa, file.path(work, "site_flanks.fa"), min_len = 1L)
    clip_tbl <- molecule_clip_table(raw, gen_site_hits = gmini, sites = site_table)
    # Gzipped for the same reason as the matches below: on a fileset with no headroom, a 250 KB
    # write truncates where a 40 KB one lands. readr reads .csv.gz transparently; zcat it to look.
    write_csv_checked(clip_tbl, file.path(out, "clip_table.csv.gz"))
    # Gzipped: it is the largest thing this stage writes (~60k rows), it is a debugging artefact
    # rather than something the report reads, and the fileset this lands on has no room to spare.
    # readr reads .csv.gz transparently, so nothing downstream needs to know.
    write_csv_checked(dplyr::filter(raw, nm <= budget),
                      file.path(out, "clip_matches.csv.gz"))
    message("=== clip table: ", nrow(clip_tbl), " molecules (",
            sum(clip_tbl$channels == "three_prime"), " 3' only, ",
            sum(clip_tbl$channels == "genomic"), " genomic only, ",
            sum(clip_tbl$channels == "both"), " both), ",
            sum(clip_tbl$clip_len < IS_DEFAULTS$assign_floor$min_len, na.rm = TRUE),
            " below the floor")

    # THE EQUALITY GATE. Filtered at the current floor and tolerance, the clip table must reproduce
    # the pipeline's own output exactly. If it does not, the table is wrong and nothing built on it
    # is interpretable - so this stops the run rather than warning.
    gate <- clip_table_gate(raw, site_table, cell_table, gen_site_hits = gmini,
                            called = if (exists("called")) called else NULL)
    print(as.data.frame(gate), row.names = FALSE)
    if (!all(gate$agrees)) {
      stop("clip table does not reproduce the pipeline output; see the gate table above")
    }
    message("=== clip table gate: PASS")

    # Per-cell detection and detection-corrected clone size. Uses the clip table's molecule counts,
    # which is the only place they are both de-duplicated across channels and free of the floor.
    det <- detection_model(cell_table, clip_tbl)
    write_csv_checked(det, file.path(out, "detection_model.csv"))

    # Every number the report quotes about floors, ambient and completeness is computed here, not
    # typed into the prose. A figure in a rendered report that cannot be regenerated from the code
    # and the data is a figure nobody can check.
    write_csv_checked(clip_floor_curve(clip_tbl, sites = site_table),
                      file.path(out, "clip_floor_curve.csv"))
    write_csv_checked(unassigned_molecules(clip_tbl, discovery_for_unassigned, cell_table),
                      file.path(out, "unassigned_molecules.csv"))
    counts_f <- file.path(work, "barcode_counts.csv")
    if (file.exists(called_f) && file.exists(counts_f)) {
      write_csv_checked(
        ambient_rate(clip_tbl, called, readr::read_csv(counts_f, show_col_types = FALSE)),
        file.path(out, "ambient_rate.csv"))
    }
    mch <- if (file.exists(called_f)) sum(called$tube %in% c("Rx_4", "Rx_6")) else NA_integer_
    vec_f <- file.path(out, "vector_expression.csv")
    vpos <- if (file.exists(vec_f)) {
      v <- readr::read_csv(vec_f, show_col_types = FALSE)
      sum(v$vector_positive[v$tube %in% c("Rx_4", "Rx_6")], na.rm = TRUE)
    } else NA_integer_
    write_csv_checked(repertoire_completeness(det, mch, vpos),
                      file.path(out, "repertoire_completeness.csv"))
    if (nrow(det)) {
      message("=== detection model (NA where the truncated mean is too close to 1 to fit)")
      print(as.data.frame(det |> dplyr::transmute(
        site_id, cells, mean_molecules = round(mean_molecules, 2), lambda = round(lambda, 3),
        detection = round(detection, 3), clone_size, clone_lo, clone_hi,
        obs_1, exp_1 = round(exp_1), obs_2, exp_2 = round(exp_2),
        var_mean = round(var_mean_ratio, 2))), row.names = FALSE)
    }
  }

  # The counts in the site table must come from the same place the cell table does.
  site_table <- update_site_counts(site_table, evidence, readLines(file.path(work, "tubes.txt")))
  write_csv_checked(site_table, file.path(out, "integration_sites.csv"))

  # What the two channels say about each other. Written as a table rather than computed in the
  # report, so the estimator and its caveats live next to the code that produced them.
  cr <- channel_capture_recapture(cell_table)
  write_csv_checked(cr, file.path(out, "channel_capture_recapture.csv"))
  if (nrow(cr)) {
    message("=== per-site capture-recapture between the channels")
    print(as.data.frame(cr), row.names = FALSE, digits = 3)
  }

  five      <- read_bwa_hits(file.path(work, "five_prime.sam"), channel = "five_prime")
  genome_fa <- if (length(args) >= 3L) args[3] else NULL
  # Optional 4th argument: the gtf, for the orphan-junction verdicts below. Optional rather than
  # required so an existing invocation of this stage keeps working - it only costs those verdicts.
  gtf_path  <- if (length(args) >= 4L) args[4] else NA_character_

  tsd <- tsd_check(site_table, five, genome_fa = genome_fa)
  write_csv_checked(tsd, file.path(out, "tsd_check.csv"))

  # One row per site with every flag that bears on trusting it, so downstream work can filter on
  # `defensible` instead of re-deriving a judgement this analysis has now made three times. Placed
  # here because it needs both tsd and the detection model, and an earlier placement referenced tsd
  # before it existed - caught by checking the statement order rather than by running it.
  write_csv_checked(
    site_evidence(site_table, tsd = tsd,
                  detection = if (exists("det")) det else NULL,
                  flank_fa = file.path(work, "site_flanks.fa")),
    file.path(out, "site_evidence.csv"))

  # Orphan 5' junctions, tested at the one position their partner must occupy. Reports the
  # denominator alongside the hits, because a single Delta = 4 match with no idea how many positions
  # were tried is not evidence of anything.
  discovery <- read_bwa_hits(file.path(work, "discovery.sam"), channel = "three_prime")
  resc <- rescue_orphan_junctions(five, discovery, site_table)
  write_csv_checked(resc, file.path(out, "orphan_rescue.csv"))
  if (nrow(resc)) {
    message("=== orphan 5' junctions tested: ", attr(resc, "candidates"),
            "   rescued at delta = 4: ", attr(resc, "rescued"))
    print(as.data.frame(resc), row.names = FALSE)
  }

  # What each orphan actually IS. "No 3' partner" is the absence of a finding rather than one, and
  # these turned out to be answerable from evidence already on disk: one is a spliced host-vector
  # fusion from a provirus already called, one reaches the cell and UMI floor once the two channels
  # are pooled, and the rest are too thin to separate. The gtf and the hashtag calls are optional, so
  # the stage still runs where they are absent - with the clonality column left NA and said so.
  if (!is.na(gtf_path) && file.exists(gtf_path) && nrow(resc)) {
    ht_path <- file.path(work, "hashtag_calls.csv.gz")
    cb_path <- file.path(work, "called_barcodes.csv.gz")
    ht <- if (file.exists(ht_path)) {
      dplyr::mutate(readr::read_csv(ht_path, show_col_types = FALSE),
                    barcode = sub("^.*_([ACGT]+-1)$", "\\1", cell))
    } else NULL
    cb <- if (file.exists(cb_path)) readr::read_csv(cb_path, show_col_types = FALSE) else NULL
    if (is.null(ht)) {
      message("!!! no hashtag calls: orphan clonality left untested ",
              "(run scripts/R/export-hashtag-calls.R)")
    }
    oev <- orphan_junction_evidence(five, discovery, site_table, gtf_path, genome_fa,
                                    hashtag_calls = ht, called_barcodes = cb)
    write_csv_checked(oev, file.path(out, "orphan_evidence.csv"))
    message("=== orphan verdicts")
    print(as.data.frame(oev[, c("chrom", "pos_5p", "pooled_cells", "pooled_umis", "meets_floor",
                                "exon_end_dist", "donor", "clip_beyond", "vector_acceptor",
                                "replicate", "replicate_p", "verdict")]), row.names = FALSE)
  }

  # Per-base pile-up either side of each junction, for the report's browser-style figure.
  write_csv_checked(
    junction_coverage(mini, file.path(work, "five_prime.sam"), site_table),
    file.path(out, "junction_coverage.csv"))

  write_csv_checked(
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
