# Provenance of a report: which code, environment, input data and pipeline steps produced it.
# Run bookkeeping, not analysis. Used by notebooks/_provenance.qmd (reports) and _targets.R.
#
#   prov <- collect_provenance()     1. collect: every fact in one list, nothing printed
#   cat(render_summary(prov))        2. display: markdown for the provenance checklist and the callouts
#
# Everything degrades gracefully: without git, targets or renv the facts become NA, never an error.

# ---- helpers ------------------------------------------------------------------------------------

# Run git; return all output lines (character(0) on failure)
git_lines <- function(...) {
  tryCatch(suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
           error = function(e) character())
}
git <- function(...) {                      # first line only, or NA
  out <- git_lines(...)
  if (length(out) == 0) NA_character_ else out[1]
}

md5_of_lines <- function(lines) {           # md5 of a character vector, via a temporary file
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(lines, f)
  unname(tools::md5sum(f))
}

human_size <- function(bytes) {
  vapply(bytes, function(b) if (is.na(b)) "" else
    format(structure(b, class = "object_size"), units = "auto"), character(1))
}

# targets reads _targets.R in a separate R process (it can't run inside knitr while this report
# renders); keep that process's console output (e.g. renv startup notes) out of the report.
quiet_callr <- list(show = FALSE)

# ---- 1. collect ---------------------------------------------------------------------------------

#' The commit the code is at, named after the nearest tag (e.g. v0.2, or v0.2-3-g8f2c1ab three commits later).
git_version <- function() {
  v <- git("describe", "--tags", "--always")
  if (is.na(v)) "not a git repository" else v
}

#' Web link to the current commit on GitHub (NA if there is no GitHub remote). Works once pushed.
git_commit_url <- function() {
  sha <- git("rev-parse", "HEAD")
  remote <- git("remote", "get-url", "origin")
  if (is.na(sha) || is.na(remote) || !grepl("github\\.com", remote)) return(NA_character_)
  repo <- sub("\\.git$", "", sub("^git@github\\.com:", "https://github.com/", remote))
  paste0(repo, "/commit/", sha)
}

#' Files that differ from the commit (modified, or new and not yet added). character(0) when clean.
uncommitted_files <- function(path = ".") {
  sub("^.. ", "", git_lines("status", "--porcelain", "--", path))
}

#' The path, if git tracks it; character(0) otherwise. Used to decide whether git can vouch for a
#' declared input that is neither data nor pipeline output - a shared config or a sourced script.
git_tracked <- function(path) {
  git_lines("ls-files", "--error-unmatch", "--", path)
}

#' One declared input against the manifest the pipeline wrote beside it.
#'
#' Read with every column as character: a commit hash is often a valid number, and the short form of
#' `8e55518` parses as 8 x 10^55518, i.e. Inf. The manifest records the full SHA for that reason.
manifest_check <- function(path, md5, manifest = "manifest.tsv") {
  man <- file.path(dirname(path), manifest)
  if (!file.exists(man)) return("no pipeline manifest")
  m <- tryCatch(read.delim(man, colClasses = "character", na.strings = character()),
                error = function(e) NULL)
  if (is.null(m) || !all(c("file", "md5") %in% names(m))) return("unreadable pipeline manifest")
  row <- m[m$file == basename(path), , drop = FALSE]
  if (nrow(row) != 1L) return("not in pipeline manifest")
  if (identical(row$md5[[1]], md5)) "matches pipeline manifest" else "DIFFERS from pipeline manifest"
}

#' TRUE / FALSE: do the installed R packages match renv.lock? NA if renv isn't active.
renv_synced <- function() {
  if (!requireNamespace("renv", quietly = TRUE) || !nzchar(Sys.getenv("RENV_PROJECT"))) return(NA)
  s <- tryCatch({ invisible(utils::capture.output(x <- renv::status())); x }, error = function(e) NULL)
  if (is.null(s)) NA else isTRUE(s$synchronized)
}

#' Path of the report being rendered, relative to the project root.
report_input <- function() {
  input <- normalizePath(knitr::current_input(dir = TRUE))
  input <- sub("\\.rmarkdown$", ".qmd", input)   # knitr sees Quarto's intermediate copy, not the .qmd
  sub(paste0("^", normalizePath(getwd()), "/"), "", input)
}

#' Name of the pipeline target that renders this report (the tar_quarto() whose command names this file).
#' (targets::tar_name() can't tell: Quarto runs the report's code in a separate R process.)
report_target <- function() {
  m <- pipeline_manifest()
  if (!nrow(m)) return(NA_character_)
  hit <- m$name[grepl(report_input(), m$command, fixed = TRUE)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' The manifest of the pipeline as defined, computed once per session (each call costs a callr
#' subprocess, and three of the functions below want it).
.manifest_cache <- new.env(parent = emptyenv())
pipeline_manifest <- function() {
  if (is.null(.manifest_cache$m)) {
    .manifest_cache$m <- tryCatch(
      targets::tar_manifest(fields = c("name", "command"), callr_arguments = quiet_callr),
      error = function(e) data.frame(name = character(), command = character()))
  }
  .manifest_cache$m
}

#' The files a report declares as inputs, taken from its target's command.
#'
#' tar_meta() is the better source and is used whenever it has the target. But a report target that
#' has never been built by tar_make() has no meta row, and a report reading the output of a pipeline
#' run outside targets - by sbatch, say - is exactly the kind that gets rendered directly, so this is
#' the common case rather than the edge one. Building the target just to populate meta is not a
#' substitute: for this report it would rebuild twenty upstream targets, the whole Seurat chain among
#' them, which is a great deal of compute to make a banner green.
#'
#' The command text is then the only record of what the report declared, so the quoted paths are
#' read out of it. Anything that is not a file on disk is dropped, which also drops the target names
#' and arguments that share the quoting.
declared_inputs <- function(target) {
  m <- pipeline_manifest()
  row <- m[m$name == target, , drop = FALSE]
  if (nrow(row) != 1L) return(character())
  quoted <- regmatches(row$command[[1]], gregexpr('"[^"]+"', row$command[[1]]))[[1]]
  paths <- gsub('"', "", quoted)
  unique(paths[file.exists(paths) & !dir.exists(paths)])
}

#' Every target `target` depends on, directly or indirectly (not the target itself).
upstream_targets <- function(target) {
  edges <- targets::tar_network(targets_only = TRUE, callr_arguments = quiet_callr)$edges
  found <- character()
  todo <- target
  while (length(todo)) {
    parents <- setdiff(edges$from[edges$to %in% todo], found)
    found <- c(found, parents)
    todo <- parents
  }
  found
}

#' md5 of a file; reuses the auto data catalogue's checksum when size and modification time match,
#' so large inputs are not re-hashed on every render.
file_md5 <- function(path, auto = "metadata/auto_catalogue_data.tsv") {
  if (startsWith(path, "data/") && file.exists(auto)) {
    a <- read.delim(auto, colClasses = "character", na.strings = character())
    row <- a[a$file == sub("^data/", "", path), , drop = FALSE]
    info <- file.info(path)
    if (nrow(row) == 1 && nzchar(row$md5) &&
        row$modified == format(info$mtime, "%Y-%m-%d %H:%M:%S") && row$size == human_size(info$size)) {
      return(row$md5)
    }
  }
  unname(tools::md5sum(path))
}

#' The input files this report's pipeline read, checked against whatever can vouch for each one.
#' Files in data/ must be in the reviewed catalogue with the same md5; files under output/ against
#' the manifest the pipeline wrote beside them; anything else tracked in git is covered by git.
#'
#' Two things here were wrong and made the check vacuous for a report whose inputs are produced
#' outside targets. Only the FIRST path of each file target was taken, and a tar_quarto() target is
#' one file target whose paths are the rendered output, the sources and every declared input - so all
#' but one input was dropped. And only UPSTREAM targets were consulted, which for such a report is
#' none of them. The report's own target is now included, and its own output and sources are removed
#' instead (they are what the Code and Pipeline checks are for), leaving exactly its inputs.
input_files <- function(meta, catalogue = "metadata/catalogue_data.tsv") {
  files <- meta[meta$format %in% "file", , drop = FALSE]
  # as.character(), because unlist() of an empty list is NULL rather than character(0) and every
  # filter below is a string operation - collect_provenance() calls this with no files at all to
  # build the empty frame, and that call is the one that reached startsWith(NULL, ...).
  paths <- unique(as.character(unlist(files$path, use.names = FALSE)))
  paths <- paths[!is.na(paths)]
  if (length(paths) == 0) {
    return(data.frame(path = character(), kind = character(), size = character(), md5 = character(),
                      check = character()))
  }
  # Inputs only: not the report's own rendered output, not its sources, and not the code and config
  # it declares alongside them (git-tracked, outside data/, metadata/ and output/) - a sourced .R
  # file and renv.lock are what the Code and R-packages checks are for, and listing them here as
  # "input data" is noise that makes the rows that matter harder to see.
  paths <- paths[!startsWith(paths, "output/reports/") & !grepl("\\.qmd$|\\.rmarkdown$", paths)]
  is_data <- startsWith(paths, "data/") | startsWith(paths, "metadata/") |
             startsWith(paths, "output/")
  paths <- paths[is_data | !vapply(paths, function(p) length(git_tracked(p)) > 0, logical(1))]
  if (length(paths) == 0) {
    return(data.frame(path = character(), kind = character(), size = character(), md5 = character(),
                      check = character()))
  }
  cat_tab <- if (file.exists(catalogue)) {
    read.delim(catalogue, colClasses = "character", na.strings = character())
  } else data.frame(file = character(), kind = character(), md5 = character())
  rows <- lapply(paths, function(p) {
    # A directory target - Cell Ranger's per-sample matrix folder is one - has no md5 of its own, and
    # md5sum() returning NA for it was once reported as "missing": the one wrong thing a provenance
    # check can say about a file sitting right there, because it reads as data loss. That became
    # "directory (not hashed)", which was honest and still an admission - eight of thirty-five
    # declared inputs unverified. It is now DIGESTED, over its contents; see dir_digest(). The cost
    # was measured rather than assumed: the eight matrix directories are 1.7 GB over 24 files and
    # take about five seconds. The BAMs beside them are 426 GB, which is what dir_digest_max_bytes()
    # refuses.
    if (dir.exists(p)) {
      # A DIRECTORY IS HASHED BY ITS CONTENTS, not skipped. These entries used to read
      # "directory (not hashed)" - amber, honest, and the one place where this report's own standard
      # (every number reproducible from data and code that is available) was asserted rather than
      # demonstrated. Eight of the thirty-five declared inputs were Cell Ranger matrix directories,
      # so a ninth of the provenance table was an admission.
      #
      # md5sum of a directory returns NA, which is why it was skipped. The fix is a digest OVER the
      # directory: every file's own md5, paired with its relative path, sorted, and hashed. Sorting
      # is what makes it reproducible - readdir order is not stable across filesystems - and
      # including the paths is what makes it sensitive to a file being renamed or moved rather than
      # only to its bytes changing.
      #
      # The cost is bounded and was measured before this was written: the eight matrix directories
      # are 1.7 GB over 24 files, about a minute. It is NOT a licence to hash anything: the BAMs
      # under data/ are 426 GB, and `max_bytes` refuses rather than quietly spending an hour. A
      # directory too large to digest still reports its file count, as before.
      fs <- list.files(p, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
      fs <- fs[!dir.exists(fs)]
      n <- length(fs)
      size_txt <- sprintf("%d file%s", n, if (n == 1L) "" else "s")
      total <- sum(file.info(fs)$size, na.rm = TRUE)
      if (!n || total > dir_digest_max_bytes()) {
        return(data.frame(path = p, kind = "directory",
                          size = size_txt, md5 = NA_character_,
                          check = if (!n) "directory (empty)" else
                            sprintf("directory (%.0f GB, over the digest limit)", total / 1e9)))
      }
      md5 <- dir_digest(p, fs)
      # The catalogue holds FILES, so a directory never has a row of its own. Verify it through its
      # contents instead: every file inside, looked up in the catalogue by its path under data/. The
      # directory passes only if every catalogued file matches, which is a stronger statement than a
      # single digest comparison would be - it names which file differs.
      rel <- sub("^data/", "", fs)
      rows <- cat_tab[cat_tab$file %in% rel, , drop = FALSE]
      with_md5 <- rows[nzchar(rows$md5), , drop = FALSE]
      check <- if (!startsWith(p, "data/")) "digested (outside data/)"
        else if (nrow(rows) == 0) "digested, not catalogued"
        else if (nrow(with_md5) == 0) "digested; catalogue has no md5 yet"
        else {
          now <- vapply(file.path("data", with_md5$file), function(f) unname(file_md5(f)),
                        character(1), USE.NAMES = FALSE)
          bad <- with_md5$file[!is.na(now) & now != with_md5$md5]
          if (length(bad)) paste0("DIFFERS from catalogue: ", paste(basename(bad), collapse = ", "))
          else if (nrow(with_md5) < length(fs))
            sprintf("matches catalogue (%d of %d files)", nrow(with_md5), length(fs))
          else "matches catalogue"
        }
      return(data.frame(path = p, kind = "directory",
                        size = sprintf("%s, %.0f MB", size_txt, total / 1e6),
                        md5 = md5, check = check))
    }
    md5 <- if (file.exists(p)) file_md5(p) else NA_character_
    if (startsWith(p, "data/")) {
      entry <- cat_tab[cat_tab$file == sub("^data/", "", p), , drop = FALSE]
      kind <- if (nrow(entry) == 1 && nzchar(entry$kind)) entry$kind else ""
      check <- if (is.na(md5)) "missing"
        else if (nrow(entry) == 0) "not catalogued"
        else if (!nzchar(entry$md5)) "catalogued (no md5)"
        else if (entry$md5 == md5) "matches catalogue"
        else "DIFFERS from catalogue"
    } else if (file.exists(file.path(dirname(p), "manifest.tsv"))) {
      # Produced by a pipeline outside targets, so neither the catalogue nor git can vouch for it.
      # Its own run can: the stage that wrote it recorded size and md5 in manifest.tsv beside it, in
      # the same call that line-counted the write. A file whose md5 is not the recorded one is not
      # the pipeline's output, which is the failure that otherwise renders as a normal report.
      #
      # Keyed on the manifest being there rather than on the path starting with output/, so that a
      # pipeline writing anywhere else is checked the same way. Testing this with a file in a
      # scratch directory is what showed the path prefix silently skipping the check.
      kind  <- "pipeline output"
      check <- if (is.na(md5)) "missing" else manifest_check(p, md5)
    } else if (startsWith(p, "metadata/") || length(git_tracked(p))) {
      kind <- if (startsWith(p, "metadata/")) "metadata" else "code"
      check <- if (is.na(md5)) "missing" else if (length(uncommitted_files(p))) "uncommitted" else "in git"
    } else {
      kind <- ""
      check <- if (is.na(md5)) "missing" else "outside data/, metadata/ and output/"
    }
    data.frame(path = p, kind = kind, size = human_size(file.info(p)$size), md5 = md5, check = check)
  })
  do.call(rbind, rows)
}

#' Everything the report's provenance shows, in one list. Prints nothing.
collect_provenance <- function(catalogue = "metadata/catalogue_data.tsv") {
  p <- list(
    version          = git_version(),
    commit_url       = git_commit_url(),
    uncommitted      = uncommitted_files(),
    rendered         = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    host             = Sys.info()[["nodename"]],
    renv_synced      = renv_synced(),
    built_by_targets = isTRUE(targets::tar_active()),
    report           = report_target(),
    catalogue        = catalogue,
    upstream         = character(),
    fingerprint      = NA_character_,
    inputs           = input_files(data.frame(format = character(), path = I(list())))
  )
  if (!is.na(p$report)) {
    p$upstream <- upstream_targets(p$report)
    meta <- targets::tar_meta(names = tidyselect::any_of(p$upstream),
                              fields = tidyselect::any_of(c("data", "format", "path")))
    # One hash of every upstream result: identical fingerprints = identical upstream state, anywhere
    p$fingerprint <- substr(md5_of_lines(sort(paste(meta$name, meta$data))), 1, 12)
    # The report's own target as well, because a report whose inputs are produced outside targets
    # declares them on its own tar_quarto() target and has no upstream file target at all. When that
    # target has never been built there is no meta row either, and the declaration is read from the
    # command instead - see declared_inputs().
    own <- targets::tar_meta(names = tidyselect::any_of(p$report),
                             fields = tidyselect::any_of(c("data", "format", "path")))
    if (!nrow(own)) {
      own <- data.frame(format = "file", path = I(list(declared_inputs(p$report))))
    }
    p$inputs <- input_files(rbind(meta[, names(own)], own), catalogue)
  }

  # One row per check: "ok", "warning" (check before trusting) or "error" (can't be trusted as rendered)
  # A digested directory is verified, or explicitly without a baseline - not an error. Without
  # these the banner reads every Cell Ranger matrix directory as a failure, which is the
  # opposite of what hashing them was for. The startsWith() covers the partial-match string,
  # which sprintf() builds with a file count in it.
  inputs_ok <- startsWith(p$inputs$check, "matches catalogue") |
    p$inputs$check %in% c("in git", "matches pipeline manifest",
                          "digested, not catalogued", "digested; catalogue has no md5 yet",
                          "digested (outside data/)")
  p$checks <- data.frame(
    check = c("Code", "R&nbsp;packages", "Pipeline", "Input&nbsp;data"),
    level = c(
      if (length(p$uncommitted)) "warning" else "ok",
      if (isTRUE(p$renv_synced)) "ok" else "warning",
      if (p$built_by_targets && !is.na(p$report)) "ok" else "warning",
      # all() of an empty vector is TRUE, so a report with no declared inputs used to come out
      # GREEN here - the one state that should never read as a pass, because nothing was checked.
      # A file that is not the one the pipeline wrote is an error, not a warning: every number in the
      # report downstream of it is then computed from something unaccounted for.
      if (any(p$inputs$check %in% c("DIFFERS from catalogue", "DIFFERS from pipeline manifest",
                                    "missing"))) "error"
        else if (nrow(p$inputs) == 0L || !all(inputs_ok)) "warning" else "ok"
    )
  )
  p$state <- if (any(p$checks$level == "error")) "error" else if (any(p$checks$level == "warning")) "warning" else "ok"
  p
}

# ---- 2. display ---------------------------------------------------------------------------------

commit_link <- function(p) {
  if (is.na(p$commit_url)) paste0("`", p$version, "`") else paste0("[`", p$version, "`](", p$commit_url, ")")
}

catalogue_link <- function(p) {
  if (is.na(p$commit_url)) return(paste0("`", p$catalogue, "`"))
  paste0("[`", p$catalogue, "`](", sub("/commit/", "/blob/", p$commit_url), "/", p$catalogue, ")")
}

#' The always-visible provenance box: one row per check with a ✔ / ✘, coloured by the overall state.
render_summary <- function(p) {
  mark <- c(ok      = '<span style="color:#198754">&#10003;</span>',   # green tick (U+2713: plain text, so it takes colour; U+2714 renders as an emoji)
            warning = '<span style="color:#fd7e14">&#10008;</span>',   # amber cross: check before trusting
            error   = '<span style="color:#dc3545">&#10008;</span>')   # red cross: can't be trusted as rendered
  code <- if (length(p$uncommitted)) {
    sprintf("%d uncommitted file(s) since %s: %s", length(p$uncommitted), commit_link(p),
            paste0("`", p$uncommitted, "`", collapse = ", "))
  } else paste("committed:", commit_link(p))
  # The analysis code is a PACKAGE now, not a file in this repository, so the commit above says
  # nothing about which version of it produced the numbers. targets cannot watch a package either,
  # which means a package change does not re-render the report. Recording the version is the least
  # a reader needs in order to know what they are looking at.
  pkg <- tryCatch(as.character(utils::packageVersion("scarab")), error = function(e) NA_character_)
  if (!is.na(pkg)) code <- paste0(code, " &middot; `scarab` ", pkg)
  renv <- if (isTRUE(p$renv_synced)) "match `renv.lock`"
          else if (isFALSE(p$renv_synced)) "differ from `renv.lock`" else "unknown (renv not active)"
  pipeline <- if (is.na(p$report)) "not a target in `_targets.R`: no pipeline information"
              else if (!p$built_by_targets) "rendered outside the pipeline: upstream results may be out of date"
              else paste0("built by targets &middot; upstream fingerprint `", p$fingerprint, "`")
  bad <- p$inputs[!p$inputs$check %in% c("matches catalogue", "in git",
                                         "matches pipeline manifest"), , drop = FALSE]
  # "0 of 0 file(s) verified" is worse than a failure: it reads as a pass when in fact nothing was
  # checked, which is what a report whose inputs are not declared as targets produces. Say so.
  inputs <- if (nrow(p$inputs) == 0L) {
    "**no input files declared** - nothing was verified, which is not the same as a pass"
  } else {
    sprintf("%d of %d file(s) verified (details below)", nrow(p$inputs) - nrow(bad), nrow(p$inputs))
  }
  if (nrow(bad)) inputs <- paste0(inputs, ": ", paste0("`", bad$path, "` ", bad$check, collapse = "; "))

  tab <- data.frame(mark = mark[p$checks$level], check = p$checks$check,
                    detail = c(code, renv, pipeline, inputs))
  callout <- c(ok = "callout-tip", warning = "callout-warning", error = "callout-important")[[p$state]]
  paste0("::: {.", callout, " collapse=\"true\"}\n## Provenance\n\n",   # foldable, starts folded; its colour stays visible
         paste(knitr::kable(tab, col.names = c("", "", ""), escape = FALSE), collapse = "\n"),
         "\n\nRendered ", p$rendered, " on `", p$host, "`. R and package versions: session info at the end.\n:::\n\n")
}

#' Input data callout body: the files this report's pipeline read, checked against the catalogue.
render_inputs <- function(p) {
  head <- sprintf("%d file(s) read by this report's pipeline; data files are described in %s.\n\n",
                  nrow(p$inputs), catalogue_link(p))
  if (nrow(p$inputs) == 0) return(head)
  tab <- transform(p$inputs, path = paste0("`", path, "`"), md5 = paste0("`", substr(md5, 1, 12), "`"))
  paste0(head, paste(knitr::kable(tab, escape = FALSE), collapse = "\n"), "\n\n")
}

#' Static Mermaid diagram of this report and its upstream targets, as a Quarto mermaid block.
pipeline_mermaid <- function(p) {
  lines <- if (is.na(p$report)) {
    targets::tar_mermaid(targets_only = TRUE, outdated = FALSE, legend = FALSE, callr_arguments = quiet_callr)
  } else {
    targets::tar_mermaid(names = tidyselect::any_of(!!p$report), targets_only = TRUE,
                         outdated = FALSE, legend = FALSE, callr_arguments = quiet_callr)
  }
  # drop the wrapping "Graph" box targets draws around everything
  lines <- lines[!grepl("^\\s*(subgraph Graph|end|direction LR|style Graph)\\b", lines)]
  paste(c("```{mermaid}", lines, "```"), collapse = "\n")
}

#' Interactive graph of this report and its upstream targets, sized to its widest layer.
pipeline_graph <- function(p) {
  g <- if (is.na(p$report)) {
    targets::tar_visnetwork(targets_only = TRUE, outdated = FALSE, callr_arguments = quiet_callr)
  } else {
    # !! injects the name now: the selection is evaluated in the separate process
    targets::tar_visnetwork(names = tidyselect::any_of(!!p$report), targets_only = TRUE,
                            outdated = FALSE, callr_arguments = quiet_callr)
  }
  per_layer <- max(table(g$x$nodes$level), 1)
  g$height <- paste0(max(220, 70 * per_layer + 60), "px")
  g
}

#' Folder for the figures and tables a report saves: artifacts/ next to its HTML in output/reports/latest/.
#' notebooks/<report>/<report>.qmd -> output/reports/latest/notebooks/<report>/artifacts/
#' Emptied at the start of each render, so it only ever holds what the latest render produced.
#' Call it once per report (in the setup chunk).
report_dir <- function() {
  dir <- file.path("output", "reports", "latest", dirname(report_input()), "artifacts")
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

#' A reproducible digest over a directory's contents.
#'
#' Every file's md5 paired with its path relative to the directory, sorted, then hashed. Sorting is
#' what makes it stable - readdir order is not guaranteed across filesystems - and including the
#' relative paths is what makes it notice a rename, which a digest of the concatenated bytes would
#' not.
#'
#' Not recursive into symlinked directories beyond what list.files() follows, and not sensitive to
#' mode or mtime: this answers "are these the same contents", not "is this the same filesystem
#' state".
dir_digest <- function(dir, files = NULL) {
  if (is.null(files)) {
    files <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    files <- files[!dir.exists(files)]
  }
  if (!length(files)) return(NA_character_)
  rel <- sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", dir), "/?"), "", files)
  o <- order(rel)
  per <- vapply(files[o], function(f) unname(tools::md5sum(f)), character(1), USE.NAMES = FALSE)
  unname(openssl::md5(paste(rel[o], per, collapse = "\n")))
}

#' How large a directory may be before it is reported rather than digested.
#'
#' A guard, not a preference. The matrix directories this exists for are 1.7 GB in total; the BAMs
#' beside them are 426 GB, and a provenance banner must not quietly spend an hour reading them.
#' Override with the SCARAB_DIR_DIGEST_MAX_GB environment variable.
dir_digest_max_bytes <- function() {
  g <- suppressWarnings(as.numeric(Sys.getenv("SCARAB_DIR_DIGEST_MAX_GB", "8")))
  if (is.na(g) || g <= 0) g <- 8
  g * 1e9
}
