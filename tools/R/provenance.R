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
  m <- tryCatch(targets::tar_manifest(fields = c("name", "command"), callr_arguments = quiet_callr),
                error = function(e) NULL)
  if (is.null(m)) return(NA_character_)
  hit <- m$name[grepl(report_input(), m$command, fixed = TRUE)]
  if (length(hit) == 0) NA_character_ else hit[1]
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

#' The input files this report's pipeline read (its upstream file targets), checked against the catalogue.
#' Files in data/ must be in the reviewed catalogue with the same md5; files in metadata/ are covered by git.
input_files <- function(meta, catalogue = "metadata/catalogue_data.tsv") {
  files <- meta[meta$format %in% "file", , drop = FALSE]
  paths <- vapply(files$path, function(p) p[[1]], character(1))
  if (length(paths) == 0) {
    return(data.frame(path = character(), kind = character(), size = character(), md5 = character(),
                      check = character()))
  }
  cat_tab <- if (file.exists(catalogue)) {
    read.delim(catalogue, colClasses = "character", na.strings = character())
  } else data.frame(file = character(), kind = character(), md5 = character())
  rows <- lapply(paths, function(p) {
    md5 <- if (file.exists(p)) file_md5(p) else NA_character_
    if (startsWith(p, "data/")) {
      entry <- cat_tab[cat_tab$file == sub("^data/", "", p), , drop = FALSE]
      kind <- if (nrow(entry) == 1 && nzchar(entry$kind)) entry$kind else ""
      check <- if (is.na(md5)) "missing"
        else if (nrow(entry) == 0) "not catalogued"
        else if (!nzchar(entry$md5)) "catalogued (no md5)"
        else if (entry$md5 == md5) "matches catalogue"
        else "DIFFERS from catalogue"
    } else if (startsWith(p, "metadata/")) {
      kind <- "metadata"
      check <- if (is.na(md5)) "missing" else if (length(uncommitted_files(p))) "uncommitted" else "in git"
    } else {
      kind <- ""
      check <- if (is.na(md5)) "missing" else "outside data/ and metadata/"
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
    p$inputs <- input_files(meta, catalogue)
  }

  # One row per check: "ok", "warning" (check before trusting) or "error" (can't be trusted as rendered)
  inputs_ok <- p$inputs$check %in% c("matches catalogue", "in git")
  p$checks <- data.frame(
    check = c("Code", "R&nbsp;packages", "Pipeline", "Input&nbsp;data"),
    level = c(
      if (length(p$uncommitted)) "warning" else "ok",
      if (isTRUE(p$renv_synced)) "ok" else "warning",
      if (p$built_by_targets && !is.na(p$report)) "ok" else "warning",
      if (any(p$inputs$check %in% c("DIFFERS from catalogue", "missing"))) "error"
        else if (!all(inputs_ok)) "warning" else "ok"
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
  renv <- if (isTRUE(p$renv_synced)) "match `renv.lock`"
          else if (isFALSE(p$renv_synced)) "differ from `renv.lock`" else "unknown (renv not active)"
  pipeline <- if (is.na(p$report)) "not a target in `_targets.R`: no pipeline information"
              else if (!p$built_by_targets) "rendered outside the pipeline: upstream results may be out of date"
              else paste0("built by targets &middot; upstream fingerprint `", p$fingerprint, "`")
  bad <- p$inputs[!p$inputs$check %in% c("matches catalogue", "in git"), , drop = FALSE]
  inputs <- sprintf("%d of %d file(s) verified (details below)", nrow(p$inputs) - nrow(bad), nrow(p$inputs))
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
