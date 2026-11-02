# This profile is only for the project's own R in ./env. Any other R (e.g. a
# system R found first on PATH) would create a second, incompatible renv
# library, so it gets a warning on stderr and a plain session instead.
if (local({
  root <- Sys.getenv("RENV_PROJECT", unset = getwd())
  ok <- identical(
    normalizePath(file.path(R.home(), "..", ".."), mustWork = FALSE),
    normalizePath(file.path(root, "env"), mustWork = FALSE)
  )
  if (!ok) message("[", basename(root), "] ", R.home(), " is not this project's R; renv not activated. Use env/bin/R.")
  ok
})) {

# Put the conda env's tools (compilers, pkg-config, pandoc) on PATH even when
# R is started directly as env/bin/R, without `mamba activate`.
local({
  prefix <- normalizePath(file.path(R.home(), "..", ".."), mustWork = FALSE)
  bins <- if (.Platform$OS.type == "windows") {
    file.path(prefix, c("Library/mingw-w64/bin", "Library/usr/bin", "Library/bin", "Scripts"))
  } else {
    file.path(prefix, "bin")
  }
  Sys.setenv(PATH = paste(c(bins, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  # conda's pkg-config searches $CONDA_PREFIX when it is set: with (base) active in the shell, packages
  # would compile against base's libraries, or not find ./env's (e.g. fs: "libuv not found")
  Sys.unsetenv("CONDA_PREFIX")

  # conda's quarto only works with the variables its activation script sets (QUARTO_SHARE_PATH, ...).
  # Apply just those, not a full activation, which would also change compiler settings.
  act <- file.path(prefix, "etc", "conda", "activate.d", "quarto.sh")
  if (file.exists(act)) {
    for (line in grep("^[[:space:]]*export [A-Z_]+=", readLines(act), value = TRUE)) {
      kv <- sub("^[[:space:]]*export ", "", line)
      do.call(Sys.setenv, stats::setNames(list(sub("^[^=]*=", "", kv)), sub("=.*$", "", kv)))
    }
  }
})

# conda R cannot use CRAN macOS binaries, so always build from source
options(
  pkgType = "source",
  repos = c(
    CRAN        = "https://cloud.r-project.org",
    vscdebugger = "https://manuelhentschel.r-universe.dev"
  )
)

# RENV_PROJECT lets kernels started in subfolders find the project
local({
  root <- Sys.getenv("RENV_PROJECT", unset = getwd())
  source(file.path(root, "renv", "activate.R"))
})

# Keep command history across sessions in .Rhistory (git-ignored). Saved after every
# command, so it survives closing the terminal; R's own save-at-exit is off (--no-save in ./.vscode/settings.json).

if (interactive()) local({
  hist <- file.path(Sys.getenv("RENV_PROJECT", unset = getwd()), ".Rhistory")
  invisible(addTaskCallback(function(...) {
    try(utils::savehistory(hist), silent = TRUE)
    TRUE                                   # TRUE = keep this callback for later commands
  }, name = "save-history"))
})

}
