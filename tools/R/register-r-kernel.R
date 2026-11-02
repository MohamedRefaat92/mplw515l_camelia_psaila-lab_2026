# Writes the project's R kernelspec into env/share/jupyter/kernels/<name>-r.
# Run from the project root with the env's R: env/bin/Rscript tools/R/register-r-kernel.R <name>
name <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(name)) stop("usage: Rscript tools/R/register-r-kernel.R <project-name>")
root   <- normalizePath(".", winslash = "/")
prefix <- normalizePath(file.path(R.home(), "..", ".."), winslash = "/")  # conda env root
r_bin  <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")

kdir <- file.path(prefix, "share", "jupyter", "kernels", paste0(name, "-r"))
dir.create(kdir, recursive = TRUE, showWarnings = FALSE)

spec <- list(
  argv = c(r_bin, "--slave", "-e", "IRkernel::main()", "--args", "{connection_file}"),
  display_name = sprintf("R (%s · renv)", name),
  language = "R",
  env = list(
    RENV_PROJECT   = root,
    R_PROFILE_USER = file.path(root, ".Rprofile")
  )
)
jsonlite::write_json(spec, file.path(kdir, "kernel.json"), auto_unbox = TRUE, pretty = TRUE)
invisible(file.copy(list.files(system.file("kernelspec", package = "IRkernel"), pattern = "^logo", full.names = TRUE), kdir, overwrite = TRUE))
message("Wrote ", file.path(kdir, "kernel.json"))

