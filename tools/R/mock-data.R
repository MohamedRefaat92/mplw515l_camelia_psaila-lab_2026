# Small, reproducible mock inputs for trying the example pipeline (not analysis).
#   <dir>/counts.tsv       input DATA: genes x samples count matrix (first column: gene); linked into data/
#   metadata/samples.csv   METADATA: the sample sheet (sample, condition); written only if it doesn't exist yet
write_mock_data <- function(dir, samples_file = "metadata/samples.csv",
                            n_genes = 200, n_per_group = 3, seed = 1) {
  set.seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  samples <- data.frame(
    sample    = sprintf("S%02d", seq_len(2 * n_per_group)),
    condition = rep(c("control", "treated"), each = n_per_group)
  )
  base <- rgamma(n_genes, shape = 2, scale = 50)                  # baseline expression per gene
  effect <- ifelse(runif(n_genes) < 0.1, 2^rnorm(n_genes, 0, 1.5), 1)  # ~10% of genes respond
  counts <- sapply(seq_len(nrow(samples)), function(j) {
    mu <- base * if (samples$condition[j] == "treated") effect else 1
    rnbinom(n_genes, mu = mu, size = 10)
  })
  dimnames(counts) <- list(sprintf("gene%03d", seq_len(n_genes)), samples$sample)
  write.table(data.frame(gene = rownames(counts), counts, check.names = FALSE),
              file.path(dir, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!file.exists(samples_file)) {          # never overwrite a real sample sheet
    dir.create(dirname(samples_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(samples, samples_file, row.names = FALSE)
  }
  invisible(c(counts = file.path(dir, "counts.tsv"), samples = samples_file))
}
