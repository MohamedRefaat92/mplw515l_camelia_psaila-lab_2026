# Example analysis functions: replace with your own. Each one is a step the pipeline (_targets.R) calls.

read_counts <- function(path) {
  as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
}

read_samples <- function(path) {
  read.csv(path)
}

#' log2 counts per million
log_cpm <- function(counts) {
  log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
}

#' Mean log-CPM per condition and the treated - control difference, one row per gene
summarise_genes <- function(logcpm, samples) {
  ctrl <- samples$sample[samples$condition == "control"]
  trt <- samples$sample[samples$condition == "treated"]
  out <- data.frame(
    gene    = rownames(logcpm),
    control = rowMeans(logcpm[, ctrl, drop = FALSE]),
    treated = rowMeans(logcpm[, trt, drop = FALSE])
  )
  out$log2_fc <- out$treated - out$control
  out[order(-abs(out$log2_fc)), ]
}

plot_genes <- function(summary) {
  ggplot2::ggplot(summary, ggplot2::aes(control, treated, colour = abs(log2_fc) > 1)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"),
                                 labels = c("|log2 FC| <= 1", "|log2 FC| > 1"), name = NULL) +
    ggplot2::labs(x = "control (mean log2 CPM)", y = "treated (mean log2 CPM)") +
    ggplot2::theme_minimal()
}
