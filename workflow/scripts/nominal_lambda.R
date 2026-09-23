# Genomic inflation (lambda) and QQ plot from sampled nominal p-values.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

read_p <- function(path) {
  p <- suppressWarnings(as.numeric(readLines(path, warn = FALSE)))
  p[!is.na(p)]
}
p <- list(sending = read_p(snakemake@input[["sending"]]),
          receiving = read_p(snakemake@input[["receiving"]]))

lambda <- data.table(direction = names(p),
                     n_sampled = vapply(p, length, integer(1)),
                     lambda = vapply(p, function(x) if (length(x)) genomic_lambda(x) else NA_real_, numeric(1)))
print(lambda)
write_tsv(lambda, snakemake@output[["table"]])

qq <- rbindlist(lapply(names(p), function(d) qq_points(p[[d]])[, direction := d]))
if (nrow(qq) == 0) {
  empty_plot(snakemake@output[["qq"]], "No nominal p-values sampled")
} else {
  labels <- setNames(sprintf("%s (lambda = %.3f)", lambda$direction, lambda$lambda), lambda$direction)
  p_qq <- ggplot(qq, aes(expected, observed)) +
    geom_abline(linetype = "dashed", colour = "red") +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~direction, labeller = labeller(direction = labels)) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = "Nominal cciQTL p-values (random sample)") +
    theme_pipeline()
  save_plot(p_qq, snakemake@output[["qq"]], width = 8, height = 4)
}
