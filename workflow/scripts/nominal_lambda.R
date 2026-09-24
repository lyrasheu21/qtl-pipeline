# Genomic inflation (lambda) and QQ plot per set from sampled nominal p-values.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

read_p <- function(path) {
  p <- suppressWarnings(as.numeric(readLines(path, warn = FALSE)))
  p[!is.na(p)]
}
sets <- unlist(snakemake@params[["sets"]])
p <- setNames(lapply(unlist(snakemake@input), read_p), sets)
lambda <- data.table(set = sets, n_sampled = vapply(p, length, integer(1)),
                     lambda = vapply(p, function(x) if (length(x)) genomic_lambda(x) else NA_real_, numeric(1)))
print(lambda)
write_tsv(lambda, snakemake@output[["table"]])

qq <- rbindlist(lapply(sets, function(s) qq_points(p[[s]])[, set := s]))
if (nrow(qq) == 0) {
  empty_plot(snakemake@output[["qq"]], "No nominal p-values sampled")
} else {
  labels <- setNames(sprintf("%s (lambda = %.3f)", lambda$set, lambda$lambda), lambda$set)
  p_qq <- ggplot(qq, aes(expected, observed)) +
    geom_abline(linetype = "dashed", colour = "red") +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~set, labeller = labeller(set = labels)) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = "Nominal QTL p-values (random sample)") +
    theme_pipeline()
  save_plot(p_qq, snakemake@output[["qq"]], width = 4 + 3 * min(3, length(sets)), height = 4 * ceiling(length(sets) / 3))
}
