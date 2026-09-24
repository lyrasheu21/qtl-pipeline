# Combine phenotype-expression correlations across sets and plot them.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

res <- rbindlist(lapply(unlist(snakemake@input), read_tsv), fill = TRUE)
write_tsv(res, snakemake@output[["table"]])
res <- res[!is.na(r)]
if (nrow(res) == 0) {
  empty_plot(snakemake@output[["plot"]], "No phenotype-expression correlations computed")
} else {
  print(res[, .(median_r = median(r), n = .N), by = set])
  p <- ggplot(res, aes(r)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    facet_wrap(~set) +
    labs(x = "Pearson r (phenotype vs its gene's expression, across donors)", y = "Phenotypes") +
    theme_pipeline()
  save_plot(p, snakemake@output[["plot"]], width = 3 + 3 * min(3, uniqueN(res$set)), height = 4)
}
