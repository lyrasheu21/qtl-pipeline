# Combine per-subclass CCI-expression correlations and plot their distribution.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

res <- rbindlist(lapply(unlist(snakemake@input), fread), fill = TRUE)
write_tsv(res, snakemake@output[["table"]])

res <- res[!is.na(r)]
if (nrow(res) == 0) {
  empty_plot(snakemake@output[["plot"]], "No CCI-expression correlations computed")
} else {
  med <- res[, .(median_r = median(r)), by = side]
  print(med)
  p <- ggplot(res, aes(r)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    facet_wrap(~side, labeller = labeller(side = c(ligand = "Ligand expression (sending subclass)",
                                                   receptor = "Receptor expression (receiving subclass)"))) +
    labs(x = "Pearson r (CCI score vs gene expression, across donors)", y = "CCIs") +
    theme_pipeline()
  save_plot(p, snakemake@output[["plot"]], width = 9, height = 4)
}
