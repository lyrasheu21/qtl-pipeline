# Heatmap of MAGMA gene-set results: -log10(p) per GWAS x cell group, one
# panel per direction. "*" = p < 0.05; "#" = BH-adjusted p < 0.05, adjusting
# across GWAS within each set (as in the original analysis).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

level <- snakemake@wildcards[["level"]]
ct <- read_cell_types(snakemake@params[["cell_types"]])
gwas_names <- unlist(snakemake@params[["gwas"]])

res <- rbindlist(Map(function(path, g) read_magma_gsa(path)[, gwas := g],
                     unlist(snakemake@input), gwas_names))
if (nrow(res)) {
  res[, direction := sub("_.*$", "", set)]
  res[, group := sub("^(sending|receiving)_", "", set)]
  res[, p_adj := p.adjust(p, method = "BH"), by = set]
  res[, label := fifelse(p_adj < 0.05, "#", fifelse(p < 0.05, "*", ""))]
}
write_tsv(res, snakemake@output[["table"]])

if (nrow(res) == 0) {
  empty_plot(snakemake@output[["plot"]], "No MAGMA results (empty gene sets)")
} else {
  group_levels <- switch(level, subclass = subclass_order(ct), class = class_order(ct),
                         superclass = superclass_order(ct))
  res[, group := factor(group, levels = union(group_levels, unique(group)))]
  res[, gwas := factor(gwas, levels = rev(gwas_names))]
  res[, direction := factor(direction, levels = DIRECTIONS)]
  p <- ggplot(res, aes(group, gwas, fill = -log10(p))) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3) +
    facet_wrap(~direction, ncol = 1) +
    scale_fill_gradient(low = "white", high = "red", name = "-log10(p)") +
    labs(x = sprintf("cciQTL genes by anchor %s", level), y = "GWAS",
         title = sprintf("MAGMA: %s cciQTL gene sets", snakemake@wildcards[["threshold"]])) +
    theme_pipeline() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, snakemake@output[["plot"]], width = max(6, 0.3 * nlevels(res$group) + 3),
            height = max(5, 0.5 * length(gwas_names) * 2 + 2))
}
