# Number of unique anchor genes with a colocalised cciQTL per GWAS, stacked by
# anchor subclass and split by direction.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

ct <- read_cell_types(snakemake@input[["cell_types"]])
master <- fread(snakemake@input[["master"]])
gwas_levels <- unlist(snakemake@params[["gwas"]])

if (nrow(master) == 0) {
  empty_plot(snakemake@output[[1]], "No colocalised cciQTLs")
} else {
  counts <- master[, .(n_genes = uniqueN(gene)), by = .(gwas, direction, subclass)]
  counts[, `:=`(gwas = factor(gwas, levels = rev(gwas_levels)),
                direction = factor(direction, levels = DIRECTIONS),
                subclass = factor(subclass, levels = union(subclass_order(ct), unique(subclass))))]
  p <- ggplot(counts, aes(gwas, n_genes, fill = subclass)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~direction) +
    scale_fill_manual(values = subclass_colors(ct), drop = TRUE) +
    labs(x = "GWAS", y = sprintf("Anchor genes with colocalised cciQTL (PP.H4 >= %g)", snakemake@params[["pp_h4"]]),
         fill = "Anchor subclass") +
    theme_pipeline()
  save_plot(p, snakemake@output[[1]], width = 10, height = max(4, 0.4 * length(gwas_levels) + 2))
}
