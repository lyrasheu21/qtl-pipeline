# Combine MAGMA results across GWAS and draw the heatmap.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
gwas_names <- unlist(snakemake@params[["gwas"]])
res <- rbindlist(Map(function(path, g) read_magma_gsa(path)[, gwas := g], unlist(snakemake@input[["results"]]), gwas_names))
write_tsv(res, snakemake@output[["table"]])
if (nrow(res) == 0) {
  empty_plot(snakemake@output[["plot"]], "No MAGMA results (empty gene sets)")
} else {
  level <- snakemake@wildcards[["level"]]
  order <- switch(level, group = group_order(gt), class = class_order(gt), superclass = superclass_order(gt))
  p <- plot_magma_heatmap(res, order, title = sprintf("MAGMA: %s QTL genes by %s", snakemake@wildcards[["threshold"]], level))
  save_plot(p, snakemake@output[["plot"]], width = max(6, 0.3 * uniqueN(res$set) + 4), height = 3 + 0.8 * length(gwas_names))
}
