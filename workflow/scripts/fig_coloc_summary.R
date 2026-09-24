# Genes with a colocalised QTL per GWAS, stacked by group, one panel per set.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
master <- read_tsv(snakemake@input[["master"]])
gwas <- unlist(snakemake@params[["gwas"]])
if (nrow(master) == 0) {
  empty_plot(snakemake@output[[1]], "No colocalised QTLs")
} else {
  p <- plot_coloc_summary(master, gt, gwas, snakemake@params[["pp_h4"]])
  save_plot(p, snakemake@output[[1]], width = 4 + 3 * min(3, uniqueN(master$set)), height = max(4, 0.4 * length(gwas) + 2))
}
