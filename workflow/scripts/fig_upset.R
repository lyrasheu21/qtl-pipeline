# UpSet plots of QTL genes shared across classes: all tested (page 1) and
# FDR-significant (page 2).
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
qt <- rbindlist(lapply(unlist(snakemake@input[["tables"]]), fread), fill = TRUE)
pdf(snakemake@output[[1]], width = 9, height = 6, onefile = TRUE)
plot_upset_classes(qt, gt, "All tested genes")
plot_upset_classes(qt[fdr_pass == TRUE], gt, "FDR QTL genes")
invisible(dev.off())
