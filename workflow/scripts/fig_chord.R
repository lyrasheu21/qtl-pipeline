# Mediation chord diagram for one anchor ("group:gene").
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
plot_chord(fread(snakemake@input[["mediation"]]), snakemake@params[["target"]], gt, snakemake@output[[1]])
