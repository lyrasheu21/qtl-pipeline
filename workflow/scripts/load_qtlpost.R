# Load the qtlpost package source for a pipeline script. Sourcing (instead of
# requiring an installed package) keeps the pipeline self-contained.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
