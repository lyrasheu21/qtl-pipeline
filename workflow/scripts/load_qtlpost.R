# Load the qtlpost package source for a pipeline script. Sourcing (instead of
# requiring an installed package) keeps the pipeline self-contained.

# Use only the conda environment's R library. Personal libraries added through
# R_LIBS / R_LIBS_USER (often built for another R version) would otherwise be
# searched first and can shadow the environment's packages with builds that
# fail to load. .libPaths() always keeps the site and base libraries.
.libPaths(.Library)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
