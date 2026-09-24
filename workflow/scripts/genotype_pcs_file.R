# Genotype PCs from a user file (samples as rows or PCs as rows, e.g. a .pca
# file), rewritten as a samples x PCs table for the covariate step.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

samples <- read_ids(snakemake@input[["samples"]])
pcs <- read_sample_factors(snakemake@input[["pcs"]], samples)
colnames(pcs) <- sub("^.*_?(PC[0-9]+)$", "\\1", colnames(pcs))
found <- intersect(samples, rownames(pcs))
message(sprintf("Genotype PCs: %d PCs for %d of %d genotyped samples", ncol(pcs), length(found), length(samples)))
if (!length(found)) stop("No genotyped sample IDs found in ", snakemake@input[["pcs"]])
write_tsv(data.table(sample = found, pcs[found, , drop = FALSE]), snakemake@output[[1]])
