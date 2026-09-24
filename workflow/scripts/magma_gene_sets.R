# MAGMA --set-annot file (gene ID, set) for one threshold and grouping level.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

qt <- rbindlist(lapply(unlist(snakemake@input[["tables"]]), fread), fill = TRUE)
sets <- magma_gene_sets(qt, read_annotation(snakemake@input[["annotation"]]),
                        threshold = snakemake@wildcards[["threshold"]], level = snakemake@wildcards[["level"]],
                        id_column = snakemake@params[["gene_id_column"]], nominal_p = snakemake@params[["nominal_p"]])
print(sets[, .N, by = set])
fwrite(sets, snakemake@output[[1]], sep = "\t", col.names = FALSE)
