# MAGMA --set-annot file: one line per (gene ID, set). Sets are
# "<direction>_<group>", where group is the anchor subclass, class or
# superclass. Genes come from FDR-significant cciQTLs (threshold = fdr) or from
# cciQTLs with adj_beta_pval <= nominal_p (threshold = nominal).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

threshold <- snakemake@wildcards[["threshold"]]
level <- snakemake@wildcards[["level"]]
id_col <- snakemake@params[["gene_id_column"]]

qt <- rbind(fread(snakemake@input[["sending"]]), fread(snakemake@input[["receiving"]]))
hits <- if (threshold == "fdr") qt[fdr_pass == TRUE] else qt[adj_beta_pval <= snakemake@params[["nominal_p"]]]

annot <- read_annotation(snakemake@input[["annotation"]])
if (!id_col %in% names(annot)) stop("Annotation has no column '", id_col, "'")
hits[, gene_id := annot[[id_col]][match(gene, annot$gene_name)]]
message(sprintf("%d of %d genes mapped to %s", uniqueN(hits[!is.na(gene_id), gene]), uniqueN(hits$gene), id_col))

sets <- unique(hits[!is.na(gene_id), .(gene_id, set = paste0(direction, "_", get(level)))])
setorder(sets, set, gene_id)
print(sets[, .N, by = set])
fwrite(sets, snakemake@output[[1]], sep = "\t", col.names = FALSE)
