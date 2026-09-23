# Build QTLtools phenotype BEDs from the CCI matrix.
#
# Each CCI is written twice: once at its ligand gene (sending) and once at its
# receptor gene (receiving). CCIs below the completeness threshold, with a gene
# missing from the annotation (e.g. multi-subunit complexes) or on a chromosome
# not being analysed are dropped; the report counts each reason.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

mat <- readRDS(snakemake@input[["cci"]])
donors <- read_ids(snakemake@input[["donors"]])
annot <- read_annotation(snakemake@input[["annotation"]])
chromosomes <- as.character(unlist(snakemake@params[["chromosomes"]]))
min_completeness <- as.numeric(snakemake@params[["min_completeness"]])

mat <- mat[, donors, drop = FALSE]
completeness <- rowMeans(!is.na(mat) & mat != 0)
if (isTRUE(snakemake@params[["zero_as_missing"]])) mat[mat == 0] <- NA
passes_completeness <- completeness >= min_completeness

ids <- parse_cci_ids(rownames(mat))
report <- list()

for (direction in DIRECTIONS) {
  gene <- if (direction == "sending") ids$ligand else ids$receptor
  coords <- annot[match(gene, gene_name)]
  annotated <- !is.na(coords$gene_name)
  on_chrom <- annotated & coords$chrom_norm %in% chromosomes
  keep <- passes_completeness & on_chrom

  bed <- data.table(
    chr = coords$chrom_norm[keep],
    start = as.integer(coords$start[keep]),
    end = as.integer(coords$end[keep]),
    pid = rownames(mat)[keep],
    gid = rownames(mat)[keep],
    strand = "+"
  )
  bed <- cbind(bed, as.data.table(mat[keep, , drop = FALSE]))
  bed <- bed[order(as.integer(factor(chr, levels = chromosomes)), start)]
  setnames(bed, "chr", "#chr")
  fwrite(bed, snakemake@output[[direction]], sep = "\t", na = "NA", quote = FALSE)

  report[[direction]] <- data.table(
    direction = direction,
    n_cci_in_matrix = nrow(mat),
    below_completeness = sum(!passes_completeness),
    gene_not_in_annotation = sum(passes_completeness & !annotated),
    gene_not_on_analysed_chromosome = sum(passes_completeness & annotated & !on_chrom),
    n_cci_in_bed = sum(keep),
    n_donors = ncol(mat),
    min_completeness = min_completeness
  )
}

report <- rbindlist(report)
print(report)
write_tsv(report, snakemake@output[["report"]])
