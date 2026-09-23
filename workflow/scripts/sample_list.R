# Donors used for QTL mapping: those present in the QC'd genotypes, the CCI
# matrix and the metadata. Order follows the CCI matrix columns.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

vcf <- read_ids(snakemake@input[["vcf_samples"]])
cci <- colnames(readRDS(snakemake@input[["cci"]]))
meta <- read_table_any(snakemake@input[["metadata"]])
meta_id <- snakemake@params[["metadata_id"]]
if (!meta_id %in% names(meta)) stop("Metadata ID column '", meta_id, "' not found")
meta_ids <- as.character(meta[[meta_id]])

donors <- cci[cci %in% vcf & cci %in% meta_ids]

report <- data.table(
  source = c("genotype_vcf", "cci_matrix", "metadata", "all_three"),
  n_donors = c(length(vcf), length(cci), length(meta_ids), length(donors))
)
print(report)
write_tsv(report, snakemake@output[["report"]])

if (length(donors) < snakemake@params[["min_donors"]]) {
  stop(sprintf("Only %d donors are shared by genotypes, CCI matrix and metadata (min_donors = %d). ",
               length(donors), snakemake@params[["min_donors"]]),
       "Check that donor IDs match across files (examples - VCF: ",
       paste(head(vcf, 3), collapse = ", "), "; CCI: ", paste(head(cci, 3), collapse = ", "), ").")
}
writeLines(donors, snakemake@output[["donors"]])
