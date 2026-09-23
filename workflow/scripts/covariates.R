# Write a QTLtools covariate file for one covariate set: metadata covariates
# (with their transforms) plus the first n ancestry PCs, one column per donor.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

donors <- read_ids(snakemake@input[["donors"]])
meta <- read_table_any(snakemake@input[["metadata"]])
n_pcs <- as.integer(snakemake@params[["n_pcs"]])

pcs <- data.table()
pcs_id <- snakemake@params[["pcs_id"]]
if (n_pcs > 0) {
  if (length(snakemake@input[["ancestry_pcs"]]) == 0) stop("n_ancestry_pcs > 0 but inputs.ancestry_pcs is empty")
  pcs <- read_table_any(snakemake@input[["ancestry_pcs"]])
}

cov <- build_covariates(
  meta = meta, meta_id = snakemake@params[["metadata_id"]],
  pcs = pcs, pc_id = pcs_id, donors = donors,
  definitions = snakemake@params[["definitions"]],
  use = unlist(snakemake@params[["use"]]),
  n_pcs = n_pcs, pc_prefix = snakemake@params[["pc_prefix"]]
)
message(sprintf("Covariates: %s", paste(cov$id, collapse = ", ")))
write_tsv(cov, snakemake@output[[1]])
