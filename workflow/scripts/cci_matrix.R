# Read the user's CCI score matrix (parquet/TSV/CSV) once and store it as an
# RDS numeric matrix (rows = CCIs, columns = donors) for the later steps.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

mat <- read_cci_matrix(snakemake@input[[1]], snakemake@params[["id_column"]])
invisible(parse_cci_ids(rownames(mat)))  # fail early on malformed IDs
message(sprintf("CCI matrix: %d CCIs x %d donors", nrow(mat), ncol(mat)))
saveRDS(mat, snakemake@output[[1]])
