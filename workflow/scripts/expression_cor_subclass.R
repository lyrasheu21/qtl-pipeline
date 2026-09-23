# For one subclass: Pearson correlation across donors between each CCI score
# and the expression of its ligand (CCIs sent by this subclass) or its receptor
# (CCIs received by this subclass).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

sc <- snakemake@wildcards[["subclass"]]
path <- snakemake@params[["expression"]]
min_donors <- snakemake@params[["min_donors"]]
empty <- data.table(cci = character(), side = character(), subclass = character(), gene = character(),
                    n_donors = integer(), r = numeric())

if (!file.exists(path)) {
  message("No expression file for ", sc, " (", path, "); writing an empty table")
  write_tsv(empty, snakemake@output[[1]])
  quit(save = "no")
}

mat <- readRDS(snakemake@input[["cci"]])
expr <- read_expression(path)
donors <- intersect(colnames(mat), colnames(expr))
ids <- parse_cci_ids(rownames(mat))

jobs <- rbind(
  ids[sending_subclass == sc, .(cci, side = "ligand", gene = ligand)],
  ids[receiving_subclass == sc, .(cci, side = "receptor", gene = receptor)]
)
jobs <- jobs[gene %in% rownames(expr)]
message(sprintf("%s: %d donors shared, %d CCI-gene pairs", sc, length(donors), nrow(jobs)))

if (length(donors) < min_donors || nrow(jobs) == 0) {
  write_tsv(empty, snakemake@output[[1]])
  quit(save = "no")
}

cci_vals <- mat[jobs$cci, donors, drop = FALSE]
gene_vals <- expr[jobs$gene, donors, drop = FALSE]
ok <- !is.na(cci_vals) & !is.na(gene_vals)
n_ok <- rowSums(ok)
r <- vapply(seq_len(nrow(jobs)), function(i) {
  if (n_ok[i] < min_donors) return(NA_real_)
  suppressWarnings(cor(cci_vals[i, ok[i, ]], gene_vals[i, ok[i, ]]))
}, numeric(1))

write_tsv(jobs[, .(cci, side, subclass = sc, gene, n_donors = n_ok, r = r)], snakemake@output[[1]])
