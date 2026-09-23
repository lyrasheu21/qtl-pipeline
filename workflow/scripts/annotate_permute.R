# Annotate QTLtools permutation results and compute q-values per anchor
# subclass.
#
# The anchor is the side the QTL was mapped at: the ligand's subclass/gene for
# sending, the receptor's for receiving. q-values use the beta-approximated
# permutation p-value (adj_beta_pval).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

direction <- snakemake@wildcards[["direction"]]
fdr <- snakemake@params[["fdr"]]
ct <- read_cell_types(snakemake@input[["cell_types"]])

qt <- read_qtltools(snakemake@input[["permute"]], "permute")
n_all <- nrow(qt)
qt <- qt[!is.na(adj_beta_pval)]
message(sprintf("%s: %d CCIs tested, %d with a cis variant", direction, n_all, nrow(qt)))

ann <- annotate_cci(qt$phe_id, direction, ct)
out <- cbind(
  qt[, .(phe_id, phe_chr, phe_from, phe_to, n_var_in_cis, var_id, var_chr, var_from,
         dist_phe_var, slope, slope_se, nom_pval, r_squared, adj_emp_pval, adj_beta_pval)],
  ann[, !"cci"]
)
out[, qval := qvalue_by_group(adj_beta_pval, subclass)]
out[, fdr_pass := !is.na(qval) & qval <= fdr]
setorder(out, qval)

message(sprintf("%s: %d CCIs at q <= %g", direction, sum(out$fdr_pass), fdr))
write_tsv(out, snakemake@output[["table"]])
writeLines(unique(out[fdr_pass == TRUE, phe_id]), snakemake@output[["fdr_ids"]])
