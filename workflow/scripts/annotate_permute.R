# Standard permutation output for one set: QTLtools results joined to the
# phenotype annotation, with q-values computed per `fdr_group_by` (per group,
# class, superclass, or one global q-value with "none").
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

set_name <- snakemake@wildcards[["set"]]
gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
qt <- read_qtltools(snakemake@input[["permute"]], "permute")
annotation <- read_tsv(snakemake@input[["annotation"]])
message(sprintf("%s: %d phenotypes tested, %d with a cis variant", set_name, nrow(qt), sum(!is.na(qt$adj_beta_pval))))

out <- annotate_permute(qt, annotation, gt, set_name, snakemake@params[["fdr_group_by"]], snakemake@params[["fdr"]],
                        fdr_method = snakemake@params[["fdr_method"]])
message(sprintf("%s: FDR method %s, strata: %s", set_name, snakemake@params[["fdr_method"]], snakemake@params[["fdr_group_by"]]))
if (nrow(out)) {
  summ <- fdr_summary(out, by = "set")
  message(sprintf("%s: %d QTLs at q <= %g (pi0 = %.3f; %s)", set_name, summ$n_fdr, snakemake@params[["fdr"]],
                  summ$pi0, pi0_diagnostic(summ$pi0)))
} else {
  message(set_name, ": no phenotypes with results")
}
write_tsv(out, snakemake@output[["table"]])
# keys as they appear in column 1 of the nominal output (the group ID for --grp-best runs)
keys <- out[fdr_pass == TRUE, fifelse(is.na(grp_id), phenotype_id, grp_id)]
writeLines(unique(keys), snakemake@output[["fdr_keys"]])
