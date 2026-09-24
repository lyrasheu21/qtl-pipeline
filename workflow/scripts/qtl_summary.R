# QC summaries of the main permutation run across all sets:
#   FDR summary per set and group (pi0, counts at FDR 1/5/10%), dot plot,
#   QQ plots per group, histograms of permutation p-values, and phenotypes
#   whose lead SNP is FDR-significant in more than one set.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

gt <- read_group_table(if (length(snakemake@input[["group_table"]])) snakemake@input[["group_table"]] else NULL)
qt <- rbindlist(lapply(unlist(snakemake@input[["tables"]]), fread), fill = TRUE)

by_group <- fdr_summary(qt, by = c("set", "group"))
by_set <- fdr_summary(qt, by = "set")[, group := "all groups"]
summary <- rbind(by_set, by_group, use.names = TRUE)
summary[, pi0_note := vapply(pi0, pi0_diagnostic, "")]
setcolorder(summary, c("set", "group"))
write_tsv(summary, snakemake@output[["summary"]])
print(by_set)

save_plot(plot_fdr_counts(by_group, gt), snakemake@output[["dotplot"]], width = 9,
          height = 2.5 + 2.5 * uniqueN(by_group$set))

pdf(snakemake@output[["qq"]], width = 9, height = 7)
pages <- if (nrow(qt)) plot_qq_pages(qt, gt) else list()
if (!length(pages)) { plot.new(); text(0.5, 0.5, "No permutation results") }
for (p in pages) print(p)
invisible(dev.off())

save_plot(plot_pvalue_hist(qt), snakemake@output[["hist"]], width = 9, height = 2 + 2 * uniqueN(qt$set))

hits <- qt[fdr_pass == TRUE, .(set, phenotype_id, var_id, adj_beta_pval, qval)]
shared <- hits[, if (uniqueN(set) > 1 && uniqueN(var_id) == 1)
  .(sets = paste(set, collapse = ";"), lead_snp = var_id[1], combined_p = if (.N == 2) fisher_combine(adj_beta_pval[1], adj_beta_pval[2]) else NA_real_),
  by = phenotype_id]
if (!nrow(shared)) shared <- data.table(phenotype_id = character(), sets = character(), lead_snp = character(), combined_p = numeric())
setorder(shared, combined_p)
message(sprintf("Phenotypes FDR-significant in several sets with the same lead SNP: %d", nrow(shared)))
write_tsv(shared, snakemake@output[["shared"]])
