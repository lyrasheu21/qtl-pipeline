# pi1 replication of FDR QTLs in an external eQTL resource, matched on
# (lead SNP, gene, group), plus a beta-vs-beta scatter.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

s <- snakemake@params[["settings"]]
qt <- rbindlist(lapply(unlist(snakemake@input[["tables"]]), fread), fill = TRUE)
hits <- qt[fdr_pass == TRUE & !is.na(gene), .(set, phenotype_id, group, gene, snp = var_id, chr = var_chr, qtl_beta = slope)]
message(sprintf("FDR QTLs to look up: %d", nrow(hits)))
res <- pi1_replication(hits, s$eqtl_pattern, has_header = isTRUE(s$has_header), columns = unlist(s$columns),
                       snp_column = s$snp_column, gene_column = s$gene_column, beta_column = s$beta_column,
                       z_column = s$z_column, p_column = s$p_column)
print(res$pi1)
write_tsv(res$pi1, snakemake@output[["pi1"]])
write_tsv(res$matched, snakemake@output[["matched"]])

plot_data <- res$matched[!is.na(eqtl_beta)]
if (nrow(plot_data) < 3) {
  empty_plot(snakemake@output[["scatter"]], "Fewer than 3 QTLs matched in the eQTL resource")
} else {
  rho <- plot_data[, .(rho = suppressWarnings(cor(eqtl_beta, qtl_beta, method = "spearman"))), by = set]
  labels <- setNames(sprintf("%s (Spearman rho = %.2f)", rho$set, rho$rho), rho$set)
  p <- ggplot(plot_data, aes(eqtl_beta, qtl_beta, colour = set)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    scale_colour_manual(values = set_colors(plot_data$set), labels = labels) +
    labs(x = "eQTL beta", y = "QTL beta") +
    theme_pipeline()
  save_plot(p, snakemake@output[["scatter"]], width = 7, height = 5)
}
