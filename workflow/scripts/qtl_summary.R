# QC summaries of the main permutation run:
#   - FDR counts per anchor subclass and direction
#   - dot plot of FDR cciQTLs vs CCIs tested per subclass
#   - QQ plots of permutation p-values per subclass, one page per superclass
#   - CCIs whose sending and receiving QTLs are both FDR-significant with the
#     same lead SNP (Fisher-combined p)

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

ct <- read_cell_types(snakemake@input[["cell_types"]])
qt <- rbind(fread(snakemake@input[["sending"]]), fread(snakemake@input[["receiving"]]))
qt[, direction := factor(direction, levels = DIRECTIONS)]
sub_levels <- union(subclass_order(ct), unique(qt$subclass))

# FDR counts -------------------------------------------------------------------
counts <- qt[, .(n_tested = .N, n_fdr = sum(fdr_pass), n_fdr_genes = uniqueN(gene[fdr_pass])),
             by = .(direction, superclass, class, subclass)]
counts[, subclass := factor(subclass, levels = sub_levels)]
setorder(counts, direction, subclass)
write_tsv(counts, snakemake@output[["counts"]])

p_dot <- ggplot(counts, aes(subclass, n_fdr, colour = subclass, size = n_tested)) +
  geom_point() +
  facet_wrap(~direction, ncol = 1) +
  scale_colour_manual(values = subclass_colors(ct), guide = "none") +
  scale_size(range = c(0.5, 4), name = "CCIs tested") +
  labs(x = "Anchor subclass", y = "FDR-significant cciQTLs") +
  theme_pipeline() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_dot, snakemake@output[["dotplot"]], width = 9, height = 7)

# QQ per subclass, one page per superclass -------------------------------------
qq <- qt[, qq_points(adj_beta_pval), by = .(direction, superclass, subclass)]
pdf(snakemake@output[["qq"]], width = 9, height = 7)
if (nrow(qq) == 0) {
  plot.new(); text(0.5, 0.5, "No permutation results")
}
for (sc in intersect(union(superclass_order(ct), unique(qq$superclass)), unique(qq$superclass))) {
  p <- ggplot(qq[superclass == sc], aes(expected, observed, colour = direction)) +
    geom_abline(linetype = "dashed") +
    geom_point(size = 0.8, alpha = 0.6) +
    facet_wrap(~subclass, scales = "free") +
    scale_colour_manual(values = DIRECTION_COLORS) +
    labs(title = sc, x = "Expected -log10(p)", y = "Observed -log10(p)") +
    theme_pipeline()
  print(p)
}
invisible(dev.off())

# Shared lead SNP on both sides -------------------------------------------------
send <- qt[direction == "sending" & fdr_pass == TRUE,
           .(phe_id, sending_snp = var_id, sending_slope = slope, sending_p = adj_beta_pval, sending_q = qval)]
recv <- qt[direction == "receiving" & fdr_pass == TRUE,
           .(phe_id, receiving_snp = var_id, receiving_slope = slope, receiving_p = adj_beta_pval, receiving_q = qval)]
shared <- merge(send, recv, by = "phe_id")[sending_snp == receiving_snp]
shared[, combined_p := fisher_combine(sending_p, receiving_p)]
shared <- cbind(shared, parse_cci_ids(shared$phe_id)[, !"cci"])
setorder(shared, combined_p)
message(sprintf("CCIs FDR-significant on both sides with the same lead SNP: %d", nrow(shared)))
write_tsv(shared, snakemake@output[["shared"]])
