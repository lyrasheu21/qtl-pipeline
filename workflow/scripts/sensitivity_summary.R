# FDR cciQTL counts for the main run and every sensitivity variant.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

info <- snakemake@params[["variant_info"]]
tables <- lapply(snakemake@input, fread)
variant_names <- unlist(snakemake@params[["variants"]])

directions <- rep(DIRECTIONS, length.out = length(tables))  # inputs are variant x direction
counts <- rbindlist(Map(function(dt, v, d) {
  # variants whose BED has no CCIs (e.g. a strict completeness filter) report zeros
  data.table(variant = v, direction = d, n_tested = nrow(dt),
             n_fdr = if (nrow(dt)) sum(dt$fdr_pass) else 0L,
             n_fdr_genes = if (nrow(dt)) uniqueN(dt$gene[dt$fdr_pass]) else 0L)
}, tables, variant_names, directions))
counts[, `:=`(
  completeness = vapply(variant, function(v) as.character(info[[v]]$completeness), ""),
  covariate_set = vapply(variant, function(v) as.character(info[[v]]$covset), ""),
  donor_fraction = vapply(variant, function(v) if (is.null(info[[v]]$donor_fraction)) NA_real_ else as.numeric(info[[v]]$donor_fraction), 0),
  n_permutations = vapply(variant, function(v) as.numeric(info[[v]]$n_perm), 0)
)]
setcolorder(counts, c("variant", "direction"))
write_tsv(counts, snakemake@output[["table"]])

counts[, variant := factor(variant, levels = unique(variant_names))]
p <- ggplot(counts, aes(variant, n_fdr, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = DIRECTION_COLORS) +
  labs(x = NULL, y = "FDR-significant cciQTLs") +
  theme_pipeline() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, snakemake@output[["plot"]], width = 9, height = 5)
