# FDR QTL counts for the main run and every sensitivity variant, plus the
# hits-vs-K curve for a hidden-factor sweep.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

info <- snakemake@params[["variant_info"]]
covsets <- snakemake@params[["covsets"]]
variants <- unlist(snakemake@params[["variants"]])
sets <- unlist(snakemake@params[["sets"]])
tables <- lapply(unlist(snakemake@input), read_tsv)

counts <- rbindlist(Map(function(dt, v, s) {
  # variants whose BED has no phenotypes (e.g. a strict nonzero filter) report zeros
  cs <- covsets[[info[[v]]$covset]]
  data.table(variant = v, set = s, n_tested = nrow(dt),
             n_fdr_01 = if (nrow(dt)) sum(dt$qval <= 0.01, na.rm = TRUE) else 0L,
             n_fdr_05 = if (nrow(dt)) sum(dt$fdr_pass) else 0L,
             n_fdr_10 = if (nrow(dt)) sum(dt$qval <= 0.10, na.rm = TRUE) else 0L,
             pi0 = if (nrow(dt)) pi0_estimate(dt$adj_beta_pval) else NA_real_,
             min_nonzero_fraction = as.character(info[[v]]$fraction %||% "set default"),
             donor_fraction = if (is.null(info[[v]]$donor_fraction)) NA_real_ else as.numeric(info[[v]]$donor_fraction),
             covariate_set = info[[v]]$covset, k = as.integer(cs$hidden_k),
             n_permutations = as.numeric(info[[v]]$n_perm))
}, tables, variants, sets))
write_tsv(counts, snakemake@output[["table"]])
print(counts)

pdf(snakemake@output[["plot"]], width = 9, height = 5)
counts[, variant := factor(variant, levels = unique(variants))]
print(ggplot(counts, aes(variant, n_fdr_05, fill = set)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = set_colors(counts$set)) +
        labs(x = NULL, y = "QTLs at FDR 5%") +
        theme_pipeline() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
sweep <- counts[variant == "main" | grepl("^k[0-9]+$", variant)]
if (uniqueN(sweep$k) > 1) {
  sweep <- sweep[, .SD[1], by = .(set, k)]
  best <- sweep[, .SD[which.max(n_fdr_05)], by = set]
  message("Most QTLs at FDR 5% per set: ", paste(sprintf("%s K=%d (%d)", best$set, best$k, best$n_fdr_05), collapse = "; "))
  print(plot_ksweep(sweep[, .(set, k, n_fdr_01, n_fdr_05, n_fdr_10)]))
}
invisible(dev.off())
