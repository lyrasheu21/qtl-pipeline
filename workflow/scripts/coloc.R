# coloc.abf between each FDR-significant cciQTL and one GWAS, with tiering of
# the hits and locus plots.
#
# Outputs
#   summary: one row per CCI tested, with locus-level posteriors; for hits
#            (PP.H4 >= pp_h4) also the tier and lead-SNP GWAS / QTL statistics
#   snps:    SNP-level results for hits (SNP.PP.H4 > snp_pp_h4)
#   plots:   locus plots (QTL, GWAS, SNP.PP.H4, genes) for Gold_Unique hits

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
suppressPackageStartupMessages({
  library(coloc)
  library(patchwork)
})

direction <- snakemake@wildcards[["direction"]]
spec <- snakemake@params[["spec"]]
settings <- snakemake@params[["settings"]]
n_qtl <- length(read_ids(snakemake@input[["donors"]]))
plot_dir <- snakemake@output[["plots"]]
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

nom <- read_qtltools(snakemake@input[["nominal"]], "nominal")
ccis <- unique(nom$phe_id)
message(sprintf("%s / %s: %d FDR-significant CCIs, %d nominal rows", spec$name, direction, length(ccis), nrow(nom)))

gwas_cols <- unlist(spec[c("snp_col", "bp_col", "beta_col", "se_col", "p_col", "n_col")])
gwas_cols <- gwas_cols[!is.na(gwas_cols) & nzchar(gwas_cols)]
gwas <- fread(snakemake@input[["gwas"]], select = unname(gwas_cols))
gwas <- standardize_gwas(gwas[get(spec$snp_col) %chin% unique(nom$var_id)], spec)
message(sprintf("GWAS SNPs overlapping cciQTL SNPs: %d", nrow(gwas)))

annot <- read_annotation(snakemake@input[["annotation"]])
ids <- parse_cci_ids(ccis)
focus_genes <- setNames(if (direction == "sending") ids$ligand else ids$receptor, ccis)

summaries <- list()
snp_rows <- list()
for (cci in ccis) {
  d <- nom[phe_id == cci]
  shared <- intersect(d$var_id, gwas$snp)
  if (length(shared) < settings$min_snps) next
  d <- d[var_id %chin% shared]
  g <- gwas[snp %chin% shared]

  res <- tryCatch(
    suppressMessages(suppressWarnings(coloc.abf(
      qtl_coloc_dataset(d, n_qtl), gwas_coloc_dataset(g, spec),
      p1 = settings$p1, p2 = settings$p2, p12 = settings$p12))),
    error = function(e) { message("coloc failed for ", cci, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next

  s <- as.list(res$summary)
  row <- data.table(gwas = spec$name, direction = direction, cci = cci,
                    n_snps = s$nsnps, PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
                    PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, ratio_H4_H3 = s$PP.H4.abf / s$PP.H3.abf,
                    tier = NA_character_, lead_snp = NA_character_, lead_snp_pp_h4 = NA_real_,
                    gwas_beta = NA_real_, gwas_p = NA_real_,
                    qtl_beta = NA_real_, qtl_p = NA_real_, qtl_is_best_hit = NA)

  if (s$PP.H4.abf >= settings$pp_h4) {
    r <- as.data.table(res$results)
    if (!"position" %in% names(r)) r[, position := d$var_from[match(snp, d$var_id)]]
    kept <- r[SNP.PP.H4 > settings$snp_pp_h4]
    lead <- r[which.max(SNP.PP.H4)]
    g_lead <- g[snp == lead$snp][1]
    d_lead <- d[var_id == lead$snp][1]
    row[, `:=`(tier = classify_coloc_tier(kept$SNP.PP.H4),
               lead_snp = lead$snp, lead_snp_pp_h4 = lead$SNP.PP.H4,
               gwas_beta = g_lead$beta, gwas_p = g_lead$p,
               qtl_beta = d_lead$slope, qtl_p = d_lead$nom_pval,
               qtl_is_best_hit = d_lead$best_hit == 1)]
    if (nrow(kept)) {
      snp_rows[[cci]] <- data.table(gwas = spec$name, direction = direction, cci = cci,
                                    snp = kept$snp, position = kept$position, SNP.PP.H4 = kept$SNP.PP.H4)
    }

    n_plotted <- length(list.files(plot_dir))
    if (row$tier == "Gold_Unique" && n_plotted < settings$max_locus_plots) {
      focus_gene <- focus_genes[[cci]]
      window <- range(d$var_from)
      pts <- rbind(
        d[, .(position = var_from, snp = var_id, value = -log10(nom_pval), track = "cciQTL -log10(p)")],
        g[, .(position = bp, snp, value = -log10(p), track = "GWAS -log10(p)")],
        r[, .(position, snp, value = SNP.PP.H4, track = "SNP.PP.H4")]
      )
      pts[, `:=`(track = factor(track, levels = c("cciQTL -log10(p)", "GWAS -log10(p)", "SNP.PP.H4")),
                 is_lead = snp == lead$snp)]
      p_assoc <- ggplot(pts, aes(position, value)) +
        geom_point(aes(colour = is_lead, size = is_lead)) +
        facet_wrap(~track, ncol = 1, scales = "free_y", strip.position = "left") +
        scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "red"), guide = "none") +
        scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2.5), guide = "none") +
        coord_cartesian(xlim = window) +
        labs(x = NULL, y = NULL,
             title = cci, subtitle = sprintf("%s, %s | PP.H4 = %.3f | lead SNP %s", spec$name, direction,
                                             s$PP.H4.abf, lead$snp)) +
        theme_pipeline() + theme(strip.placement = "outside")
      p_genes <- gene_track(annot, d$var_chr[1], window[1], window[2], highlight = focus_gene) +
        coord_cartesian(xlim = window) + labs(x = sprintf("chr%s position", d$var_chr[1]))
      fig <- p_assoc / p_genes + plot_layout(heights = c(3, 0.6))
      save_plot(fig, file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_.-]+", "_", sub("->", "__to__", cci, fixed = TRUE)), ".pdf")), width = 7, height = 9)
    }
  }
  summaries[[cci]] <- row
}

summary <- rbindlist(summaries)
if (nrow(summary) == 0) {
  summary <- data.table(gwas = character(), direction = character(), cci = character(), n_snps = numeric(),
                        PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(),
                        PP.H4 = numeric(), ratio_H4_H3 = numeric(), tier = character(), lead_snp = character(),
                        lead_snp_pp_h4 = numeric(), gwas_beta = numeric(), gwas_p = numeric(),
                        qtl_beta = numeric(), qtl_p = numeric(), qtl_is_best_hit = logical())
}
snps <- rbindlist(snp_rows)
if (nrow(snps) == 0) {
  snps <- data.table(gwas = character(), direction = character(), cci = character(),
                     snp = character(), position = numeric(), SNP.PP.H4 = numeric())
}
message(sprintf("CCIs tested: %d; colocalised (PP.H4 >= %g): %d",
                nrow(summary), settings$pp_h4, sum(summary$PP.H4 >= settings$pp_h4)))
write_tsv(summary, snakemake@output[["summary"]])
write_tsv(snps, snakemake@output[["snps"]])
