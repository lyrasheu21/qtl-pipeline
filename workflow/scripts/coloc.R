# coloc.abf between the FDR-significant QTLs of one set and one GWAS, with
# locus plots for Gold_Unique hits (qtlpost::run_coloc / plot_locus).
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

set_name <- snakemake@wildcards[["set"]]
spec <- snakemake@params[["spec"]]
settings <- snakemake@params[["settings"]]
plot_dir <- snakemake@output[["plots"]]
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

nom <- fread(snakemake@input[["nominal"]])
qt <- fread(snakemake@input[["permute"]])
n_qtl <- length(read_ids(snakemake@input[["donors"]]))
message(sprintf("%s / %s: %d FDR phenotypes, %d nominal rows, N = %d", spec$name, set_name,
                uniqueN(nom$phe_id), nrow(nom), n_qtl))

cols <- unlist(spec[c("snp_col", "bp_col", "beta_col", "se_col", "p_col", "n_col")])
cols <- cols[!is.na(cols) & nzchar(cols)]
gwas <- fread(snakemake@input[["gwas"]], select = unname(cols))
gwas <- standardize_gwas(gwas[get(spec$snp_col) %chin% unique(nom$var_id)], spec)
message(sprintf("GWAS SNPs overlapping QTL SNPs: %d", nrow(gwas)))

res <- run_coloc(nom, gwas, spec, n_qtl, settings, label = list(set = set_name))
write_tsv(res$summary, snakemake@output[["summary"]])
write_tsv(res$snps, snakemake@output[["snps"]])

annot <- if (length(snakemake@input[["annotation"]])) read_annotation(snakemake@input[["annotation"]]) else NULL
gold <- head(res$summary[tier == "Gold_Unique"][order(-PP.H4)], settings$max_locus_plots)
for (pid in gold$phenotype_id) {
  d <- nom[phe_id == pid & var_id %chin% gwas$snp]
  g <- gwas[snp %chin% d$var_id]
  p <- plot_locus(d, g, res$results[[pid]], title = pid,
                  subtitle = sprintf("%s, %s | PP.H4 = %.3f | lead SNP %s", spec$name, set_name,
                                     gold[phenotype_id == pid, PP.H4], gold[phenotype_id == pid, lead_snp]),
                  annot = annot, focus_gene = qt[phenotype_id == pid, gene][1])
  fname <- gsub("[^A-Za-z0-9_.-]+", "_", sub("->", "__to__", pid, fixed = TRUE))
  save_plot(p, file.path(plot_dir, paste0(fname, ".pdf")), width = 7, height = 9)
}
message(sprintf("Tested: %d; colocalised (PP.H4 >= %g): %d; locus plots: %d", nrow(res$summary), settings$pp_h4,
                sum(res$summary$PP.H4 >= settings$pp_h4), nrow(gold)))
