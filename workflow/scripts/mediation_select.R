# Colocalised QTLs strong enough for mediation analysis, and their lead SNPs.
# Only phenotypes with a partner (e.g. the other side of a CCI) are used.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

s <- snakemake@params[["settings"]]
gold <- fread(snakemake@input[[1]])
if (!"partner_gene" %in% names(gold)) gold[, `:=`(partner_gene = NA_character_, partner_group = NA_character_)]
hits <- gold[!is.na(lead_snp_pp_h4) & lead_snp_pp_h4 > s$snp_pp_h4 & PP.H4 > s$pp_h4 &
             !is.na(gwas_p) & gwas_p < s$gwas_p & !is.na(qtl_p) & qtl_p < s$qtl_p &
             !is.na(partner_gene) & !is.na(partner_group)]
hits <- hits[, .(gwas, set, phenotype_id, snp = lead_snp, group, gene, partner_group, partner_gene,
                 PP.H4, lead_snp_pp_h4, gwas_p, qtl_p)]
setorder(hits, gwas, set, -PP.H4)
message(sprintf("Colocalised hits passing mediation filters: %d (from %d Gold_Unique)", nrow(hits), nrow(gold)))
write_tsv(hits, snakemake@output[["hits"]])
writeLines(unique(hits$snp), snakemake@output[["snps"]])
