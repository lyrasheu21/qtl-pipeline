# Pick the colocalised cciQTLs strong enough for mediation analysis
# and list their lead SNPs.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

s <- snakemake@params[["settings"]]
gold <- fread(snakemake@input[[1]])

hits <- gold[!is.na(lead_snp_pp_h4) & lead_snp_pp_h4 > s$snp_pp_h4 &
             PP.H4 > s$pp_h4 &
             !is.na(gwas_p) & gwas_p < s$gwas_p &
             !is.na(qtl_p) & qtl_p < s$qtl_p]
hits <- hits[, .(gwas, direction, cci, snp = lead_snp, subclass, gene, partner_subclass, partner_gene,
                 PP.H4, lead_snp_pp_h4, gwas_p, qtl_p)]
setorder(hits, gwas, direction, -PP.H4)
message(sprintf("Colocalised hits passing mediation filters: %d (from %d Gold_Unique)", nrow(hits), nrow(gold)))

write_tsv(hits, snakemake@output[["hits"]])
writeLines(unique(hits$snp), snakemake@output[["snps"]])
