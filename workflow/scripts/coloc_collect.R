# Combine coloc results across GWAS and directions.
#   master: every colocalised CCI (PP.H4 >= pp_h4), with anchor annotations
#   counts: hits per GWAS x direction x tier
#   gold:   Gold_Unique hits (exactly one SNP with SNP.PP.H4 > 0.5)

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])

all <- rbindlist(lapply(unlist(snakemake@input[["summaries"]]), fread), fill = TRUE)
master <- all[PP.H4 >= snakemake@params[["pp_h4"]]]

# Anchor columns are added even when there are no hits, so downstream steps
# always see the same columns.
ids <- parse_cci_ids(as.character(master$cci))
master[, `:=`(sending_subclass = ids$sending_subclass, ligand = ids$ligand,
              receiving_subclass = ids$receiving_subclass, receptor = ids$receptor)]
master[, `:=`(
  subclass = fifelse(direction == "sending", sending_subclass, receiving_subclass),
  gene = fifelse(direction == "sending", ligand, receptor),
  partner_subclass = fifelse(direction == "sending", receiving_subclass, sending_subclass),
  partner_gene = fifelse(direction == "sending", receptor, ligand)
)]
setorder(master, gwas, direction, -PP.H4)

counts <- master[, .N, by = .(gwas, direction, tier)]
print(counts)
write_tsv(master, snakemake@output[["master"]])
write_tsv(counts, snakemake@output[["counts"]])
write_tsv(master[tier == "Gold_Unique"], snakemake@output[["gold"]])
