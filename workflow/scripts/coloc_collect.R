# Combine coloc results across GWAS and sets, joined to the phenotype
# annotation (gene, group, partner).
#   master: every colocalised phenotype (PP.H4 >= pp_h4)
#   counts: hits per GWAS x set x tier
#   gold:   Gold_Unique hits (exactly one SNP with SNP.PP.H4 > 0.5)
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

all <- rbindlist(lapply(unlist(snakemake@input[["summaries"]]), fread), fill = TRUE)
qt <- rbindlist(lapply(unlist(snakemake@input[["permute"]]), fread), fill = TRUE)
keep_cols <- intersect(c("set", "phenotype_id", "gene", "group", "class", "superclass",
                         "partner_gene", "partner_group"), names(qt))
ann <- unique(qt[, ..keep_cols], by = c("set", "phenotype_id"))
for (col in c("partner_gene", "partner_group")) if (!col %in% names(ann)) ann[, (col) := NA_character_]

master <- merge(all[PP.H4 >= snakemake@params[["pp_h4"]]], ann, by = c("set", "phenotype_id"), all.x = TRUE)
setcolorder(master, c("gwas", "set", "phenotype_id", "gene", "group"))
setorder(master, gwas, set, -PP.H4)
counts <- master[, .N, by = .(gwas, set, tier)]
print(counts)
write_tsv(master, snakemake@output[["master"]])
write_tsv(counts, snakemake@output[["counts"]])
write_tsv(master[tier == "Gold_Unique"], snakemake@output[["gold"]])
