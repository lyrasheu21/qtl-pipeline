# Generate the small synthetic datasets used by tests/configs/*.yaml.
#
#   Rscript tests/data/generate_toy_data.R
#
# Shared: 60 donors, 2 chromosomes, 6 genes, GWAS, eQTL and expression files.
#   CCI    (tests/data/cci_matrix.tsv.gz): ~20 CCIs for the cci adapter
#   CAGE   (tests/data/cage/): peak BED (5 peaks per gene, gid = gene) for the
#          bed adapter, PEER-like factors in an .RData keyed by assay sample
#          IDs (one replicate pair) + sample map, genotype PCs as a .pca file
#   matrix (tests/data/matrix/): expression matrix with assay sample IDs,
#          coordinates table and sample map for the matrix adapter
# Two causal SNPs are planted:
#   rs_causal_1 (chr1, near LIGA) raises every LIGA-containing CCI score, LIGA
#     expression, the TOY_CC GWAS signal and the LIGA eQTL; RECA expression
#     partly follows LIGA (a mediation signal).
#   rs_causal_2 (chr2, near RECB) does the same for RECB and the TOY_QUANT GWAS.
# Only base R and data.table are needed.

suppressPackageStartupMessages(library(data.table))
set.seed(20240313)

out_dir <- "tests/data"
dir.create(file.path(out_dir, "eqtl"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "expression"), showWarnings = FALSE, recursive = TRUE)

n_donors <- 60
donors <- sprintf("D%03d", seq_len(n_donors))
subclasses <- c("Astro", "OPC", "Micro", "EN_L2_3_IT", "IN_VIP")

# Genes -----------------------------------------------------------------------
genes <- data.table(
  gene_name = c("LIGA", "RECA", "LIGB", "RECB", "LIGC", "RECC"),
  gene_id = sprintf("ENSG%011d", 1:6),
  chrom = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"),
  start = c(1000000, 1600000, 3000000, 2000000, 4000000, 5000000),
  strand = "+"
)
genes[, end := start + 20000]
setcolorder(genes, c("gene_name", "gene_id", "chrom", "start", "end", "strand"))
fwrite(genes, file.path(out_dir, "gene_annotation.tsv"), sep = "\t")

# Genotypes -------------------------------------------------------------------
snps <- rbindlist(lapply(seq_len(nrow(genes)), function(i) {
  data.table(chrom = sub("chr", "", genes$chrom[i]),
             pos = sort(sample(seq(genes$start[i] - 200000, genes$start[i] + 200000, by = 97), 40)))
}))
snps <- unique(snps, by = c("chrom", "pos"))[order(as.integer(chrom), pos)]
snps[, id := sprintf("rs%d", 100000 + .I)]
snps[, maf := runif(.N, 0.1, 0.5)]
causal <- c(snps[chrom == "1"][which.min(abs(pos - 1005000)), id],
            snps[chrom == "2"][which.min(abs(pos - 2010000)), id])
snps[id == causal[1], `:=`(id = "rs_causal_1", maf = 0.4)]
snps[id == causal[2], `:=`(id = "rs_causal_2", maf = 0.4)]
causal <- c("rs_causal_1", "rs_causal_2")

geno <- t(vapply(snps$maf, function(m) rbinom(n_donors, 2, m), numeric(n_donors)))
dimnames(geno) <- list(snps$id, donors)
gt <- matrix(c("0/0", "0/1", "1/1")[geno + 1], nrow = nrow(geno), dimnames = dimnames(geno))
# QC targets: two high-missingness SNPs (fail --geno) and one monomorphic SNP (fails --maf)
gt[3, sample(n_donors, 15)] <- "./."
gt[50, sample(n_donors, 15)] <- "./."
gt[20, ] <- "0/0"

vcf_path <- file.path(out_dir, "genotypes.vcf.gz")
con <- gzfile(vcf_path, "w")
writeLines(c("##fileformat=VCFv4.2",
             "##contig=<ID=1,length=10000000>", "##contig=<ID=2,length=10000000>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", donors), collapse = "\t")),
           con)
body <- cbind(snps$chrom, snps$pos, snps$id, "A", "G", ".", "PASS", ".", "GT", gt)
writeLines(apply(body, 1, paste, collapse = "\t"), con)
close(con)

# Donor 60 is excluded by the keep list (tests plink --keep).
fwrite(data.table(donors[-60], donors[-60]), file.path(out_dir, "keep.txt"), sep = "\t", col.names = FALSE)

# CCI matrix -----------------------------------------------------------------
lr <- data.table(ligand = c("LIGA", "LIGB", "LIGC", "LIGA"), receptor = c("RECA", "RECB", "RECC", "RECB"))
ccis <- rbindlist(lapply(seq_len(nrow(lr)), function(i) {
  pairs <- CJ(send = subclasses, recv = subclasses)[sample(.N, 5)]
  pairs[, cci := sprintf("%s:%s->%s:%s", send, lr$ligand[i], recv, lr$receptor[i])]
}))$cci
ccis <- c(ccis, "Astro:LIGA->OPC:ITGA1_ITGB1")  # complex receptor: absent from annotation, dropped

g1 <- geno["rs_causal_1", ]
g2 <- geno["rs_causal_2", ]
extra_donors <- c("D901", "D902")  # in the CCI matrix but not genotyped
all_cci_donors <- c(donors, extra_donors)
score <- t(vapply(ccis, function(id) {
  s <- rnorm(length(all_cci_donors), 5, 1)
  g1x <- c(g1, rbinom(2, 2, 0.4))
  g2x <- c(g2, rbinom(2, 2, 0.4))
  if (grepl(":LIGA->", id, fixed = TRUE)) s <- s + 1.2 * g1x
  if (grepl(":RECB$", id)) s <- s + 1.2 * g2x
  s[sample(length(s), 3)] <- 0  # a few zero scores
  s
}, numeric(length(all_cci_donors))))
colnames(score) <- all_cci_donors
cci_dt <- data.table(cci = ccis, round(score, 4))
fwrite(cci_dt, file.path(out_dir, "cci_matrix.tsv.gz"), sep = "\t")

# Metadata and ancestry PCs --------------------------------------------------
meta <- data.table(SubID = donors,
                   PMI = round(runif(n_donors, 5, 30), 1),
                   Age = round(runif(n_donors, 40, 95)),
                   Sex = sample(c("M", "F"), n_donors, replace = TRUE),
                   prep = sample(c("A", "B"), n_donors, replace = TRUE),
                   Brain_bank = sample(c("BB1", "BB2", "BB3"), n_donors, replace = TRUE),
                   avg_n_counts = round(runif(n_donors, 2000, 9000)),
                   PC1 = rnorm(n_donors), PC2 = rnorm(n_donors), PC3 = rnorm(n_donors))
fwrite(meta, file.path(out_dir, "metadata.tsv"), sep = "\t")

pcs <- data.table(SubID = donors, matrix(rnorm(n_donors * 10), ncol = 10, dimnames = list(NULL, paste0("PC", 1:10))))
fwrite(pcs, file.path(out_dir, "ancestry_pcs.tsv"), sep = "\t")

# GWAS -------------------------------------------------------------------------
make_gwas <- function(signal_snp, z_signal) {
  z <- rnorm(nrow(snps))
  z[snps$id == signal_snp] <- z_signal
  se <- 0.02
  data.table(SNP = snps$id, CHR = snps$chrom, BP = snps$pos, beta = z * se, se = se,
             p = 2 * pnorm(-abs(z)), N = 20000)
}
fwrite(make_gwas("rs_causal_1", 9), file.path(out_dir, "gwas_toy_cc.tsv.gz"), sep = "\t")
fwrite(make_gwas("rs_causal_2", -8), file.path(out_dir, "gwas_toy_quant.tsv.gz"), sep = "\t")

# Expression residuals ---------------------------------------------------------
for (sc in subclasses) {
  e <- matrix(rnorm(nrow(genes) * n_donors), nrow = nrow(genes), dimnames = list(genes$gene_name, donors))
  e["LIGA", ] <- 0.9 * scale(g1)[, 1] + rnorm(n_donors, 0, 0.6)
  e["RECA", ] <- 0.8 * e["LIGA", ] + rnorm(n_donors, 0, 0.6)
  e["RECB", ] <- 0.9 * scale(g2)[, 1] + rnorm(n_donors, 0, 0.6)
  e["LIGB", ] <- 0.8 * e["RECB", ] + rnorm(n_donors, 0, 0.6)
  fwrite(data.table(gene = rownames(e), round(e, 4)), file.path(out_dir, "expression", paste0(sc, ".tsv.gz")), sep = "\t")
}

# eQTL resource (headerless: chr start end snp ref alt gene beta se z) ---------
for (sc in subclasses) for (chr in c("1", "2")) {
  s <- snps[chrom == chr]
  gs <- genes[sub("chr", "", chrom) == chr, gene_name]
  rows <- CJ(snp = s$id, gene = gs)
  rows[, pos := s$pos[match(snp, s$id)]]
  rows[, z := rnorm(.N)]
  rows[snp == "rs_causal_1" & gene == "LIGA", z := 7]
  rows[snp == "rs_causal_2" & gene == "RECB", z := -7]
  rows[, se := 0.05]
  out <- rows[, .(chr, pos, pos + 1, snp, "A", "G", gene, round(z * se, 5), se, round(z, 4))]
  fwrite(out, file.path(out_dir, "eqtl", sprintf("eqtl_%s_chr%s.tsv.gz", sc, chr)), sep = "\t", col.names = FALSE)
}

# Cell types used by the toy data ----------------------------------------------
ct <- fread("config/cell_types.tsv")[subclass %in% subclasses]
fwrite(ct, file.path(out_dir, "cell_types.tsv"), sep = "\t")



# Extra groups used by the CAGE ("DLPFC") and matrix ("BulkBrain") examples -----
set.seed(163)
for (grp in c("DLPFC", "BulkBrain")) {
  e <- matrix(rnorm(nrow(genes) * n_donors), nrow = nrow(genes), dimnames = list(genes$gene_name, donors))
  e["LIGA", ] <- 0.9 * scale(g1)[, 1] + rnorm(n_donors, 0, 0.6)
  e["RECB", ] <- 0.9 * scale(g2)[, 1] + rnorm(n_donors, 0, 0.6)
  fwrite(data.table(gene = rownames(e), round(e, 4)), file.path(out_dir, "expression", paste0(grp, ".tsv.gz")), sep = "\t")
  for (chr in c("1", "2")) {
    s_chr <- snps[chrom == chr]
    gs <- genes[sub("chr", "", chrom) == chr, gene_name]
    rows <- CJ(snp = s_chr$id, gene = gs)
    rows[, pos := s_chr$pos[match(snp, s_chr$id)]]
    rows[, z := rnorm(.N)]
    rows[snp == "rs_causal_1" & gene == "LIGA", z := 7]
    rows[snp == "rs_causal_2" & gene == "RECB", z := -7]
    rows[, se := 0.05]
    out <- rows[, .(chr, pos, pos + 1, snp, "A", "G", gene, round(z * se, 5), se, round(z, 4))]
    fwrite(out, file.path(out_dir, "eqtl", sprintf("eqtl_%s_chr%s.tsv.gz", grp, chr)), sep = "\t", col.names = FALSE)
  }
}

# CAGE-seq style peaks (bed adapter) --------------------------------------------
set.seed(1631)
cage_dir <- file.path(out_dir, "cage")
dir.create(cage_dir, showWarnings = FALSE)
hidden <- matrix(rnorm(n_donors * 3), ncol = 3, dimnames = list(donors, NULL))  # unmeasured confounders
peaks <- rbindlist(lapply(seq_len(nrow(genes)), function(i) {
  data.table(chr = sub("chr", "", genes$chrom[i]), start = genes$start[i] + c(-300, 0, 150, 400, 900),
             gene = genes$gene_name[i], k = 1:5)
}))
peaks[, `:=`(end = start + 100, pid = sprintf("peak_%s_%d", gene, k), strand = "+")]
loadings <- matrix(rnorm(nrow(peaks) * 3, sd = 1.5), ncol = 3)
values <- 10 + loadings %*% t(hidden) + matrix(rnorm(nrow(peaks) * n_donors), nrow = nrow(peaks))
values[peaks$gene == "LIGA" & peaks$k <= 2, ] <- values[peaks$gene == "LIGA" & peaks$k <= 2, ] + rep(1.5 * g1, each = 2)
values[peaks$gene == "RECB" & peaks$k <= 3, ] <- values[peaks$gene == "RECB" & peaks$k <= 3, ] + rep(1.5 * g2, each = 3)
colnames(values) <- donors
cage_bed <- cbind(peaks[, .(`#chr` = chr, start, end, pid, gid = gene, strand)], as.data.table(round(values, 4)))
setorder(cage_bed, `#chr`, start)
fwrite(cage_bed, file.path(cage_dir, "phenotypes.bed.gz"), sep = "\t")

# PEER-like hidden factors keyed by CAGE sample IDs, one donor sequenced twice
cage_ids <- sprintf("CAGE_S%03d", seq_len(n_donors))
sample_map <- data.table(sample_id = c(cage_ids, "CAGE_S005_rep"), donor_id = c(donors, "D005"))
fwrite(sample_map, file.path(cage_dir, "sample_map.csv"))
f_true <- rbind(hidden, hidden["D005", , drop = FALSE])
hidden_factors <- cbind(f_true + matrix(rnorm(length(f_true), sd = 0.1), ncol = 3),
                        matrix(rnorm(nrow(f_true) * 2), ncol = 2))
rownames(hidden_factors) <- sample_map$sample_id
save(hidden_factors, file = file.path(cage_dir, "peer_factors.RData"))

# Genotype PCs with PCs as rows and donors as columns (like a .pca file)
pca <- matrix(rnorm(5 * n_donors), nrow = 5, dimnames = list(sprintf("subset_toy_PC%d", 1:5), donors))
fwrite(data.table(SampleID = rownames(pca), round(pca, 5)), file.path(cage_dir, "genotype.pca"), sep = "\t")

# eQTL-style expression matrix (matrix adapter) ---------------------------------
set.seed(1632)
mat_dir <- file.path(out_dir, "matrix")
dir.create(mat_dir, showWarnings = FALSE)
extra <- data.table(feature_id = sprintf("GENE%02d", 7:16), chrom = rep(c("chr1", "chr2"), each = 5),
                    start = c(seq(1200000, 3800000, length.out = 5), seq(2300000, 5600000, length.out = 5)))
coords <- rbind(genes[, .(feature_id = gene_name, chrom, start)], extra)
coords[, `:=`(end = start + 15000, strand = "+", gene = feature_id)]
fwrite(coords, file.path(mat_dir, "coordinates.tsv"), sep = "\t")
hidden2 <- matrix(rnorm(n_donors * 2), ncol = 2)
expr <- 5 + matrix(rnorm(nrow(coords) * 2, sd = 1.2), ncol = 2) %*% t(hidden2) +
  matrix(rnorm(nrow(coords) * n_donors), nrow = nrow(coords))
expr[coords$feature_id == "LIGA", ] <- expr[coords$feature_id == "LIGA", ] + 1.5 * g1
expr[coords$feature_id == "RECB", ] <- expr[coords$feature_id == "RECB", ] + 1.5 * g2
rna_ids <- sprintf("RNA_%s", donors)
colnames(expr) <- rna_ids
expr <- cbind(expr, RNA_D010_rep = expr[, "RNA_D010"] + rnorm(nrow(coords), sd = 0.2))
fwrite(data.table(feature_id = coords$feature_id, round(expr, 4)), file.path(mat_dir, "expression.tsv.gz"), sep = "\t")
fwrite(data.table(rna_sample = c(rna_ids, "RNA_D010_rep"), donor = c(donors, "D010")), file.path(mat_dir, "sample_map.tsv"), sep = "\t")

message("Toy data written to ", out_dir, ": ", nrow(snps), " SNPs, ", length(ccis), " CCIs, ",
        nrow(peaks), " CAGE peaks, ", nrow(coords), " expression features, ", n_donors, " donors")
