# Phenotype adapters: build a QTLtools BED and the phenotype annotation table
# for one phenotype set.
#
#   bed     an existing QTLtools BED (e.g. CAGE-seq peaks); gene from an
#           optional annotation table or the BED gid column
#   matrix  feature x sample table + coordinates table (eQTL/caQTL/pQTL)
#   cci     donor x CCI matrix; each CCI placed at its ligand (sending) or
#           receptor (receiving) gene
#
# Shared steps: optional sample map (assay ID -> donor, replicates averaged),
# restriction to genotyped donors, min_nonzero_fraction filter, chromosome
# filter, sorting. The report counts what was dropped at each step.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

set_name <- snakemake@wildcards[["set"]]
adapter <- snakemake@params[["adapter"]]
spec <- snakemake@params[["spec"]]
min_fraction <- as.numeric(snakemake@params[["min_fraction"]])
chromosomes <- as.character(unlist(snakemake@params[["chromosomes"]]))
vcf_samples <- read_ids(snakemake@input[["vcf_samples"]])
opt <- function(name, default = NULL) if (is.null(spec[[name]]) || identical(spec[[name]], "")) default else spec[[name]]

# ---- adapter-specific: phenotype matrix + coordinates + annotation -----------
if (adapter == "bed") {
  bed <- read_bed_matrix(snakemake@input[["bed"]])
  mat <- bed$mat
  coords <- data.table(phenotype_id = bed$info$pid, chrom = bed$info$chr, start = bed$info$start,
                       end = bed$info$end, gid = bed$info$gid, strand = bed$info$strand)
  if (!is.null(snakemake@input[["annotation"]])) {
    ann <- fread(snakemake@input[["annotation"]])
    ann <- ann[match(coords$phenotype_id, ann$phenotype_id)]
    coords[, gene := ann$gene]
    coords[, group := if ("group" %in% names(ann)) ann$group else opt("group", "all")]
    for (col in c("partner_gene", "partner_group")) if (col %in% names(ann)) coords[, (col) := ann[[col]]]
  } else {
    coords[, gene := fifelse(gid != phenotype_id, gid, NA_character_)]
    coords[, group := opt("group", "all")]
  }
} else if (adapter == "matrix") {
  mat <- read_feature_matrix(snakemake@input[["matrix"]], opt("id_column", "feature_id"))
  cd <- fread(snakemake@input[["coordinates"]])
  id_col <- opt("coordinates_id_column", "feature_id")
  if (!id_col %in% names(cd)) stop("Coordinates table has no column '", id_col, "'")
  cd <- cd[match(rownames(mat), cd[[id_col]])]
  coords <- data.table(phenotype_id = rownames(mat), chrom = cd$chrom, start = cd$start, end = cd$end,
                       gid = rownames(mat), strand = if ("strand" %in% names(cd)) cd$strand else "+",
                       gene = if ("gene" %in% names(cd)) cd$gene else rownames(mat),
                       group = if ("group" %in% names(cd)) cd$group else opt("group", "all"))
} else if (adapter == "cci") {
  mat <- read_cci_matrix(snakemake@input[["matrix"]], opt("id_column", "cci"))
  direction <- snakemake@params[["direction"]]
  ann <- cci_annotation(rownames(mat), direction)
  genes <- read_annotation(snakemake@input[["gene_annotation"]])
  g <- genes[match(ann$gene, gene_name)]
  coords <- data.table(phenotype_id = ann$phenotype_id, chrom = g$chrom, start = g$start, end = g$end,
                       gid = ann$phenotype_id, strand = "+", gene = ann$gene, group = ann$group,
                       partner_gene = ann$partner_gene, partner_group = ann$partner_group,
                       ligand = ann$ligand, receptor = ann$receptor,
                       sending_group = ann$sending_group, receiving_group = ann$receiving_group)
} else {
  stop("Unknown adapter: ", adapter)
}
n_input <- nrow(mat)
n_samples_input <- ncol(mat)

# ---- shared steps ------------------------------------------------------------
smap <- spec[["sample_map"]]
if (!is.null(smap) && length(smap) && !is.null(snakemake@input[["sample_map"]])) {
  mat <- apply_sample_map(mat, read_table_any(snakemake@input[["sample_map"]]), smap$from, smap$to)
}
donors <- intersect(colnames(mat), vcf_samples)
if (!length(donors)) {
  stop("No phenotype samples match the genotype VCF (examples - phenotype: ",
       paste(head(colnames(mat), 3), collapse = ", "), "; VCF: ", paste(head(vcf_samples, 3), collapse = ", "),
       "). Use a sample_map if the IDs differ.")
}
mat <- mat[, donors, drop = FALSE]

nonzero <- rowMeans(!is.na(mat) & mat != 0)
if (isTRUE(opt("zero_as_missing", FALSE))) mat[mat == 0] <- NA
coords[, chrom_norm := sub("^chr", "", as.character(chrom), ignore.case = TRUE)]
has_coords <- !is.na(coords$chrom) & !is.na(coords$start)
on_chrom <- has_coords & coords$chrom_norm %in% chromosomes
pass_fraction <- nonzero >= min_fraction
keep <- pass_fraction & on_chrom

out <- data.table(chr = coords$chrom_norm, start = as.integer(coords$start), end = as.integer(coords$end),
                  pid = coords$phenotype_id, gid = coords$gid, strand = coords$strand)[keep]
out <- cbind(out, as.data.table(mat[keep, , drop = FALSE]))
out <- out[order(match(chr, chromosomes), start)]
setnames(out, "chr", "#chr")
fwrite(out, snakemake@output[["bed"]], sep = "\t", na = "NA", quote = FALSE)

annotation <- coords[, !c("gid", "chrom_norm")]
annotation[, set := set_name]
setcolorder(annotation, c("phenotype_id", "set", "chrom", "start", "end", "strand", "gene", "group"))
write_tsv(annotation, snakemake@output[["annotation"]])

report <- data.table(
  set = set_name, adapter = adapter, n_phenotypes_input = n_input, n_samples_input = n_samples_input,
  n_donors_genotyped = length(donors), below_min_nonzero_fraction = sum(!pass_fraction),
  no_coordinates = sum(pass_fraction & !has_coords),
  not_on_analysed_chromosome = sum(pass_fraction & has_coords & !on_chrom),
  n_phenotypes_in_bed = sum(keep), min_nonzero_fraction = min_fraction
)
print(report)
write_tsv(report, snakemake@output[["report"]])
