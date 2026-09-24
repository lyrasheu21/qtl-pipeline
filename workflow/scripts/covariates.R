# QTLtools covariate file for one phenotype set and covariate set. Blocks:
#   metadata covariates (with transforms), genotype PCs, hidden factors
#   (a factors file such as exported PEER factors, or phenotype PCs).
# Donors missing any covariate are dropped and reported: QTLtools itself only
# uses samples present in all of VCF, BED and covariate file.
source(file.path(snakemake@scriptdir, "load_qtlpost.R"))

cs <- snakemake@params[["covset"]]
hidden <- snakemake@params[["hidden"]]
bed <- read_bed_matrix(snakemake@input[["bed"]])
donors <- colnames(bed$mat)
blocks <- list()

use <- unlist(cs$use)
if (length(use)) {
  meta <- read_table_any(snakemake@input[["metadata"]])
  blocks$metadata <- metadata_covariates(meta, snakemake@params[["metadata_id"]], donors,
                                         snakemake@params[["definitions"]], use)
}

n_gpc <- as.integer(cs$genotype_pcs)
if (n_gpc > 0) {
  if (is.null(snakemake@input[["pcs"]])) stop("covariates.genotype_pcs > 0 but genotypes.pcs.source is 'none'")
  pcs <- read_tsv(snakemake@input[["pcs"]])
  pc_mat <- as.matrix(pcs[, -1])
  rownames(pc_mat) <- as.character(pcs[[1]])
  blocks$genotype_pcs <- factor_covariates(pc_mat, donors, n_gpc, "genotype_PC")
}

k <- as.integer(cs$hidden_k)
if (k > 0) {
  source_type <- hidden$source %||% "none"
  if (source_type == "file") {
    smap <- hidden$sample_map
    lookup <- donors
    map <- NULL
    if (!is.null(smap) && length(smap) && !is.null(snakemake@input[["factor_map"]])) {
      map <- read_table_any(snakemake@input[["factor_map"]])
      lookup <- as.character(map[[smap$from]])
    }
    f <- read_sample_factors(snakemake@input[["factors"]], lookup,
                             object = if (nzchar(hidden$object %||% "")) hidden$object else NULL)
    if (!is.null(map)) f <- t(apply_sample_map(t(f), map, smap$from, smap$to))  # averages replicates
    blocks$hidden <- factor_covariates(f, donors, k, "hidden_")
  } else if (source_type == "phenotype_pca") {
    blocks$hidden <- factor_covariates(phenotype_pcs(bed$mat, k), donors, k, "phenotype_PC")
  } else {
    stop("hidden_k = ", k, " but covariates.hidden_factors.source is '", source_type, "'")
  }
}

res <- assemble_covariates(blocks, donors)
if (length(res$dropped)) {
  message(sprintf("Dropped %d donors with missing covariates: %s", length(res$dropped),
                  paste(head(res$dropped, 20), collapse = ", ")))
}
if (length(res$donors) < snakemake@params[["min_donors"]]) {
  stop(sprintf("Only %d donors have phenotypes, genotypes and all covariates (min_donors = %d)",
               length(res$donors), snakemake@params[["min_donors"]]))
}
if (nrow(res$cov) == 0) {
  # the QTLtools rules always pass --cov, which needs at least one covariate
  stop("No covariates selected for covariate set '", snakemake@wildcards[["covset"]],
       "'. Add metadata covariates, genotype PCs or hidden factors.")
}
write_tsv(res$cov, snakemake@output[["cov"]])
writeLines(res$donors, snakemake@output[["donors"]])
report <- data.table(set = snakemake@wildcards[["set"]], covset = snakemake@wildcards[["covset"]],
                     n_donors_in_bed = length(donors), n_donors_with_covariates = length(res$donors),
                     n_dropped = length(res$dropped), n_covariates = nrow(res$cov),
                     covariates = paste(res$cov$id, collapse = ","))
print(report)
write_tsv(report, snakemake@output[["report"]])
