# Reading and annotating the pipeline's standard outputs
# (see docs/output_contract.md in the pipeline repository).

#' Annotate QTLtools permutation results with the phenotype annotation and
#' compute q-values.
#'
#' @param qt output of read_qtltools(..., "permute").
#' @param annotation phenotype annotation table (phenotype_id, set, gene, group, ...).
#' @param gt group table (read_group_table).
#' @param fdr_group_by "group", "class", "superclass" or "none" (one global q-value).
annotate_permute <- function(qt, annotation, gt, set_name, fdr_group_by = "group", fdr = 0.05) {
  qt <- qt[!is.na(adj_beta_pval)]
  ann <- copy(annotation)[, c("chrom", "start", "end", "strand") := NULL]
  ann[, set := NULL]
  out <- merge(qt[, .(phenotype_id = phe_id, grp_id, n_phe_in_grp, phe_chr, phe_from, phe_to,
                      n_var_in_cis, var_id, var_chr, var_from, dist_phe_var, slope, slope_se,
                      nom_pval, r_squared, adj_emp_pval, adj_beta_pval)],
               ann, by = "phenotype_id", all.x = TRUE, sort = FALSE)
  out[, set := set_name]
  out[is.na(group), group := "all"]
  out[, `:=`(class = map_group(group, gt, "class"), superclass = map_group(group, gt, "superclass"))]
  strata <- if (fdr_group_by == "none") rep("all", nrow(out)) else out[[fdr_group_by]]
  out[, qval := qvalue_by_group(adj_beta_pval, strata)]
  out[, fdr_pass := !is.na(qval) & qval <= fdr]
  setcolorder(out, c("set", "phenotype_id", "gene", "group", "class", "superclass"))
  setorder(out, qval)
  out[]
}

#' Per-set / per-group FDR summary with pi0 and counts at several thresholds.
fdr_summary <- function(qt, by = c("set", "group")) {
  qt[, .(n_tested = .N, pi0 = pi0_estimate(adj_beta_pval),
         n_fdr_01 = sum(qval <= 0.01, na.rm = TRUE), n_fdr = sum(fdr_pass),
         n_fdr_10 = sum(qval <= 0.10, na.rm = TRUE), n_fdr_genes = uniqueN(gene[fdr_pass & !is.na(gene)])),
     by = by]
}

#' Warnings about pi0, as in the original CAGE-seq q-value script.
pi0_diagnostic <- function(pi0) {
  if (!length(pi0) || is.na(pi0)) return("pi0 could not be estimated (too few phenotypes)")
  if (pi0 > 0.95) return("pi0 > 0.95: weak signal - over-correction (too many covariates), normalisation issue, or few true QTLs")
  if (pi0 < 0.05) return("pi0 < 0.05: nearly everything significant - possible under-correction or inflation")
  "ok"
}

#' Read the standard outputs of one pipeline run.
#' @param results_dir the pipeline's results directory.
#' @param sets phenotype sets to read (default: all found).
read_qtl_results <- function(results_dir, sets = NULL) {
  files <- list.files(file.path(results_dir, "qtl", "main"), pattern = "\\.permute\\.tsv\\.gz$", full.names = TRUE)
  found <- sub("\\.permute\\.tsv\\.gz$", "", basename(files))
  sets <- sets %||% found
  permute <- rbindlist(lapply(sets, function(s) fread(file.path(results_dir, "qtl", "main", paste0(s, ".permute.tsv.gz")))),
                       fill = TRUE)
  annotation <- rbindlist(lapply(sets, function(s) fread(file.path(results_dir, "phenotypes", paste0(s, ".annotation.tsv")))),
                          fill = TRUE)
  list(permute = permute, annotation = annotation, sets = sets)
}
