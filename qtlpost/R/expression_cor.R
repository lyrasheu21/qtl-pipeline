# Correlation of phenotype values with the expression of their gene.

#' Pearson correlation across donors between each phenotype and the expression
#' of its `gene` in its `group` (for CCIs: the CCI score vs ligand expression
#' in the sending subclass, or vs receptor expression in the receiving one).
#'
#' @param pheno numeric matrix, phenotypes x donors.
#' @param annotation data.table with phenotype_id, group, gene.
#' @param get_expression function(group) returning a genes x donors matrix (or NULL).
phenotype_expression_cor <- function(pheno, annotation, get_expression, min_donors = 4) {
  annotation <- annotation[phenotype_id %in% rownames(pheno) & !is.na(gene)]
  out <- lapply(unique(annotation$group), function(g) {
    expr <- get_expression(g)
    if (is.null(expr)) return(NULL)
    jobs <- annotation[group == g & gene %in% rownames(expr)]
    donors <- intersect(colnames(pheno), colnames(expr))
    if (!nrow(jobs) || length(donors) < min_donors) return(NULL)
    x <- pheno[jobs$phenotype_id, donors, drop = FALSE]
    y <- expr[jobs$gene, donors, drop = FALSE]
    ok <- !is.na(x) & !is.na(y)
    n_ok <- rowSums(ok)
    r <- vapply(seq_len(nrow(jobs)), function(i) {
      if (n_ok[i] < min_donors) return(NA_real_)
      suppressWarnings(cor(x[i, ok[i, ]], y[i, ok[i, ]]))
    }, numeric(1))
    jobs[, .(phenotype_id, group, gene, n_donors = n_ok, r = r)]
  })
  res <- rbindlist(out)
  if (!nrow(res)) res <- data.table(phenotype_id = character(), group = character(), gene = character(),
                                    n_donors = integer(), r = numeric())
  res
}
