# Covariate assembly for QTLtools. A covariate file is built from optional
# blocks: metadata covariates (with transforms), genotype PCs, and hidden
# factors (e.g. PEER factors or phenotype PCs).

apply_transform <- function(x, transform, name) {
  switch(transform,
    none = x,
    scale = as.numeric(scale(as.numeric(x))),
    log = {
      x <- as.numeric(x)
      if (any(x <= 0, na.rm = TRUE)) stop("Covariate '", name, "': log transform of non-positive values")
      log(x)
    },
    factor = as.character(x),
    stop("Unknown transform '", transform, "' for covariate '", name, "'")
  )
}

#' Metadata covariates as a list of per-donor vectors (in `donors` order).
#' @param definitions named list: output name -> list(column, transform).
#' @param use names of `definitions` to include.
metadata_covariates <- function(meta, meta_id, donors, definitions, use) {
  unknown <- setdiff(use, names(definitions))
  if (length(unknown)) stop("Covariates not defined: ", paste(unknown, collapse = ", "))
  meta <- meta[match(donors, as.character(meta[[meta_id]]))]
  rows <- list()
  for (name in use) {
    def <- definitions[[name]]
    if (!def$column %in% names(meta)) stop("Metadata column not found: ", def$column)
    rows[[name]] <- apply_transform(meta[[def$column]], def$transform, name)
  }
  rows
}

#' First `k` columns of a samples x factors matrix as per-donor vectors.
factor_covariates <- function(mat, donors, k, prefix) {
  if (k <= 0) return(list())
  if (ncol(mat) < k) stop(sprintf("Requested %d %s factors but only %d are available", k, prefix, ncol(mat)))
  rows <- list()
  for (i in seq_len(k)) rows[[paste0(prefix, i)]] <- as.numeric(mat[match(donors, rownames(mat)), i])
  rows
}

#' Top `k` principal components of a phenotype matrix (rows = phenotypes,
#' columns = donors), after rank-normalising each phenotype. Returns a
#' donors x k matrix.
phenotype_pcs <- function(mat, k) {
  if (k <= 0) return(matrix(numeric(0), nrow = ncol(mat), ncol = 0, dimnames = list(colnames(mat), NULL)))
  x <- rank_normal_rows(mat)
  x[is.na(x)] <- 0
  x <- x[apply(x, 1, sd) > 0, , drop = FALSE]
  k_use <- min(k, ncol(x) - 1, nrow(x))
  if (k_use < k) warning(sprintf("Only %d phenotype PCs can be computed (requested %d)", k_use, k))
  pc <- prcomp(t(x), center = TRUE, scale. = FALSE, rank. = k_use)$x
  rownames(pc) <- colnames(mat)
  colnames(pc) <- paste0("PC", seq_len(ncol(pc)))
  pc
}

#' Combine covariate blocks into a QTLtools covariate table.
#'
#' Donors with a missing value in any covariate are dropped (QTLtools only
#' uses samples present in all input files, so this is reported, not fatal).
#' @param blocks list of named per-donor vectors (from metadata_covariates /
#'   factor_covariates), all in `donors` order.
#' @return list(cov = data.table with `id` + one column per kept donor,
#'   donors = kept donors, dropped = dropped donors)
assemble_covariates <- function(blocks, donors) {
  rows <- do.call(c, unname(blocks))
  if (!length(rows)) {
    return(list(cov = data.table(id = character()), donors = donors, dropped = character()))
  }
  complete <- Reduce(`&`, lapply(rows, function(v) !is.na(v)))
  kept <- donors[complete]
  out <- as.data.table(do.call(rbind, lapply(rows, function(v) as.character(v[complete]))))
  setnames(out, kept)
  out[, id := names(rows)]
  setcolorder(out, "id")
  list(cov = out, donors = kept, dropped = donors[!complete])
}

#' Build a covariate table from metadata and ancestry/genotype PCs (kept for
#' the original CCI interface; see assemble_covariates for the general case).
#' Unlike assemble_covariates, any missing value is an error.
build_covariates <- function(meta, meta_id, pcs, pc_id, donors, definitions, use,
                             n_pcs, pc_prefix = "PC") {
  blocks <- list(metadata = metadata_covariates(meta, meta_id, donors, definitions, use))
  if (n_pcs > 0) {
    pc_cols <- paste0(pc_prefix, seq_len(n_pcs))
    missing <- setdiff(pc_cols, names(pcs))
    if (length(missing)) stop("Ancestry PC columns not found: ", paste(missing, collapse = ", "))
    pc_mat <- as.matrix(pcs[, ..pc_cols])
    rownames(pc_mat) <- as.character(pcs[[pc_id]])
    blocks$pcs <- factor_covariates(pc_mat, donors, n_pcs, "ancestry_PC")
  }
  rows <- do.call(c, unname(blocks))
  has_na <- vapply(rows, function(v) any(is.na(v)), logical(1))
  if (any(has_na)) {
    bad <- names(rows)[has_na]
    donors_na <- unique(unlist(lapply(rows[bad], function(v) donors[is.na(v)])))
    stop("Missing covariate values for ", paste(bad, collapse = ", "), " in donors: ",
         paste(head(donors_na, 20), collapse = ", "))
  }
  assemble_covariates(blocks, donors)$cov
}
