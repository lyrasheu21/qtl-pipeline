# Building the QTLtools covariate file.

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

#' Build a covariate matrix: rows = covariates, columns = donors (in `donors` order).
#'
#' @param meta,pcs data.tables with donor ID columns `meta_id` / `pc_id`.
#' @param definitions named list: output name -> list(column, transform).
#' @param use names of `definitions` to include.
#' @param n_pcs number of ancestry PCs (columns `<pc_prefix>1..n`), 0 for none.
#' @return data.table with an `id` column followed by one column per donor.
build_covariates <- function(meta, meta_id, pcs, pc_id, donors, definitions, use,
                             n_pcs, pc_prefix = "PC") {
  unknown <- setdiff(use, names(definitions))
  if (length(unknown)) stop("Covariates not defined: ", paste(unknown, collapse = ", "))

  meta <- meta[match(donors, as.character(meta[[meta_id]]))]
  rows <- list()
  for (name in use) {
    def <- definitions[[name]]
    if (!def$column %in% names(meta)) stop("Metadata column not found: ", def$column)
    rows[[name]] <- apply_transform(meta[[def$column]], def$transform, name)
  }

  if (n_pcs > 0) {
    pc_cols <- paste0(pc_prefix, seq_len(n_pcs))
    missing <- setdiff(pc_cols, names(pcs))
    if (length(missing)) stop("Ancestry PC columns not found: ", paste(missing, collapse = ", "))
    pcs <- pcs[match(donors, as.character(pcs[[pc_id]]))]
    for (i in seq_len(n_pcs)) rows[[paste0("ancestry_PC", i)]] <- as.numeric(pcs[[pc_cols[i]]])
  }

  has_na <- vapply(rows, function(v) any(is.na(v)), logical(1))
  if (any(has_na)) {
    bad <- names(rows)[has_na]
    donors_na <- unique(unlist(lapply(rows[bad], function(v) donors[is.na(v)])))
    stop("Missing covariate values for ", paste(bad, collapse = ", "), " in donors: ",
         paste(head(donors_na, 20), collapse = ", "))
  }

  out <- as.data.table(do.call(rbind, lapply(rows, as.character)))
  setnames(out, donors)
  out[, id := names(rows)]
  setcolorder(out, "id")
  out[]
}
