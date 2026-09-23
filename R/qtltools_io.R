# Column layouts of QTLtools `cis` output (v1.3). The pipeline always runs
# with --std-err; the layouts without slope_se are kept so that output from
# the original analysis (no --std-err) can still be read.

QTLTOOLS_COLS <- list(
  permute = c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
              "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
              "var_to", "dof1", "dof2", "bml1", "bml2", "nom_pval", "r_squared",
              "slope", "slope_se", "adj_emp_pval", "adj_beta_pval"),
  nominal = c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
              "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
              "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit")
)

#' Read QTLtools cis output (one or more files, plain or gzipped) with named columns.
read_qtltools <- function(paths, type = c("permute", "nominal")) {
  type <- match.arg(type)
  full <- QTLTOOLS_COLS[[type]]
  no_se <- setdiff(full, "slope_se")
  pieces <- lapply(paths, function(p) {
    if (!file.exists(p) || file.size(p) == 0) return(NULL)
    # a gzipped empty file has a non-zero size, so also catch fread's "File is empty"
    dt <- tryCatch(fread(p, header = FALSE), error = function(e) {
      if (grepl("empty", conditionMessage(e), ignore.case = TRUE)) return(NULL)
      stop(e)
    })
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt <- dt[V1 != "phe_id"]  # drop header lines written by chunk 0
    if (ncol(dt) == length(full)) {
      setnames(dt, full)
    } else if (ncol(dt) == length(no_se)) {
      setnames(dt, no_se)
      dt[, slope_se := NA_real_]
      setcolorder(dt, full)
    } else {
      stop(sprintf("%s: %d columns, expected %d or %d for QTLtools %s output",
                   p, ncol(dt), length(full), length(no_se), type))
    }
    dt
  })
  out <- rbindlist(pieces)
  if (nrow(out) == 0) {
    out <- as.data.table(setNames(replicate(length(full), character(0), simplify = FALSE), full))
  }
  num_cols <- setdiff(full, c("phe_id", "phe_chr", "phe_strd", "var_id", "var_chr"))
  out[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  out[]
}
