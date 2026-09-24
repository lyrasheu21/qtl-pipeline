# Column layouts of QTLtools `cis` output (v1.3.1, checked against the source
# in src/mode_cis/cis_permutation_pass.cpp and cis_nominal_pass.cpp).
#
#   --std-err adds slope_se right after slope.
#   --grp-best writes the group ID in column 1 and inserts two columns after
#   the strand: the best phenotype in the group and the group size.
#
# read_qtltools() always returns the same columns: phe_id (the phenotype; for
# grouped runs the group's best phenotype), grp_id and n_phe_in_grp (NA when
# ungrouped), then the standard columns.

QTLTOOLS_HEAD <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd")
QTLTOOLS_TAIL <- list(
  permute = c("n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to",
              "dof1", "dof2", "bml1", "bml2", "nom_pval", "r_squared", "slope", "slope_se",
              "adj_emp_pval", "adj_beta_pval"),
  nominal = c("n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from", "var_to",
              "nom_pval", "r_squared", "slope", "slope_se", "best_hit")
)
QTLTOOLS_OUT <- lapply(QTLTOOLS_TAIL, function(tail) c(QTLTOOLS_HEAD, "grp_id", "n_phe_in_grp", tail))

#' Possible raw layouts for one output type: with/without --std-err, with/without --grp-best.
qtltools_layouts <- function(type) {
  tail <- QTLTOOLS_TAIL[[type]]
  no_se <- setdiff(tail, "slope_se")
  list(
    list(cols = c(QTLTOOLS_HEAD, tail), grouped = FALSE),
    list(cols = c(QTLTOOLS_HEAD, no_se), grouped = FALSE),
    list(cols = c(QTLTOOLS_HEAD, "best_phe_id", "n_phe_in_grp", tail), grouped = TRUE),
    list(cols = c(QTLTOOLS_HEAD, "best_phe_id", "n_phe_in_grp", no_se), grouped = TRUE)
  )
}

#' Read QTLtools cis output (one or more files, plain or gzipped) with named columns.
read_qtltools <- function(paths, type = c("permute", "nominal")) {
  type <- match.arg(type)
  out_cols <- QTLTOOLS_OUT[[type]]
  layouts <- qtltools_layouts(type)
  pieces <- lapply(paths, function(p) {
    if (!file.exists(p) || file.size(p) == 0) return(NULL)
    # QTLtools writes space-separated output; the separator is fixed because
    # phenotype and variant IDs can contain ";" (e.g. "chr1:100-100;+").
    # A gzipped empty file has a non-zero size, so also catch fread's "File is empty".
    dt <- tryCatch(fread(p, header = FALSE, sep = " "), error = function(e) {
      if (grepl("empty", conditionMessage(e), ignore.case = TRUE)) return(NULL)
      stop(e)
    })
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt <- dt[V1 != "phe_id" & V1 != "grp_id"]  # header lines written by chunk 0
    layout <- Filter(function(l) length(l$cols) == ncol(dt), layouts)
    if (!length(layout)) {
      stop(sprintf("%s: %d columns, not a known QTLtools %s layout (expected %s)", p, ncol(dt), type,
                   paste(unique(vapply(layouts, function(l) length(l$cols), 0L)), collapse = "/")))
    }
    layout <- layout[[1]]
    setnames(dt, layout$cols)
    if (layout$grouped) {
      setnames(dt, "phe_id", "grp_id")
      setnames(dt, "best_phe_id", "phe_id")
    } else {
      dt[, `:=`(grp_id = NA_character_, n_phe_in_grp = NA_real_)]
    }
    if (!"slope_se" %in% names(dt)) dt[, slope_se := NA_real_]
    setcolorder(dt, out_cols)
    dt
  })
  out <- rbindlist(pieces)
  if (nrow(out) == 0) {
    out <- as.data.table(setNames(replicate(length(out_cols), character(0), simplify = FALSE), out_cols))
  }
  chr_cols <- c("phe_id", "grp_id", "phe_chr", "phe_strd", "var_id", "var_chr")
  num_cols <- setdiff(out_cols, chr_cols)
  out[, (chr_cols) := lapply(.SD, as.character), .SDcols = chr_cols]
  out[, (num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = num_cols]
  out[]
}
