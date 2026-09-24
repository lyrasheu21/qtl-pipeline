# Input/output helpers.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Send stdout and messages to a log file (used by the pipeline scripts).
start_logging <- function(log_path) {
  if (length(log_path) == 0 || !nzchar(log_path)) return(invisible(NULL))
  con <- file(log_path, open = "wt")
  sink(con)
  sink(con, type = "message")
  invisible(con)
}

#' Read a TSV/CSV (optionally gzipped) or parquet file into a data.table.
read_table_any <- function(path, ...) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("Reading parquet needs the 'arrow' package: ", path)
    return(as.data.table(arrow::read_parquet(path)))
  }
  fread(path, ...)
}

#' Read a wide feature x sample table as a numeric matrix (rows = features).
#' @param id_column column holding the feature (phenotype) IDs.
read_feature_matrix <- function(path, id_column) {
  dt <- read_table_any(path)
  if (!id_column %in% names(dt)) stop("ID column '", id_column, "' not found in ", path)
  ids <- as.character(dt[[id_column]])
  if (anyDuplicated(ids)) {
    stop("Duplicated IDs in ", path, ": ", paste(head(unique(ids[duplicated(ids)])), collapse = ", "))
  }
  dt[, (id_column) := NULL]
  non_numeric <- names(dt)[!vapply(dt, is.numeric, logical(1))]
  if (length(non_numeric)) stop("Non-numeric sample columns in ", path, ": ", paste(head(non_numeric), collapse = ", "))
  mat <- as.matrix(dt)
  rownames(mat) <- ids
  mat
}

#' Read a tab-separated table written by the pipeline (or documented as TSV).
#' The separator is fixed because IDs can contain ";" or ",", which would
#' otherwise confuse fread's separator detection.
read_tsv <- function(path, ...) fread(path, sep = "\t", ...)

#' Read the CCI score matrix (rows = CCIs, columns = donors).
read_cci_matrix <- function(path, id_column) read_feature_matrix(path, id_column)

#' Read a QTLtools phenotype BED (plain or bgzipped) into a list with the six
#' annotation columns and the numeric sample matrix.
read_bed_matrix <- function(path) {
  dt <- fread(path, sep = "\t")
  setnames(dt, 1:6, c("chr", "start", "end", "pid", "gid", "strand"))
  mat <- as.matrix(dt[, -(1:6), with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- dt$pid
  list(info = dt[, 1:6], mat = mat)
}

#' Rename sample columns with a sample map and average replicates.
#'
#' @param mat numeric matrix, samples in columns.
#' @param map data.frame with the assay sample ID in `from` and the donor ID in `to`.
#' @return matrix with one column per donor; unmapped samples are dropped.
apply_sample_map <- function(mat, map, from, to) {
  map <- as.data.table(map)
  donor <- as.character(map[[to]])[match(colnames(mat), as.character(map[[from]]))]
  keep <- !is.na(donor)
  if (!any(keep)) stop("No sample IDs matched the '", from, "' column of the sample map")
  mat <- mat[, keep, drop = FALSE]
  donor <- donor[keep]
  if (anyDuplicated(donor)) {
    message(sprintf("Averaging replicates: %d samples -> %d donors", length(donor), length(unique(donor))))
    groups <- split(seq_along(donor), factor(donor, levels = unique(donor)))
    avg <- vapply(groups, function(idx) rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE), numeric(nrow(mat)))
    mat <- matrix(avg, nrow = nrow(mat), dimnames = list(rownames(mat), names(groups)))
  } else {
    colnames(mat) <- donor
  }
  mat
}

#' Read a per-sample factor table (genotype PCs or hidden factors) and return
#' a matrix with samples as rows and factors as columns.
#'
#' Accepts TSV/CSV/whitespace tables in either orientation (samples as rows,
#' or factors as rows as in some `.pca` files), plink2 `.eigenvec`, or an
#' `.RData` / `.rda` file containing a matrix `object`. The orientation is
#' chosen so that the sample IDs match `samples` best.
read_sample_factors <- function(path, samples, object = NULL, id_column = NULL) {
  if (grepl("\\.(rdata|rda)$", path, ignore.case = TRUE)) {
    env <- new.env()
    load(path, envir = env)
    obj <- object %||% ls(env)[1]
    if (!exists(obj, envir = env)) stop("Object '", obj, "' not found in ", path)
    mat <- as.matrix(get(obj, envir = env))
  } else {
    dt <- fread(path)
    id_col <- id_column %||% names(dt)[1]
    if (identical(id_col, "#FID") && "IID" %in% names(dt)) id_col <- "IID"  # plink2 eigenvec
    if (!id_col %in% names(dt)) stop("ID column '", id_col, "' not found in ", path)
    ids <- as.character(dt[[id_col]])
    num <- dt[, vapply(dt, is.numeric, logical(1)), with = FALSE]
    mat <- as.matrix(num)
    rownames(mat) <- ids
  }
  mode(mat) <- "numeric"
  rows_hit <- sum(rownames(mat) %in% samples)
  cols_hit <- sum(colnames(mat) %in% samples)
  if (cols_hit > rows_hit) mat <- t(mat)
  if (is.null(colnames(mat))) colnames(mat) <- paste0("F", seq_len(ncol(mat)))
  mat
}

#' Read the gene annotation table and add `chrom_norm` (chromosome without a
#' "chr" prefix, matching plink-written VCFs).
read_annotation <- function(path) {
  annot <- fread(path)
  required <- c("gene_name", "gene_id", "chrom", "start", "end", "strand")
  missing <- setdiff(required, names(annot))
  if (length(missing)) stop("Gene annotation is missing columns: ", paste(missing, collapse = ", "))
  annot[, chrom_norm := sub("^chr", "", as.character(chrom), ignore.case = TRUE)]
  annot[!duplicated(gene_name)]
}

#' Read a residualised expression file: first column gene symbol, then donors.
read_expression <- function(path) {
  dt <- fread(path)
  genes <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  mat[!duplicated(genes), , drop = FALSE]
}

#' Fill `{name}` placeholders in a path pattern, e.g. fill_pattern(p, group = "Astro").
fill_pattern <- function(pattern, ...) {
  values <- list(...)
  for (key in names(values)) {
    pattern <- gsub(paste0("{", key, "}"), values[[key]], pattern, fixed = TRUE)
  }
  pattern
}

#' Write a data.table as TSV; gzip when the path ends in .gz.
write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}

#' Read a text file of IDs, one per line, ignoring blank lines.
read_ids <- function(path) {
  x <- readLines(path, warn = FALSE)
  x[nzchar(x)]
}
