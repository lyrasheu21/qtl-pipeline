# Input/output helpers shared by every pipeline script.

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Send stdout and messages to the rule's log file.
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
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Reading parquet needs the 'arrow' package: ", path)
    }
    return(as.data.table(arrow::read_parquet(path)))
  }
  fread(path, ...)
}

#' Read the CCI score matrix as a numeric matrix (rows = CCIs, columns = donors).
read_cci_matrix <- function(path, id_column) {
  dt <- read_table_any(path)
  if (!id_column %in% names(dt)) {
    stop("CCI ID column '", id_column, "' not found in ", path)
  }
  ids <- as.character(dt[[id_column]])
  if (anyDuplicated(ids)) {
    stop("Duplicated CCI IDs in ", path, ": ", paste(head(unique(ids[duplicated(ids)])), collapse = ", "))
  }
  dt[, (id_column) := NULL]
  non_numeric <- names(dt)[!vapply(dt, is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric donor columns in CCI matrix: ", paste(head(non_numeric), collapse = ", "))
  }
  mat <- as.matrix(dt)
  rownames(mat) <- ids
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

#' Fill `{subclass}` (and any other named placeholder) in a path pattern.
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
