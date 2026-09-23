# Cell-type lookups: subclass -> class -> superclass, plotting order and colours.
# The table lives in config/cell_types.tsv; row order is the plotting order.

read_cell_types <- function(path) {
  ct <- fread(path)
  required <- c("subclass", "class", "superclass", "color", "class_color")
  missing <- setdiff(required, names(ct))
  if (length(missing)) stop("cell_types table is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(ct$subclass)) stop("Duplicated subclasses in cell_types table")
  ct
}

#' Map subclasses to a column of the cell-type table. Unknown subclasses keep
#' their own name for `class` and become "Other" for `superclass`.
map_subclass <- function(subclass, ct, to = c("class", "superclass")) {
  to <- match.arg(to)
  out <- ct[[to]][match(subclass, ct$subclass)]
  unknown <- is.na(out) & !is.na(subclass)
  if (any(unknown)) {
    warning("Subclasses not in cell_types table: ",
            paste(unique(subclass[unknown]), collapse = ", "), call. = FALSE)
    out[unknown] <- if (to == "class") subclass[unknown] else "Other"
  }
  out
}

subclass_order <- function(ct) ct$subclass
class_order <- function(ct) unique(ct$class)
superclass_order <- function(ct) unique(ct$superclass)

subclass_colors <- function(ct) setNames(ct$color, ct$subclass)

class_colors <- function(ct) {
  x <- ct[!duplicated(class)]
  setNames(x$class_color, x$class)
}
