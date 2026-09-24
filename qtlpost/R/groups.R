# Group lookups. A "group" is whatever q-values and plots are stratified by:
# the anchor cell subclass for CCIs, a tissue or cell type for CAGE / eQTL.
# The optional group table (e.g. config/cell_types.tsv) maps each group to a
# class and superclass and sets plotting order and colours. Without a table,
# class and superclass are the group itself.

#' Read the group table. Accepts `group` or (for cell types) `subclass` as the
#' first key column; `class`, `superclass`, `color`, `class_color` are optional.
read_group_table <- function(path = NULL) {
  empty <- data.table(group = character(), class = character(), superclass = character(),
                      color = character(), class_color = character())
  if (is.null(path) || !nzchar(path)) return(empty)
  gt <- fread(path)
  if (!"group" %in% names(gt)) {
    if (!"subclass" %in% names(gt)) stop("Group table needs a 'group' (or 'subclass') column: ", path)
    setnames(gt, "subclass", "group")
  }
  if (anyDuplicated(gt$group)) stop("Duplicated groups in ", path)
  if (!"class" %in% names(gt)) gt[, class := group]
  if (!"superclass" %in% names(gt)) gt[, superclass := class]
  if (!"color" %in% names(gt)) gt[, color := NA_character_]
  if (!"class_color" %in% names(gt)) gt[, class_color := NA_character_]
  gt[, .(group = as.character(group), class = as.character(class), superclass = as.character(superclass),
         color, class_color)]
}

#' Map groups to their class or superclass. Groups missing from a non-empty
#' table keep their own name as class and become "Other" as superclass.
map_group <- function(group, gt, to = c("class", "superclass")) {
  to <- match.arg(to)
  if (nrow(gt) == 0) return(as.character(group))
  out <- gt[[to]][match(group, gt$group)]
  unknown <- is.na(out) & !is.na(group)
  if (any(unknown)) {
    warning("Groups not in the group table: ", paste(unique(group[unknown]), collapse = ", "), call. = FALSE)
    out[unknown] <- if (to == "class") group[unknown] else "Other"
  }
  out
}

group_order <- function(gt, observed = character()) union(gt$group, sort(unique(observed)))
class_order <- function(gt, observed = character()) union(unique(gt$class), sort(unique(observed)))
superclass_order <- function(gt, observed = character()) union(unique(gt$superclass), sort(unique(observed)))

#' Named colour vector for groups (or classes); falls back to a generated
#' palette for anything without a colour.
group_colors <- function(gt, groups = character(), level = c("group", "class")) {
  level <- match.arg(level)
  key <- if (level == "group") "group" else "class"
  col <- if (level == "group") "color" else "class_color"
  x <- gt[!duplicated(get(key))]
  cols <- setNames(x[[col]], x[[key]])
  cols <- cols[!is.na(cols)]
  missing <- setdiff(unique(c(groups)), names(cols))
  if (length(missing)) {
    fill <- grDevices::hcl.colors(max(length(missing), 2), "Dark 3")[seq_along(missing)]
    cols <- c(cols, setNames(fill, missing))
  }
  cols
}
