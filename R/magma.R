# MAGMA output parsing.

#' Read a MAGMA gene-set analysis result (.gsa.out, or .sets.out from older
#' versions). Returns columns set, ngenes, beta, beta_std, se, p.
read_magma_gsa <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) {
    return(data.table(set = character(), ngenes = numeric(), beta = numeric(),
                      beta_std = numeric(), se = numeric(), p = numeric()))
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  dt <- fread(text = paste(lines, collapse = "\n"), header = TRUE)
  name_col <- intersect(c("FULL_NAME", "VARIABLE", "SET"), names(dt))[1]
  if (is.na(name_col)) stop("Unrecognised MAGMA output: ", path)
  data.table(set = as.character(dt[[name_col]]), ngenes = dt$NGENES, beta = dt$BETA,
             beta_std = dt$BETA_STD, se = dt$SE, p = dt$P)
}
