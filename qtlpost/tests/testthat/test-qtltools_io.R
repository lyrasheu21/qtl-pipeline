write_rows <- function(n_cols, grouped = FALSE, header = FALSE, n_rows = 2) {
  path <- tempfile(fileext = ".txt")
  rows <- vapply(seq_len(n_rows), function(i) {
    x <- rep(as.character(i), n_cols)
    x[1] <- if (grouped) sprintf("GENE%d", i) else sprintf("A:L%d->B:R", i)
    x[2] <- "1"; x[5] <- "+"
    off <- if (grouped) 2 else 0
    if (grouped) { x[6] <- sprintf("peak%d", i); x[7] <- "3" }
    x[8 + off] <- sprintf("rs%d", i)
    x[9 + off] <- "1"
    paste(x, collapse = " ")
  }, "")
  if (header) rows <- c("phe_id x", rows)
  writeLines(rows, path)
  path
}

test_that("permute output with --std-err (21 columns) is named", {
  dt <- read_qtltools(write_rows(21), "permute")
  expect_equal(names(dt), QTLTOOLS_OUT$permute)
  expect_equal(dt$var_id, c("rs1", "rs2"))
  expect_true(all(is.na(dt$grp_id)))
  expect_type(dt$adj_beta_pval, "double")
})

test_that("permute output without --std-err (20 columns) is read", {
  dt <- read_qtltools(write_rows(20), "permute")
  expect_true(all(is.na(dt$slope_se)))
  expect_equal(dt$adj_beta_pval, c(1, 2))
})

test_that("--grp-best output puts the group in grp_id and the best phenotype in phe_id", {
  dt <- read_qtltools(write_rows(23, grouped = TRUE), "permute")
  expect_equal(dt$grp_id, c("GENE1", "GENE2"))
  expect_equal(dt$phe_id, c("peak1", "peak2"))
  expect_equal(dt$n_phe_in_grp, c(3, 3))
  expect_equal(dt$var_id, c("rs1", "rs2"))
  nom <- read_qtltools(write_rows(18, grouped = TRUE), "nominal")
  expect_equal(nom$phe_id, c("peak1", "peak2"))
})

test_that("nominal output and header lines", {
  dt <- read_qtltools(write_rows(16, header = TRUE), "nominal")
  expect_equal(nrow(dt), 2)
  expect_equal(names(dt), QTLTOOLS_OUT$nominal)
})

test_that("empty files give an empty, correctly named table", {
  path <- tempfile(); file.create(path)
  dt <- read_qtltools(path, "nominal")
  expect_equal(nrow(dt), 0)
  expect_equal(names(dt), QTLTOOLS_OUT$nominal)
})

test_that("unexpected column counts are rejected", {
  expect_error(read_qtltools(write_rows(12), "permute"), "not a known")
})
