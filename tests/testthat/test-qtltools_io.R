write_rows <- function(n_cols, n_rows = 2, header = FALSE) {
  path <- tempfile(fileext = ".txt")
  rows <- vapply(seq_len(n_rows), function(i) {
    x <- c(sprintf("A:L%d->B:R", i), "1", rep(as.character(i), n_cols - 2))
    x[8] <- sprintf("rs%d", i)
    x[9] <- "1"
    paste(x, collapse = " ")
  }, "")
  if (header) rows <- c("phe_id x", rows)
  writeLines(rows, path)
  path
}

test_that("permute output with --std-err (21 columns) is named", {
  dt <- read_qtltools(write_rows(21), "permute")
  expect_equal(names(dt), QTLTOOLS_COLS$permute)
  expect_equal(dt$var_id, c("rs1", "rs2"))
  expect_type(dt$adj_beta_pval, "double")
})

test_that("legacy permute output without --std-err (20 columns) is read", {
  dt <- read_qtltools(write_rows(20), "permute")
  expect_true(all(is.na(dt$slope_se)))
  expect_equal(dt$adj_beta_pval, c(1, 2))
})

test_that("nominal output and header lines", {
  dt <- read_qtltools(write_rows(16, header = TRUE), "nominal")
  expect_equal(nrow(dt), 2)
  expect_equal(names(dt), QTLTOOLS_COLS$nominal)
})

test_that("empty files give an empty, correctly named table", {
  path <- tempfile(); file.create(path)
  dt <- read_qtltools(path, "nominal")
  expect_equal(nrow(dt), 0)
  expect_equal(names(dt), QTLTOOLS_COLS$nominal)
})

test_that("unexpected column counts are rejected", {
  expect_error(read_qtltools(write_rows(12), "permute"), "expected")
})
