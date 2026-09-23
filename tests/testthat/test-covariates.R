meta <- data.table(SubID = c("d3", "d1", "d2"), PMI = c(10, 20, 30), Sex = c("M", "F", "M"),
                   counts = c(100, 1000, 10000))
pcs <- data.table(SubID = c("d1", "d2", "d3"), PC1 = c(0.1, 0.2, 0.3), PC2 = c(1, 2, 3))
defs <- list(scaled_PMI = list(column = "PMI", transform = "scale"),
             Sex = list(column = "Sex", transform = "factor"),
             log_counts = list(column = "counts", transform = "log"))

test_that("covariates follow donor order and apply transforms", {
  cov <- build_covariates(meta, "SubID", pcs, "SubID", donors = c("d1", "d2", "d3"),
                          definitions = defs, use = c("scaled_PMI", "Sex", "log_counts"), n_pcs = 2)
  expect_equal(names(cov), c("id", "d1", "d2", "d3"))
  expect_equal(cov$id, c("scaled_PMI", "Sex", "log_counts", "ancestry_PC1", "ancestry_PC2"))
  expect_equal(unlist(cov[id == "Sex", -1]), c(d1 = "F", d2 = "M", d3 = "M"))
  expect_equal(as.numeric(unlist(cov[id == "log_counts", -1])), log(c(1000, 10000, 100)))
  expect_equal(as.numeric(unlist(cov[id == "ancestry_PC2", -1])), c(1, 2, 3))
})

test_that("only the requested covariates are used", {
  cov <- build_covariates(meta, "SubID", pcs, "SubID", c("d1", "d2", "d3"), defs, use = "Sex", n_pcs = 0)
  expect_equal(cov$id, "Sex")
})

test_that("missing values and unknown covariates are errors", {
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1", "d4"), defs, "Sex", 0), "Missing covariate")
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1"), defs, "Age", 0), "not defined")
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1"), defs, "Sex", 5), "PC columns not found")
})
