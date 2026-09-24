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

test_that("build_covariates treats missing values and unknown covariates as errors", {
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1", "d4"), defs, "Sex", 0), "Missing covariate")
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1"), defs, "Age", 0), "not defined")
  expect_error(build_covariates(meta, "SubID", pcs, "SubID", c("d1"), defs, "Sex", 5), "PC columns not found")
})

test_that("assemble_covariates drops donors with missing values instead of failing", {
  f <- matrix(c(1, 2, 3), ncol = 1, dimnames = list(c("d1", "d2", "d3"), NULL))
  blocks <- list(metadata_covariates(meta, "SubID", c("d1", "d4", "d3"), defs, "Sex"),
                 factor_covariates(f, c("d1", "d4", "d3"), 1, "PEER_"))
  res <- assemble_covariates(blocks, c("d1", "d4", "d3"))
  expect_equal(res$donors, c("d1", "d3"))
  expect_equal(res$dropped, "d4")
  expect_equal(res$cov$id, c("Sex", "PEER_1"))
})

test_that("asking for more factors than exist is an error", {
  f <- matrix(1:6, ncol = 2, dimnames = list(c("a", "b", "c"), NULL))
  expect_error(factor_covariates(f, c("a", "b"), 3, "PEER_"), "only 2")
})

test_that("phenotype PCs have one row per donor", {
  set.seed(1)
  mat <- matrix(rnorm(200), nrow = 20, dimnames = list(NULL, paste0("d", 1:10)))
  pc <- phenotype_pcs(mat, 3)
  expect_equal(dim(pc), c(10L, 3L))
  expect_equal(rownames(pc), paste0("d", 1:10))
})
