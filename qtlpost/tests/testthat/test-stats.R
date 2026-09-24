test_that("q-values are computed within each subclass, not pooled", {
  skip_if_not_installed("qvalue")
  set.seed(1)
  # Subclass A: strong signal; subclass B: nulls only.
  p <- c(runif(200)^8, runif(200))
  g <- rep(c("A", "B"), each = 200)
  q_grouped <- qvalue_by_group(p, g)
  q_pooled <- qvalue::qvalue(p)$qvalues
  expect_equal(q_grouped[g == "B"], qvalue::qvalue(p[g == "B"])$qvalues)
  expect_false(isTRUE(all.equal(q_grouped, q_pooled)))
})

test_that("tiny groups fall back to Benjamini-Hochberg instead of failing", {
  skip_if_not_installed("qvalue")
  q <- qvalue_by_group(c(0.001, 0.2, 0.9), c("A", "A", "A"))
  expect_equal(q, p.adjust(c(0.001, 0.2, 0.9), "BH"), tolerance = 1e-8)
})

test_that("NA p-values stay NA", {
  skip_if_not_installed("qvalue")
  q <- qvalue_by_group(c(0.01, NA, 0.5), c("A", "A", "A"))
  expect_true(is.na(q[2]))
})

test_that("lambda is about 1 under the null", {
  set.seed(2)
  expect_equal(genomic_lambda(runif(1e5)), 1, tolerance = 0.03)
})

test_that("Fisher's method matches a direct calculation", {
  expect_equal(fisher_combine(0.01, 0.02), pchisq(-2 * log(0.01 * 0.02), 4, lower.tail = FALSE))
})

test_that("the base-R Storey fallback agrees with qvalue where qvalue works", {
  skip_if_not_installed("qvalue")
  set.seed(4)
  p <- c(runif(300)^4, runif(700))
  expect_equal(storey_qvalue(p)$qvalues, qvalue::qvalue(p)$qvalues, tolerance = 1e-6)
})

test_that("skewed p-values that break qvalue() still get Storey q-values", {
  skip_if_not_installed("qvalue")
  set.seed(1)
  p <- runif(200)^8
  expect_message(q <- qvalue_by_group(p, rep("A", 200)), "base-R Storey")
  expect_equal(q, storey_qvalue(p)$qvalues)
})

test_that("method = 'bh' gives Benjamini-Hochberg within each group", {
  set.seed(5)
  p <- c(runif(50)^3, runif(50))
  g <- rep(c("A", "B"), each = 50)
  q <- qvalue_by_group(p, g, method = "bh")
  expect_equal(q[g == "A"], p.adjust(p[g == "A"], "BH"))
  expect_equal(q[g == "B"], p.adjust(p[g == "B"], "BH"))
})
