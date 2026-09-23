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

test_that("small groups fall back to pi0 = 1 instead of failing", {
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

test_that("partial eta-squared is between 0 and 1", {
  set.seed(3)
  d <- data.frame(x = rnorm(100), g = rbinom(100, 2, 0.4))
  d$y <- d$x * d$g + rnorm(100)
  e <- partial_eta_sq(lm(y ~ x * g, data = d), "x:g")
  expect_true(e > 0 && e < 1)
})
