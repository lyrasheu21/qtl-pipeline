test_that("MAGMA 1.08+ .gsa.out is parsed, using FULL_NAME when present", {
  path <- tempfile()
  writeLines(c("# MEAN_SAMPLE_SIZE = 1000",
               "# TOTAL_GENES = 18000",
               "VARIABLE TYPE NGENES BETA BETA_STD SE P FULL_NAME",
               "sending_Astr SET 40 0.21 0.012 0.08 0.004 sending_Astro",
               "receiving_EN SET 12 -0.05 -0.002 0.1 0.6 receiving_EN_L2_3_IT"), path)
  x <- read_magma_gsa(path)
  expect_equal(x$set, c("sending_Astro", "receiving_EN_L2_3_IT"))
  expect_equal(x$p, c(0.004, 0.6))
})

test_that("older .sets.out layout is parsed", {
  path <- tempfile()
  writeLines(c("# comment", "SET NGENES BETA BETA_STD SE P", "sending_Astro 40 0.21 0.012 0.08 0.004"), path)
  expect_equal(read_magma_gsa(path)$set, "sending_Astro")
})

test_that("an empty result file gives an empty table", {
  path <- tempfile(); file.create(path)
  expect_equal(nrow(read_magma_gsa(path)), 0)
})
