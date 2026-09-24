test_that("MAGMA 1.08+ .gsa.out is parsed, using FULL_NAME when present", {
  path <- tempfile()
  writeLines(c("# MEAN_SAMPLE_SIZE = 1000",
               "VARIABLE TYPE NGENES BETA BETA_STD SE P FULL_NAME",
               "cci_sending__As SET 40 0.21 0.012 0.08 0.004 cci_sending__Astro",
               "cci_receiving__ SET 12 -0.05 -0.002 0.1 0.6 cci_receiving__EN_L2_3_IT"), path)
  x <- read_magma_gsa(path)
  expect_equal(x$set, c("cci_sending__Astro", "cci_receiving__EN_L2_3_IT"))
  expect_equal(x$p, c(0.004, 0.6))
})

test_that("older .sets.out layout is parsed", {
  path <- tempfile()
  writeLines(c("# comment", "SET NGENES BETA BETA_STD SE P", "cage__DLPFC 40 0.21 0.012 0.08 0.004"), path)
  expect_equal(read_magma_gsa(path)$set, "cage__DLPFC")
})

test_that("an empty result file gives an empty table", {
  path <- tempfile(); file.create(path)
  expect_equal(nrow(read_magma_gsa(path)), 0)
})

test_that("gene sets map symbols to IDs and name sets by phenotype set and level", {
  qt <- data.table(set = c("cage", "cage", "cage"), gene = c("A", "B", "C"), group = "DLPFC",
                   fdr_pass = c(TRUE, TRUE, FALSE), adj_beta_pval = c(0.001, 0.01, 0.04))
  annot <- data.table(gene_name = c("A", "B", "C"), gene_id = c("ENSG1", "ENSG2", "ENSG3"))
  fdr <- magma_gene_sets(qt, annot, "fdr", "group")
  expect_equal(fdr$gene_id, c("ENSG1", "ENSG2"))
  expect_equal(unique(fdr$set), "cage__DLPFC")
  expect_equal(nrow(magma_gene_sets(qt, annot, "nominal", "group")), 3)
})
