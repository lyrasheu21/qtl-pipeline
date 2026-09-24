spec <- list(name = "T", snp_col = "SNP", bp_col = "BP", beta_col = "b", se_col = "se", p_col = "p",
             n_col = NULL, n = 1000, type = "cc", s = 0.3)

test_that("GWAS columns are standardised and bad rows removed", {
  g <- data.table(SNP = c("rs1", "rs2", "rs2", "rs3"), BP = 1:4, b = c(0.1, 0.2, 0.2, 0.3),
                  se = c(0.01, 0.02, 0.02, 0), p = c(0.01, 0.02, 0.02, 0.03))
  out <- standardize_gwas(g, spec)
  expect_equal(out$snp, c("rs1", "rs2"))   # duplicate and se = 0 dropped
  expect_equal(out$n, c(1000, 1000))
})

test_that("case-control GWAS need a case proportion", {
  g <- standardize_gwas(data.table(SNP = "rs1", BP = 1, b = 0.1, se = 0.01, p = 0.01), spec)
  expect_equal(gwas_coloc_dataset(g, spec)$s, 0.3)
  bad <- spec; bad$s <- NULL
  expect_error(gwas_coloc_dataset(g, bad), "case proportion")
})

test_that("QTL SE comes from slope_se, or from the p-value when absent", {
  nom <- data.table(var_id = c("rs1", "rs2"), var_from = c(10, 20), slope = c(0.5, -0.2),
                    slope_se = c(0.1, 0.05), nom_pval = c(1e-6, 0.01))
  ds <- qtl_coloc_dataset(nom, 100)
  expect_equal(unname(ds$varbeta), c(0.01, 0.0025))
  nom$slope_se <- NA_real_
  ds2 <- qtl_coloc_dataset(nom, 100)
  t2 <- abs(qt(0.01 / 2, df = 99))
  expect_equal(unname(ds2$varbeta[2]), (0.2 / t2)^2)
})

test_that("coloc tiers follow the original rules", {
  expect_equal(classify_coloc_tier(c(0.8, 0.15)), "Gold_Unique")
  expect_equal(classify_coloc_tier(c(0.4, 0.35)), "Gold_LD_Pair")
  expect_equal(classify_coloc_tier(c(0.2, 0.15)), "Silver_Multiple")
  expect_equal(classify_coloc_tier(0.2), "Bronze_Single")
  expect_equal(classify_coloc_tier(numeric(0)), "Diffuse")
})
