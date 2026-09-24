fake_permute <- function() {
  set.seed(3)
  n <- 60
  data.table(phe_id = paste0("p", 1:n), grp_id = NA_character_, n_phe_in_grp = NA_real_,
             phe_chr = "1", phe_from = 1, phe_to = 2, n_var_in_cis = 10, var_id = paste0("rs", 1:n),
             var_chr = "1", var_from = 1, dist_phe_var = 0, slope = 0.1, slope_se = 0.05,
             nom_pval = 0.01, r_squared = 0.1, adj_emp_pval = 0.1,
             adj_beta_pval = c(runif(30)^6, runif(30)))
}
fake_annotation <- function() {
  data.table(phenotype_id = paste0("p", 1:60), set = "s", chrom = "1", start = 1, end = 2, strand = "+",
             gene = paste0("G", 1:60), group = rep(c("Astro", "Micro"), each = 30))
}

test_that("fdr_group_by = none gives one global q-value; group stratifies", {
  skip_if_not_installed("qvalue")
  gt <- toy_group_table()
  global <- annotate_permute(fake_permute(), fake_annotation(), gt, "s", fdr_group_by = "none")
  by_group <- annotate_permute(fake_permute(), fake_annotation(), gt, "s", fdr_group_by = "group")
  setorder(global, phenotype_id); setorder(by_group, phenotype_id)
  expect_equal(global$qval, qvalue_by_group(global$adj_beta_pval, rep(1, 60)))
  expect_false(isTRUE(all.equal(global$qval, by_group$qval)))
  expect_equal(unique(by_group$superclass), "Glia")
})

test_that("phenotypes missing from the annotation get group 'all'", {
  skip_if_not_installed("qvalue")
  ann <- fake_annotation()[1:50]
  out <- annotate_permute(fake_permute(), ann, read_group_table(NULL), "s", fdr_group_by = "group")
  expect_equal(out[phenotype_id == "p55", group], "all")
})

test_that("pi0 diagnostics flag extremes", {
  expect_match(pi0_diagnostic(0.99), "weak signal")
  expect_match(pi0_diagnostic(0.01), "under-correction")
  expect_equal(pi0_diagnostic(0.5), "ok")
})
