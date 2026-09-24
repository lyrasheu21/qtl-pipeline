test_that("current and legacy CCI ID formats parse the same way", {
  x <- parse_cci_ids(c("Astro:APOE->Micro:TREM2", "L_Astro:APOE^Micro:TREM2", "R_VLMC:COL4A1^Micro:TREM2"))
  expect_equal(x$sending_group, c("Astro", "Astro", "VLMC"))
  expect_equal(x$ligand, c("APOE", "APOE", "COL4A1"))
  expect_equal(x$receiving_group, rep("Micro", 3))
  expect_equal(x$receptor, rep("TREM2", 3))
})

test_that("subclass names with underscores and complex receptors are kept intact", {
  x <- parse_cci_ids("EN_L3_5_IT_2:NRG1->IN_PVALB_CHC:ITGA1_ITGB1")
  expect_equal(x$sending_group, "EN_L3_5_IT_2")
  expect_equal(x$receiving_group, "IN_PVALB_CHC")
  expect_equal(x$receptor, "ITGA1_ITGB1")
})

test_that("gene names containing a hyphen are not split", {
  expect_equal(parse_cci_ids("Micro:HLA-A->Astro:LILRB2")$ligand, "HLA-A")
})

test_that("malformed IDs fail loudly", {
  expect_error(parse_cci_ids("Astro:APOE"), "Unrecognised")
  expect_error(parse_cci_ids("APOE->TREM2"), "Subclass:Gene")
})

test_that("empty input gives an empty table", {
  expect_equal(nrow(parse_cci_ids(character(0))), 0)
})

test_that("the anchor follows the direction", {
  s <- cci_annotation("OPC:PDGFA->Astro:PDGFRA", "sending")
  r <- cci_annotation("OPC:PDGFA->Astro:PDGFRA", "receiving")
  expect_equal(c(s$group, s$gene, s$partner_group, s$partner_gene), c("OPC", "PDGFA", "Astro", "PDGFRA"))
  expect_equal(c(r$group, r$gene, r$partner_group, r$partner_gene), c("Astro", "PDGFRA", "OPC", "PDGFA"))
  expect_equal(s$phenotype_id, "OPC:PDGFA->Astro:PDGFRA")
})
