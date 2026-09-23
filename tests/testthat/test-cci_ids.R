test_that("current and legacy CCI ID formats parse the same way", {
  x <- parse_cci_ids(c("Astro:APOE->Micro:TREM2", "L_Astro:APOE^Micro:TREM2", "R_VLMC:COL4A1^Micro:TREM2"))
  expect_equal(x$sending_subclass, c("Astro", "Astro", "VLMC"))
  expect_equal(x$ligand, c("APOE", "APOE", "COL4A1"))
  expect_equal(x$receiving_subclass, rep("Micro", 3))
  expect_equal(x$receptor, rep("TREM2", 3))
})

test_that("subclass names with underscores and complex receptors are kept intact", {
  x <- parse_cci_ids("EN_L3_5_IT_2:NRG1->IN_PVALB_CHC:ITGA1_ITGB1")
  expect_equal(x$sending_subclass, "EN_L3_5_IT_2")
  expect_equal(x$receiving_subclass, "IN_PVALB_CHC")
  expect_equal(x$receptor, "ITGA1_ITGB1")
})

test_that("gene names containing a hyphen are not split", {
  x <- parse_cci_ids("Micro:HLA-A->Astro:LILRB2")
  expect_equal(x$ligand, "HLA-A")
})

test_that("malformed IDs fail loudly", {
  expect_error(parse_cci_ids("Astro:APOE"), "Unrecognised")
  expect_error(parse_cci_ids("APOE->TREM2"), "Subclass:Gene")
})

test_that("empty input gives an empty table", {
  expect_equal(nrow(parse_cci_ids(character(0))), 0)
})

test_that("anchor follows the direction", {
  ct <- toy_cell_types()
  s <- annotate_cci("OPC:PDGFA->Astro:PDGFRA", "sending", ct)
  r <- annotate_cci("OPC:PDGFA->Astro:PDGFRA", "receiving", ct)
  expect_equal(c(s$subclass, s$gene, s$partner_gene), c("OPC", "PDGFA", "PDGFRA"))
  expect_equal(c(r$subclass, r$gene, r$partner_gene), c("Astro", "PDGFRA", "PDGFA"))
  expect_equal(r$superclass, "Glia")
  expect_equal(s$cell_pair, "OPC->Astro")
})

test_that("unknown subclasses are flagged, not dropped", {
  expect_warning(x <- annotate_cci("Unknown:A->Astro:B", "sending", toy_cell_types()), "not in cell_types")
  expect_equal(x$class, "Unknown")
  expect_equal(x$superclass, "Other")
})
