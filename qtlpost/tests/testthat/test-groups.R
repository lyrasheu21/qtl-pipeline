test_that("groups map to class and superclass", {
  gt <- toy_group_table()
  expect_equal(map_group(c("EN_L2_3_IT", "Micro"), gt, "class"), c("EN", "Immune"))
  expect_equal(map_group("Micro", gt, "superclass"), "Glia")
})

test_that("unknown groups are flagged, not dropped", {
  expect_warning(x <- map_group("Unknown", toy_group_table(), "superclass"), "not in the group table")
  expect_equal(x, "Other")
})

test_that("without a group table, class and superclass are the group", {
  gt <- read_group_table(NULL)
  expect_equal(map_group("DLPFC", gt, "class"), "DLPFC")
  expect_equal(map_group("DLPFC", gt, "superclass"), "DLPFC")
})

test_that("a cell_types table keyed by 'subclass' is accepted", {
  path <- tempfile(fileext = ".tsv")
  fwrite(data.table(subclass = c("A", "B"), class = c("X", "X")), path, sep = "\t")
  gt <- read_group_table(path)
  expect_equal(gt$group, c("A", "B"))
  expect_equal(gt$superclass, c("X", "X"))
})

test_that("colours are filled in for groups without one", {
  cols <- group_colors(toy_group_table(), c("Astro", "NewGroup"))
  expect_equal(cols[["Astro"]], "#C19A6B")
  expect_true(!is.na(cols[["NewGroup"]]))
})
