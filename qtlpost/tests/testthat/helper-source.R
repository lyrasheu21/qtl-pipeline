# Load the package source for every test file (works without installing).
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
for (f in list.files(file.path(testthat::test_path(), "..", "..", "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}
toy_group_table <- function() {
  data.table(group = c("Astro", "OPC", "EN_L2_3_IT", "Micro"),
             class = c("Astro", "OPC", "EN", "Immune"),
             superclass = c("Glia", "Glia", "EN", "Glia"),
             color = c("#C19A6B", "#FFF380", "#659EC7", "#F75D59"),
             class_color = c("#197EC0", "#FED439", "#1A9993", "#FD8CC1"))
}
