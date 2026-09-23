# Load the pipeline's shared R helpers for every test file.
for (f in list.files(file.path(testthat::test_path(), "..", "..", "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}
toy_cell_types <- function() {
  data.table(subclass = c("Astro", "OPC", "EN_L2_3_IT", "Micro"),
             class = c("Astro", "OPC", "EN", "Immune"),
             superclass = c("Glia", "Glia", "EN", "Glia"),
             color = c("#C19A6B", "#FFF380", "#659EC7", "#F75D59"),
             class_color = c("#197EC0", "#FED439", "#1A9993", "#FD8CC1"))
}
