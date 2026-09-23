# UpSet plots of anchor genes shared across cell classes: all tested cciQTL
# genes (page 1) and FDR-significant cciQTL genes (page 2).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
suppressPackageStartupMessages(library(UpSetR))

ct <- read_cell_types(snakemake@input[["cell_types"]])
qt <- rbind(fread(snakemake@input[["sending"]]), fread(snakemake@input[["receiving"]]))
cols <- class_colors(ct)

draw_upset <- function(dt, title) {
  membership <- unique(dt[, .(gene, class)])
  if (uniqueN(membership$class) < 2) {
    plot.new(); text(0.5, 0.5, paste(title, "- fewer than two classes, no UpSet plot"))
    return(invisible())
  }
  wide <- as.data.frame.matrix(table(membership$gene, membership$class))
  wide[wide > 0] <- 1L
  sets <- intersect(names(cols), names(wide))
  sets <- c(sets, setdiff(names(wide), sets))
  bar_cols <- ifelse(sets %in% names(cols), cols[sets], "grey50")
  print(upset(wide, sets = rev(sets), keep.order = TRUE, order.by = "freq", nintersects = 20,
              main.bar.color = "black", sets.bar.color = rev(bar_cols),
              mainbar.y.label = paste(title, "- genes in intersection"), sets.x.label = "Genes per class"))
}

pdf(snakemake@output[[1]], width = 9, height = 6, onefile = TRUE)
draw_upset(qt, "All tested cciQTL genes")
draw_upset(qt[fdr_pass == TRUE], "FDR cciQTL genes")
invisible(dev.off())
