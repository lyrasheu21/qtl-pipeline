# Chord diagram for one anchor gene ("Subclass:Gene"): every tested partner
# links to the anchor; links are red where mediation is significant
# (ACME p < 0.05) and grey otherwise. Built from mediation_results.tsv.

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
suppressPackageStartupMessages(library(circlize))

target <- snakemake@params[["target"]]
ct <- read_cell_types(snakemake@input[["cell_types"]])
med <- fread(snakemake@input[["mediation"]])[anchor == target]
out <- snakemake@output[[1]]

if (nrow(med) == 0) {
  empty_plot(out, paste("No mediation results for", target))
  quit(save = "no")
}

# Best (lowest ACME p) result per partner
med <- med[order(acme_p)][!duplicated(partner)]
med[, partner_subclass := sub(":.*$", "", partner)]
med[, partner_subclass := factor(partner_subclass, levels = union(subclass_order(ct), unique(partner_subclass)))]
setorder(med, partner_subclass)

links <- data.frame(from = med$partner, to = target, value = 1)
sub_cols <- subclass_colors(ct)
grid_col <- c(setNames(ifelse(as.character(med$partner_subclass) %in% names(sub_cols),
                              sub_cols[as.character(med$partner_subclass)], "grey60"), med$partner),
              setNames(if (sub(":.*$", "", target) %in% names(sub_cols)) sub_cols[[sub(":.*$", "", target)]] else "grey60",
                       target))
link_col <- ifelse(med$acme_p < 0.05, "#D7191CCC", "#99999966")

pdf(out, width = 7, height = 7)
circos.clear()
circos.par(start.degree = 90)
chordDiagram(links, grid.col = grid_col, col = link_col, directional = 1,
             direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.3),
             order = c(med$partner, target))
circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[1] + 0.1,
              get.cell.meta.data("sector.index"), facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.6)
})
title(target)
legend("bottomleft", legend = c("Mediation ACME p < 0.05", "Not significant"),
       fill = c("#D7191C", "#999999"), border = NA, bty = "n", cex = 0.8)
circos.clear()
invisible(dev.off())
