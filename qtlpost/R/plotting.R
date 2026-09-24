# Plotting helpers. Every plot takes the standard pipeline tables.

theme_pipeline <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold", colour = "black"),
          strip.background = element_rect(fill = "#eeeeee", colour = "#eeeeee"),
          panel.border = element_rect(linewidth = 0.8, colour = "black", fill = NA),
          legend.position = "right")
}

#' Colours for phenotype sets (CCI directions keep their usual colours).
set_colors <- function(sets) {
  sets <- unique(as.character(sets))
  fixed <- c(sending = "#28827A", receiving = "#C84D4C")
  cols <- setNames(grDevices::hcl.colors(max(length(sets), 2), "Dark 3")[seq_along(sets)], sets)
  for (s in sets) {
    hit <- names(fixed)[endsWith(s, names(fixed))]
    if (length(hit)) cols[[s]] <- fixed[[hit[1]]]
  }
  cols
}

save_plot <- function(plot, path, width = 7, height = 5) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height, limitsize = FALSE)
}

#' Write a one-page PDF with a message, for figures that have no data to show.
empty_plot <- function(path, message) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = message) + theme_void()
  save_plot(p, path, width = 5, height = 3)
}

#' Gene track for a genomic window, drawn from the annotation table.
gene_track <- function(annot, chrom, from, to, highlight = NULL) {
  # local copies: `chrom`, `start` and `end` are also column names in `annot`
  win_chrom <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  win_from <- from
  win_to <- to
  g <- annot[chrom_norm == win_chrom & end >= win_from & start <= win_to]
  if (nrow(g) == 0) return(ggplot() + theme_void() + labs(y = "Genes"))
  g[, row := (seq_len(.N) - 1) %% 3]
  g[, is_target := gene_name %in% highlight]
  ggplot(g) +
    geom_segment(aes(x = pmax(start, win_from), xend = pmin(end, win_to), y = row, yend = row,
                     colour = is_target), linewidth = 2) +
    geom_text(aes(x = (pmax(start, win_from) + pmin(end, win_to)) / 2, y = row + 0.4,
                  label = gene_name, colour = is_target), size = 2.5) +
    scale_colour_manual(values = c(`FALSE` = "navy", `TRUE` = "red"), guide = "none") +
    scale_y_continuous(limits = c(-0.5, 3)) +
    labs(y = "Genes") +
    theme_pipeline() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

#' FDR hits vs phenotypes tested per group, one panel per set.
plot_fdr_counts <- function(counts, gt) {
  counts <- copy(counts)
  counts[, group := factor(group, levels = group_order(gt, group))]
  ggplot(counts, aes(group, n_fdr, colour = group, size = n_tested)) +
    geom_point() +
    facet_wrap(~set, ncol = 1) +
    scale_colour_manual(values = group_colors(gt, counts$group), guide = "none") +
    scale_size(range = c(0.5, 4), name = "Phenotypes tested") +
    labs(x = "Group", y = "FDR-significant QTLs") +
    theme_pipeline() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' QQ plots of permutation p-values per group (one plot per superclass).
plot_qq_pages <- function(qt, gt) {
  qq <- qt[, qq_points(adj_beta_pval), by = .(set, superclass, group)]
  lapply(intersect(superclass_order(gt, qq$superclass), unique(qq$superclass)), function(sc) {
    ggplot(qq[superclass == sc], aes(expected, observed, colour = set)) +
      geom_abline(linetype = "dashed") +
      geom_point(size = 0.8, alpha = 0.6) +
      facet_wrap(~group, scales = "free") +
      scale_colour_manual(values = set_colors(qq$set)) +
      labs(title = sc, x = "Expected -log10(p)", y = "Observed -log10(p)") +
      theme_pipeline()
  })
}

#' Histograms of empirical and beta-approximated permutation p-values per set.
plot_pvalue_hist <- function(qt) {
  long <- rbind(qt[, .(set, type = "adj_emp_pval (permutation)", p = adj_emp_pval)],
                qt[, .(set, type = "adj_beta_pval (beta approximation)", p = adj_beta_pval)])
  ggplot(long[!is.na(p)], aes(p)) +
    geom_histogram(binwidth = 0.02, boundary = 0, fill = "steelblue", colour = "white") +
    facet_grid(set ~ type, scales = "free_y") +
    labs(x = "p-value", y = "Phenotypes") +
    theme_pipeline()
}

#' Hits vs number of hidden factors (K) from the sweep summary.
plot_ksweep <- function(sweep) {
  long <- melt(sweep, id.vars = c("set", "k"), measure.vars = c("n_fdr_01", "n_fdr_05", "n_fdr_10"),
               variable.name = "threshold", value.name = "n")
  long[, threshold := factor(threshold, levels = c("n_fdr_01", "n_fdr_05", "n_fdr_10"),
                             labels = c("FDR 1%", "FDR 5%", "FDR 10%"))]
  best <- sweep[, .SD[which.max(n_fdr_05)], by = set]
  ggplot(long, aes(k, n, colour = threshold)) +
    geom_line() + geom_point() +
    geom_vline(data = best, aes(xintercept = k), linetype = "dashed", colour = "grey50") +
    facet_wrap(~set, scales = "free_y") +
    labs(x = "Hidden factors (K)", y = "Significant QTLs", colour = NULL,
         title = "Hidden-factor sweep (dashed line: most hits at FDR 5%)") +
    theme_pipeline()
}

#' UpSet plot of genes with QTLs shared across classes (draws on the current device).
plot_upset_classes <- function(qt, gt, title) {
  membership <- unique(qt[!is.na(gene), .(gene, class)])
  if (uniqueN(membership$class) < 2) {
    graphics::plot.new(); graphics::text(0.5, 0.5, paste(title, "- fewer than two classes, no UpSet plot"))
    return(invisible())
  }
  wide <- as.data.frame.matrix(table(membership$gene, membership$class))
  wide[wide > 0] <- 1L
  sets <- intersect(class_order(gt, membership$class), names(wide))
  cols <- group_colors(gt, sets, level = "class")[sets]
  print(UpSetR::upset(wide, sets = rev(sets), keep.order = TRUE, order.by = "freq", nintersects = 20,
                      main.bar.color = "black", sets.bar.color = rev(unname(cols)),
                      mainbar.y.label = paste(title, "- genes in intersection"), sets.x.label = "Genes per class"))
}

#' Number of unique genes with a colocalised QTL per GWAS, stacked by group.
plot_coloc_summary <- function(master, gt, gwas_levels, pp_h4) {
  counts <- master[, .(n_genes = uniqueN(gene)), by = .(gwas, set, group)]
  counts[, `:=`(gwas = factor(gwas, levels = rev(gwas_levels)),
                group = factor(group, levels = group_order(gt, group)))]
  ggplot(counts, aes(gwas, n_genes, fill = group)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~set) +
    scale_fill_manual(values = group_colors(gt, counts$group)) +
    labs(x = "GWAS", y = sprintf("Genes with colocalised QTL (PP.H4 >= %g)", pp_h4), fill = "Group") +
    theme_pipeline()
}

#' Chord diagram for one anchor ("group:gene") from mediation results: every
#' tested partner links to the anchor; red links = ACME p < 0.05.
plot_chord <- function(med, target, gt, path) {
  med <- med[anchor == target]
  if (nrow(med) == 0) return(empty_plot(path, paste("No mediation results for", target)))
  med <- med[order(acme_p)][!duplicated(partner)]
  med[, partner_group := sub(":.*$", "", partner)]
  med[, partner_group := factor(partner_group, levels = group_order(gt, partner_group))]
  setorder(med, partner_group)
  cols <- group_colors(gt, c(as.character(med$partner_group), sub(":.*$", "", target)))
  grid_col <- c(setNames(cols[as.character(med$partner_group)], med$partner),
                setNames(cols[[sub(":.*$", "", target)]], target))
  link_col <- ifelse(med$acme_p < 0.05, "#D7191CCC", "#99999966")
  grDevices::pdf(path, width = 7, height = 7)
  on.exit(grDevices::dev.off())
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90)
  circlize::chordDiagram(data.frame(from = med$partner, to = target, value = 1), grid.col = grid_col,
                         col = link_col, directional = 1, direction.type = c("diffHeight", "arrows"),
                         link.arr.type = "big.arrow", annotationTrack = "grid",
                         preAllocateTracks = list(track.height = 0.3), order = c(med$partner, target))
  circlize::circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    circlize::circos.text(mean(circlize::get.cell.meta.data("xlim")), circlize::get.cell.meta.data("ylim")[1] + 0.1,
                          circlize::get.cell.meta.data("sector.index"), facing = "clockwise", niceFacing = TRUE,
                          adj = c(0, 0.5), cex = 0.6)
  })
  graphics::title(target)
  graphics::legend("bottomleft", legend = c("Mediation ACME p < 0.05", "Not significant"),
                   fill = c("#D7191C", "#999999"), border = NA, bty = "n", cex = 0.8)
  circlize::circos.clear()
  invisible(path)
}
