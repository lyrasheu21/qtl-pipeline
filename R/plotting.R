# Shared ggplot theme and plotting helpers (replaces the copies of good_plot()).

suppressPackageStartupMessages(library(ggplot2))

theme_pipeline <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold", colour = "black"),
          strip.background = element_rect(fill = "#eeeeee", colour = "#eeeeee"),
          panel.border = element_rect(linewidth = 0.8, colour = "black", fill = NA),
          legend.position = "right")
}

DIRECTION_COLORS <- c(sending = "#28827A", receiving = "#C84D4C")

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
  if (nrow(g) == 0) {
    return(ggplot() + theme_void() + labs(y = "Genes"))
  }
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
