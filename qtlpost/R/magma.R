# MAGMA gene-set analysis helpers.

#' Read a MAGMA gene-set analysis result (.gsa.out, or .sets.out from older
#' versions). Returns columns set, ngenes, beta, beta_std, se, p.
read_magma_gsa <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) {
    return(data.table(set = character(), ngenes = numeric(), beta = numeric(),
                      beta_std = numeric(), se = numeric(), p = numeric()))
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  dt <- fread(text = paste(lines, collapse = "\n"), header = TRUE)
  name_col <- intersect(c("FULL_NAME", "VARIABLE", "SET"), names(dt))[1]
  if (is.na(name_col)) stop("Unrecognised MAGMA output: ", path)
  data.table(set = as.character(dt[[name_col]]), ngenes = dt$NGENES, beta = dt$BETA,
             beta_std = dt$BETA_STD, se = dt$SE, p = dt$P)
}

#' MAGMA --set-annot table (gene ID, set) from annotated QTL results.
#' Sets are "<set>__<level value>", e.g. "cci_sending__Astro".
#' @param qt annotated permutation results (several sets rbind-ed).
#' @param threshold "fdr" (fdr_pass) or "nominal" (adj_beta_pval <= nominal_p).
#' @param level column to group genes by: "group", "class" or "superclass".
magma_gene_sets <- function(qt, annot, threshold = c("fdr", "nominal"), level = "group",
                            id_column = "gene_id", nominal_p = 0.05) {
  threshold <- match.arg(threshold)
  hits <- if (threshold == "fdr") qt[fdr_pass == TRUE] else qt[adj_beta_pval <= nominal_p]
  if (!id_column %in% names(annot)) stop("Annotation has no column '", id_column, "'")
  hits <- hits[!is.na(gene)]
  hits[, gene_id := annot[[id_column]][match(gene, annot$gene_name)]]
  message(sprintf("%d of %d genes mapped to %s", uniqueN(hits[!is.na(gene_id), gene]), uniqueN(hits$gene), id_column))
  sets <- unique(hits[!is.na(gene_id), .(gene_id, set = paste0(set, "__", get(level)))])
  setorder(sets, set, gene_id)
  sets[]
}

#' Heatmap of MAGMA results (-log10 p) per GWAS x gene set, one panel per
#' phenotype set. "*" = p < 0.05; "#" = BH-adjusted p < 0.05 across GWAS
#' within each gene set.
plot_magma_heatmap <- function(res, group_levels = NULL, title = NULL) {
  res <- copy(res)
  res[, qtl_set := sub("__.*$", "", set)]
  res[, level_value := sub("^.*?__", "", set)]
  res[, p_adj := p.adjust(p, method = "BH"), by = set]
  res[, label := fifelse(p_adj < 0.05, "#", fifelse(p < 0.05, "*", ""))]
  res[, level_value := factor(level_value, levels = union(group_levels, sort(unique(level_value))))]
  res[, gwas := factor(gwas, levels = rev(unique(gwas)))]
  ggplot(res, aes(level_value, gwas, fill = -log10(p))) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3) +
    facet_wrap(~qtl_set, ncol = 1) +
    scale_fill_gradient(low = "white", high = "red", name = "-log10(p)") +
    labs(x = NULL, y = "GWAS", title = title) +
    theme_pipeline() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
