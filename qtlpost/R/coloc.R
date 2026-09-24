# Colocalisation helpers: build coloc.abf input datasets for the cciQTL and
# GWAS sides, and classify how sharply a colocalisation points at one SNP.

#' Standardise a GWAS table to columns snp, bp, beta, se, p, n using the
#' column mapping from config/gwas.tsv.
standardize_gwas <- function(g, spec) {
  cols <- c(snp = spec$snp_col, bp = spec$bp_col, beta = spec$beta_col,
            se = spec$se_col, p = spec$p_col)
  missing <- setdiff(cols, names(g))
  if (length(missing)) stop("GWAS '", spec$name, "' is missing columns: ", paste(missing, collapse = ", "))
  out <- data.table(snp = as.character(g[[cols[["snp"]]]]),
                    bp = as.numeric(g[[cols[["bp"]]]]),
                    beta = as.numeric(g[[cols[["beta"]]]]),
                    se = as.numeric(g[[cols[["se"]]]]),
                    p = as.numeric(g[[cols[["p"]]]]))
  if (!is.null(spec$n_col) && !is.na(spec$n_col) && nzchar(spec$n_col)) {
    out[, n := as.numeric(g[[spec$n_col]])]
  } else if (!is.null(spec$n) && !is.na(spec$n)) {
    out[, n := as.numeric(spec$n)]
  } else {
    stop("GWAS '", spec$name, "' needs either n_col or n")
  }
  out <- out[!is.na(snp) & is.finite(beta) & is.finite(se) & se > 0 & !is.na(p) & p > 0]
  out[!duplicated(snp)]
}

#' coloc dataset for the cciQTL side, from QTLtools nominal rows of one CCI.
#' Uses slope_se from --std-err; falls back to back-calculating the SE from the
#' p-value when slope_se is missing (e.g. output produced without --std-err).
qtl_coloc_dataset <- function(nom, n_samples) {
  d <- nom[!is.na(slope) & !is.na(nom_pval) & !duplicated(var_id)]
  se <- d$slope_se
  if (all(is.na(se))) {
    t_stat <- abs(qt(pmax(d$nom_pval, .Machine$double.xmin) / 2, df = n_samples - 1))
    se <- abs(d$slope / t_stat)
  }
  keep <- is.finite(se) & se > 0
  d <- d[keep]
  se <- se[keep]
  list(beta = setNames(d$slope, d$var_id),
       varbeta = setNames(se^2, d$var_id),
       snp = d$var_id,
       position = setNames(d$var_from, d$var_id),
       N = n_samples,
       type = "quant",
       sdY = 1)  # QTLtools --normal rank-normalises each phenotype
}

#' coloc dataset for the GWAS side, from a standardised GWAS table.
gwas_coloc_dataset <- function(g, spec) {
  ds <- list(beta = setNames(g$beta, g$snp),
             varbeta = setNames(g$se^2, g$snp),
             snp = g$snp,
             position = setNames(g$bp, g$snp),
             N = max(g$n, na.rm = TRUE),
             type = spec$type)
  if (spec$type == "cc") {
    if (is.null(spec$s) || is.na(spec$s)) stop("GWAS '", spec$name, "' is type cc but has no s (case proportion)")
    ds$s <- as.numeric(spec$s)
  } else {
    ds$sdY <- 1
  }
  ds
}

#' Classify how sharply a colocalisation points at one SNP, from the SNP.PP.H4
#' values kept for a hit (the SNPs above the snp_pp_h4 threshold).
classify_coloc_tier <- function(snp_pp_h4) {
  n_strong <- sum(snp_pp_h4 > 0.5)
  n_near <- sum(snp_pp_h4 > 0.3)
  if (n_strong == 1) return("Gold_Unique")
  if (n_near > 1) return("Gold_LD_Pair")
  if (length(snp_pp_h4) > 1) return("Silver_Multiple")
  if (length(snp_pp_h4) == 1) return("Bronze_Single")
  "Diffuse"
}

#' Colocalisation of every phenotype in a nominal QTL table with one GWAS.
#'
#' @param nominal QTLtools nominal rows (read_qtltools(..., "nominal")) for the
#'   phenotypes to test.
#' @param gwas standardised GWAS table (standardize_gwas), ideally pre-filtered
#'   to the SNPs in `nominal`.
#' @param spec GWAS spec: list(name, type = "cc"/"quant", s).
#' @param n_qtl number of donors in the QTL analysis.
#' @param settings list(p1, p2, p12, min_snps, pp_h4, snp_pp_h4).
#' @param label extra columns to add to every output row (e.g. list(set = "x")).
#' @return list(summary = one row per phenotype tested, snps = SNP-level rows
#'   for colocalised phenotypes, results = full coloc SNP results per hit)
run_coloc <- function(nominal, gwas, spec, n_qtl, settings, label = list()) {
  summaries <- list()
  snp_rows <- list()
  full <- list()
  for (pid in unique(nominal$phe_id)) {
    d <- nominal[phe_id == pid]
    shared <- intersect(d$var_id, gwas$snp)
    if (length(shared) < settings$min_snps) next
    d <- d[var_id %chin% shared]
    g <- gwas[snp %chin% shared]
    res <- tryCatch(
      suppressMessages(suppressWarnings(coloc::coloc.abf(
        qtl_coloc_dataset(d, n_qtl), gwas_coloc_dataset(g, spec),
        p1 = settings$p1, p2 = settings$p2, p12 = settings$p12))),
      error = function(e) { message("coloc failed for ", pid, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) next

    s <- as.list(res$summary)
    row <- data.table(gwas = spec$name, phenotype_id = pid, n_snps = s$nsnps,
                      PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
                      PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, ratio_H4_H3 = s$PP.H4.abf / s$PP.H3.abf,
                      tier = NA_character_, lead_snp = NA_character_, lead_snp_pp_h4 = NA_real_,
                      gwas_beta = NA_real_, gwas_p = NA_real_, qtl_beta = NA_real_, qtl_p = NA_real_,
                      qtl_is_best_hit = NA)
    if (s$PP.H4.abf >= settings$pp_h4) {
      r <- as.data.table(res$results)
      if (!"position" %in% names(r)) r[, position := d$var_from[match(snp, d$var_id)]]
      kept <- r[SNP.PP.H4 > settings$snp_pp_h4]
      lead <- r[which.max(SNP.PP.H4)]
      g_lead <- g[snp == lead$snp][1]
      d_lead <- d[var_id == lead$snp][1]
      row[, `:=`(tier = classify_coloc_tier(kept$SNP.PP.H4),
                 lead_snp = lead$snp, lead_snp_pp_h4 = lead$SNP.PP.H4,
                 gwas_beta = g_lead$beta, gwas_p = g_lead$p,
                 qtl_beta = d_lead$slope, qtl_p = d_lead$nom_pval,
                 qtl_is_best_hit = d_lead$best_hit == 1)]
      if (nrow(kept)) {
        snp_rows[[pid]] <- data.table(gwas = spec$name, phenotype_id = pid, snp = kept$snp,
                                      position = kept$position, SNP.PP.H4 = kept$SNP.PP.H4)
      }
      full[[pid]] <- r
    }
    summaries[[pid]] <- row
  }
  summary <- rbindlist(summaries)
  if (nrow(summary) == 0) {
    summary <- data.table(gwas = character(), phenotype_id = character(), n_snps = numeric(),
                          PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(),
                          PP.H4 = numeric(), ratio_H4_H3 = numeric(), tier = character(),
                          lead_snp = character(), lead_snp_pp_h4 = numeric(), gwas_beta = numeric(),
                          gwas_p = numeric(), qtl_beta = numeric(), qtl_p = numeric(), qtl_is_best_hit = logical())
  }
  snps <- rbindlist(snp_rows)
  if (nrow(snps) == 0) {
    snps <- data.table(gwas = character(), phenotype_id = character(), snp = character(),
                       position = numeric(), SNP.PP.H4 = numeric())
  }
  for (nm in names(label)) {
    summary[, (nm) := label[[nm]]]
    snps[, (nm) := label[[nm]]]
  }
  list(summary = summary, snps = snps, results = full)
}

#' Locus plot for one colocalised phenotype: QTL and GWAS -log10(p),
#' SNP.PP.H4 from coloc, and a gene track.
#' @param d nominal QTL rows for the phenotype; g standardised GWAS rows;
#'   r coloc SNP results (run_coloc()$results[[id]]); annot from read_annotation (or NULL).
plot_locus <- function(d, g, r, title, subtitle = NULL, annot = NULL, focus_gene = NULL) {
  lead <- r$snp[which.max(r$SNP.PP.H4)]
  window <- range(d$var_from)
  pts <- rbind(
    d[, .(position = var_from, snp = var_id, value = -log10(nom_pval), track = "QTL -log10(p)")],
    g[, .(position = bp, snp, value = -log10(p), track = "GWAS -log10(p)")],
    r[, .(position, snp, value = SNP.PP.H4, track = "SNP.PP.H4")]
  )
  pts[, `:=`(track = factor(track, levels = c("QTL -log10(p)", "GWAS -log10(p)", "SNP.PP.H4")),
             is_lead = snp == lead)]
  p_assoc <- ggplot(pts, aes(position, value)) +
    geom_point(aes(colour = is_lead, size = is_lead)) +
    facet_wrap(~track, ncol = 1, scales = "free_y", strip.position = "left") +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "red"), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 0.8, `TRUE` = 2.5), guide = "none") +
    coord_cartesian(xlim = window) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_pipeline() + theme(strip.placement = "outside")
  if (is.null(annot)) return(p_assoc)
  p_genes <- gene_track(annot, d$var_chr[1], window[1], window[2], highlight = focus_gene) +
    coord_cartesian(xlim = window) + labs(x = sprintf("chr%s position", d$var_chr[1]))
  patchwork::wrap_plots(p_assoc, p_genes, ncol = 1, heights = c(3, 0.6))
}
