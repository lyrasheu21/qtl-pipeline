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
