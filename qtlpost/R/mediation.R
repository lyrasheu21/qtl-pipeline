# Mediation and interaction models for colocalised QTLs:
#   genotype -> anchor expression (gene in group) -> partner expression.

#' Genotype strings (0/0, 0|1, ...) to allele dosage.
gt_to_dosage <- function(gt) {
  gt <- gsub("|", "/", gt, fixed = TRUE)
  out <- rep(NA_real_, length(gt))
  out[gt == "0/0"] <- 0
  out[gt %in% c("0/1", "1/0")] <- 1
  out[gt == "1/1"] <- 2
  out
}

#' Run mediation (mediation::mediate, anchor as mediator) and an interaction
#' model (partner ~ anchor * genotype) for each hit.
#'
#' @param hits data.table with snp, group, gene, partner_group, partner_gene, gwas.
#' @param genotypes table from `bcftools query` (CHROM POS ID REF ALT + one GT column per donor).
#' @param get_expression function(group) returning a genes x donors matrix (or NULL).
#' @return list(mediation = ..., interaction = ...)
run_mediation <- function(hits, genotypes, get_expression, sims = 1000, min_donors = 10) {
  empty_med <- data.table(anchor = character(), snp = character(), partner = character(), gwas = character(),
                          n = integer(), acme = numeric(), acme_p = numeric(), ade = numeric(), ade_p = numeric(),
                          total = numeric(), total_p = numeric(), prop_mediated = numeric(), prop_mediated_p = numeric())
  empty_int <- data.table(anchor = character(), snp = character(), partner = character(), gwas = character(),
                          n = integer(), interaction_beta = numeric(), interaction_p = numeric(),
                          partial_eta_sq = numeric())
  tests <- hits[!is.na(partner_gene) & !is.na(partner_group),
                .(gwas = paste(sort(unique(gwas)), collapse = ";")),
                by = .(group, gene, snp, partner_group, partner_gene)]
  message(sprintf("Mediation tests: %d", nrow(tests)))
  donor_cols <- setdiff(names(genotypes), c("CHROM", "POS", "ID", "REF", "ALT"))

  med_rows <- list()
  int_rows <- list()
  for (i in seq_len(nrow(tests))) {
    t <- tests[i]
    label <- sprintf("%s:%s %s -> %s:%s", t$group, t$gene, t$snp, t$partner_group, t$partner_gene)
    g_row <- genotypes[ID == t$snp]
    anchor_expr <- get_expression(t$group)
    partner_expr <- get_expression(t$partner_group)
    if (nrow(g_row) == 0 || is.null(anchor_expr) || is.null(partner_expr) ||
        !t$gene %in% rownames(anchor_expr) || !t$partner_gene %in% rownames(partner_expr)) {
      message("Skipped (missing genotype or expression): ", label)
      next
    }
    donors <- Reduce(intersect, list(donor_cols, colnames(anchor_expr), colnames(partner_expr)))
    df <- data.table(geno = gt_to_dosage(unlist(g_row[1, ..donors])),
                     anchor = as.numeric(anchor_expr[t$gene, donors]),
                     partner = as.numeric(partner_expr[t$partner_gene, donors]))
    df <- df[complete.cases(df)]
    if (nrow(df) < min_donors || uniqueN(df$geno) < 2) {
      message("Skipped (too few donors or monomorphic): ", label)
      next
    }
    key <- list(anchor = paste0(t$group, ":", t$gene), snp = t$snp,
                partner = paste0(t$partner_group, ":", t$partner_gene), gwas = t$gwas, n = nrow(df))

    med_fit <- lm(anchor ~ geno, data = df)
    out_fit <- lm(partner ~ anchor + geno, data = df)
    med <- tryCatch(mediation::mediate(med_fit, out_fit, treat = "geno", mediator = "anchor",
                                       robustSE = TRUE, sims = sims),
                    error = function(e) { message("mediate failed for ", label, ": ", conditionMessage(e)); NULL })
    if (!is.null(med)) {
      med_rows[[length(med_rows) + 1]] <- as.data.table(c(key, list(
        acme = med$d0, acme_p = med$d0.p, ade = med$z0, ade_p = med$z0.p,
        total = med$tau.coef, total_p = med$tau.p, prop_mediated = med$n0, prop_mediated_p = med$n0.p)))
    }
    int_fit <- lm(partner ~ anchor * geno, data = df)
    coefs <- summary(int_fit)$coefficients
    if ("anchor:geno" %in% rownames(coefs)) {
      int_rows[[length(int_rows) + 1]] <- as.data.table(c(key, list(
        interaction_beta = coefs["anchor:geno", "Estimate"],
        interaction_p = coefs["anchor:geno", "Pr(>|t|)"],
        partial_eta_sq = partial_eta_sq(int_fit, "anchor:geno"))))
    }
  }
  mediation_res <- if (length(med_rows)) rbindlist(med_rows) else empty_med
  interaction_res <- if (length(int_rows)) rbindlist(int_rows) else empty_int
  setorder(mediation_res, acme_p)
  setorder(interaction_res, interaction_p)
  list(mediation = mediation_res, interaction = interaction_res)
}
