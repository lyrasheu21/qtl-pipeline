# Mediation and interaction models for colocalised cciQTLs.
#
# For each anchor (subclass:gene) and lead SNP, and each interaction partner
# (the other side of every CCI sharing that anchor and SNP):
#   mediation:   geno -> anchor expression -> partner expression
#                (mediation::mediate, anchor as mediator)
#   interaction: partner ~ anchor * geno; reports the anchor:geno term and
#                its partial eta-squared
# Expression is the residualised per-subclass expression (inputs.expression_pattern).

for (f in list.files(snakemake@params[["rdir"]], pattern = "\\.R$", full.names = TRUE)) source(f)
start_logging(snakemake@log[[1]])
suppressPackageStartupMessages(library(mediation))

set.seed(snakemake@params[["seed"]])
sims <- snakemake@params[["sims"]]
pattern <- snakemake@params[["expression_pattern"]]

hits <- fread(snakemake@input[["hits"]])
geno_raw <- fread(snakemake@input[["genotypes"]], colClasses = "character")

empty_mediation <- data.table(anchor = character(), snp = character(), partner = character(), gwas = character(),
                              n = integer(), acme = numeric(), acme_p = numeric(), ade = numeric(), ade_p = numeric(),
                              total = numeric(), total_p = numeric(), prop_mediated = numeric(), prop_mediated_p = numeric())
empty_interaction <- data.table(anchor = character(), snp = character(), partner = character(), gwas = character(),
                                n = integer(), interaction_beta = numeric(), interaction_p = numeric(),
                                partial_eta_sq = numeric())

gt_to_dosage <- function(gt) {
  gt <- gsub("|", "/", gt, fixed = TRUE)
  out <- rep(NA_real_, length(gt))
  out[gt == "0/0"] <- 0
  out[gt %in% c("0/1", "1/0")] <- 1
  out[gt == "1/1"] <- 2
  out
}

# One test per (anchor subclass, anchor gene, SNP, partner subclass, partner gene).
tests <- hits[, .(gwas = paste(sort(unique(gwas)), collapse = ";")),
              by = .(subclass, gene, snp, partner_subclass, partner_gene)]
message(sprintf("Mediation tests: %d", nrow(tests)))

expr_cache <- new.env()
get_expr <- function(sc) {
  if (!exists(sc, envir = expr_cache, inherits = FALSE)) {
    path <- fill_pattern(pattern, subclass = sc)
    assign(sc, if (file.exists(path)) read_expression(path) else NULL, envir = expr_cache)
    if (!file.exists(path)) message("Expression file missing: ", path)
  }
  get(sc, envir = expr_cache)
}

med_rows <- list()
int_rows <- list()
for (i in seq_len(nrow(tests))) {
  t <- tests[i]
  label <- sprintf("%s:%s %s -> %s:%s", t$subclass, t$gene, t$snp, t$partner_subclass, t$partner_gene)
  g_row <- geno_raw[ID == t$snp]
  anchor_expr <- get_expr(t$subclass)
  partner_expr <- get_expr(t$partner_subclass)
  if (nrow(g_row) == 0 || is.null(anchor_expr) || is.null(partner_expr) ||
      !t$gene %in% rownames(anchor_expr) || !t$partner_gene %in% rownames(partner_expr)) {
    message("Skipped (missing genotype or expression): ", label)
    next
  }
  donor_cols <- setdiff(names(geno_raw), c("CHROM", "POS", "ID", "REF", "ALT"))
  donors <- Reduce(intersect, list(donor_cols, colnames(anchor_expr), colnames(partner_expr)))
  df <- data.table(geno = gt_to_dosage(unlist(g_row[1, ..donors])),
                   anchor = as.numeric(anchor_expr[t$gene, donors]),
                   partner = as.numeric(partner_expr[t$partner_gene, donors]))
  df <- df[complete.cases(df)]
  if (nrow(df) < 10 || uniqueN(df$geno) < 2) {
    message("Skipped (too few donors or monomorphic): ", label)
    next
  }

  key <- list(anchor = paste0(t$subclass, ":", t$gene), snp = t$snp,
              partner = paste0(t$partner_subclass, ":", t$partner_gene), gwas = t$gwas, n = nrow(df))

  med_fit <- lm(anchor ~ geno, data = df)
  out_fit <- lm(partner ~ anchor + geno, data = df)
  med <- tryCatch(mediate(med_fit, out_fit, treat = "geno", mediator = "anchor", robustSE = TRUE, sims = sims),
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

mediation_res <- if (length(med_rows)) rbindlist(med_rows) else empty_mediation
interaction_res <- if (length(int_rows)) rbindlist(int_rows) else empty_interaction
setorder(mediation_res, acme_p)
setorder(interaction_res, interaction_p)
message(sprintf("Mediation results: %d (ACME p < 0.05: %d)", nrow(mediation_res), sum(mediation_res$acme_p < 0.05)))
write_tsv(mediation_res, snakemake@output[["mediation"]])
write_tsv(interaction_res, snakemake@output[["interaction"]])
