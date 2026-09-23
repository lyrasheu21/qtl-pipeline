# Statistical helpers.

#' q-values computed separately within each group (e.g. anchor subclass).
#'
#' This fixes the original `group_by("subclass")`, which grouped by a string
#' constant and so computed q-values across all CCIs at once. Groups too
#' small for qvalue's pi0 estimate fall back to pi0 = 1 (Benjamini-Hochberg).
qvalue_by_group <- function(p, group) {
  stopifnot(length(p) == length(group))
  q <- rep(NA_real_, length(p))
  for (g in unique(group)) {
    idx <- which(group == g & !is.na(p))
    if (!length(idx)) next
    q[idx] <- tryCatch(
      qvalue::qvalue(p[idx])$qvalues,
      error = function(e) {
        message("qvalue pi0 estimate failed for group '", g, "' (n = ", length(idx),
                "); using pi0 = 1 (Benjamini-Hochberg): ", conditionMessage(e))
        p.adjust(p[idx], method = "BH")  # q-values with pi0 = 1 are BH-adjusted p-values
      }
    )
  }
  q
}

#' Proportion of true alternatives, pi1 = 1 - pi0 (limma's estimator, as in the
#' original analysis).
pi1_estimate <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) < 2) return(NA_real_)
  1 - limma::propTrueNull(p)
}

#' Genomic inflation factor from p-values.
genomic_lambda <- function(p) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
}

#' Fisher's method for combining two independent p-values.
fisher_combine <- function(p1, p2) {
  pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
}

#' Expected and observed -log10 p for a QQ plot.
qq_points <- function(p) {
  p <- sort(p[!is.na(p) & p > 0])
  data.table(expected = -log10(ppoints(length(p))), observed = -log10(p))
}

#' Partial eta-squared of one term in a linear model (type I sums of squares).
partial_eta_sq <- function(fit, term) {
  a <- anova(fit)
  if (!term %in% rownames(a)) return(NA_real_)
  ss <- a[term, "Sum Sq"]
  ss / (ss + a["Residuals", "Sum Sq"])
}
