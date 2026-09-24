# Statistical helpers.

#' Storey q-values in base R (smoothed pi0, as the qvalue package defaults).
#' Used when qvalue::qvalue() fails on unusual p-value distributions.
storey_qvalue <- function(p, lambda = seq(0.05, 0.95, 0.05)) {
  m <- length(p)
  pi0_vec <- vapply(lambda, function(l) mean(p >= l) / (1 - l), numeric(1))
  pi0 <- if (m >= 4) {
    sp <- stats::smooth.spline(lambda, pi0_vec, df = 3)
    stats::predict(sp, x = max(lambda))$y
  } else 1
  pi0 <- min(max(pi0, 0), 1)
  if (pi0 <= 0) pi0 <- 1
  o <- order(p, decreasing = TRUE)
  q <- pmin(1, cummin(pi0 * m * p[o] / seq(m, 1)))
  list(qvalues = q[order(o)], pi0 = pi0)
}

#' q-values computed separately within each group (e.g. anchor subclass).
#'
#' Tries qvalue::qvalue(); if that fails, uses the base-R Storey
#' implementation for groups of at least `min_storey` p-values and
#' Benjamini-Hochberg (pi0 = 1) for smaller ones, where pi0 cannot be
#' estimated reliably. Pass a constant `group` for one global q-value.
#' @param method "storey" (q-values, default) or "bh" (Benjamini-Hochberg
#'   adjusted p-values, i.e. pi0 = 1; more conservative).
qvalue_by_group <- function(p, group, min_storey = 20, method = c("storey", "bh")) {
  method <- match.arg(method)
  stopifnot(length(p) == length(group))
  q <- rep(NA_real_, length(p))
  for (g in unique(group)) {
    idx <- which(group == g & !is.na(p))
    if (!length(idx)) next
    if (method == "bh") {
      q[idx] <- p.adjust(p[idx], method = "BH")
      next
    }
    q[idx] <- tryCatch(
      qvalue::qvalue(p[idx])$qvalues,
      error = function(e) {
        if (length(idx) >= min_storey) {
          message("qvalue failed for group '", g, "' (n = ", length(idx), "): ", conditionMessage(e),
                  "; using base-R Storey q-values")
          return(storey_qvalue(p[idx])$qvalues)
        }
        message("Group '", g, "' has only ", length(idx), " p-values; using Benjamini-Hochberg (pi0 = 1)")
        p.adjust(p[idx], method = "BH")
      }
    )
  }
  q
}

#' Storey's pi0 for a vector of p-values (NA when it cannot be estimated).
pi0_estimate <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) < 4) return(NA_real_)
  tryCatch(qvalue::qvalue(p)$pi0, error = function(e) storey_qvalue(p)$pi0)
}

#' Proportion of true alternatives, pi1 = 1 - pi0 (limma's estimator).
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

#' Rank-based inverse normal transform of each row (as QTLtools --normal).
rank_normal_rows <- function(mat) {
  t(apply(mat, 1, function(x) {
    out <- rep(NA_real_, length(x))
    ok <- !is.na(x)
    out[ok] <- qnorm((rank(x[ok]) - 0.5) / sum(ok))
    out
  }))
}
