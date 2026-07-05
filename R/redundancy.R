# Topological-overlap node redundancy (Hallquist, Wright & Molenaar 2021),
# ported from networktools::goldbricker. Two items are redundant when their
# correlations with every OTHER item are statistically indistinguishable AND the
# two items are themselves strongly correlated. The pairwise comparison is the
# Hittner-May-Silver (2003) test for two dependent, overlapping correlations,
# reimplemented in base R (validated against cocor).

# Two-sided p-value for H0: r_jk == r_jh (dependent, overlapping; shared index).
# r_kh is the correlation between the two non-shared variables. Vectorised over
# r_jk / r_jh; returns NA where a correlation is +/-1 (e.g. comparing a variable
# with itself).
.psn_hittner2003 <- function(r_jk, r_jh, r_kh, n) {
  # Invalid inputs (|r| == 1, e.g. a variable compared with itself) make the
  # transform/denominator non-finite; those become NA. Suppress the expected
  # NaN warnings, exactly as cocor/goldbricker do.
  suppressWarnings({
    z <- function(r) atanh(r)
    r_bt <- tanh((z(r_jk) + z(r_jh)) / 2)
    covm <- (r_kh * (1 - 2 * r_bt^2) -
             0.5 * r_bt^2 * (1 - 2 * r_bt^2 - r_kh^2)) / (1 - r_bt^2)^2
    stat <- sqrt(n - 3) * (z(r_jk) - z(r_jh)) / sqrt(2 - 2 * covm)
    p <- 2 * stats::pnorm(-abs(stat))
  })
  p[!is.finite(stat)] <- NA_real_
  p
}

#' Detect redundant node pairs ("goldbricker")
#'
#' Flags pairs of items that behave redundantly in a network: their correlations
#' with all other items are mostly statistically indistinguishable (a small
#' proportion of significantly different correlations) and the two items are
#' themselves strongly correlated. Ported from `networktools::goldbricker`.
#'
#' @param data A numeric data frame or matrix (rows = observations).
#' @param p Significance level for each pairwise correlation-difference test.
#'   Default 0.05.
#' @param threshold Maximum proportion of significantly different correlations
#'   for a pair to be flagged redundant. Default 0.25.
#' @param cor_min Minimum correlation between the two items themselves. Default
#'   0.5.
#' @param cor_method Correlation type: `"auto"` (default, [cor_auto()] -
#'   polychoric/polyserial as appropriate, matching goldbricker), `"pearson"`,
#'   `"spearman"`, or `"kendall"`.
#' @return A tidy `data.frame` (class `psychnet_redundancy`), one row per flagged
#'   pair (sorted most-redundant first), with columns `item1`, `item2`,
#'   `proportion` (share of significantly different correlations) and
#'   `correlation`. Zero rows when nothing is flagged. The full proportion matrix
#'   is in `attr(x, "proportion_matrix")`.
#' @references Hallquist, M. N., Wright, A. G. C., & Molenaar, P. C. M. (2021).
#'   Problems with centrality measures in psychopathology networks.
#'   *Multivariate Behavioral Research*, 56(2), 199-223.
#' @examples
#' redundancy(SRL_Claude)
#' @export
redundancy <- function(data, p = 0.05, threshold = 0.25, cor_min = 0.5,
                       cor_method = c("auto", "pearson", "spearman", "kendall")) {
  cor_method <- match.arg(cor_method)
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  labs <- colnames(mat)
  if (is.null(labs)) labs <- paste0("V", seq_len(ncol(mat)))
  R <- if (cor_method == "auto") cor_auto(mat) else stats::cor(mat, method = cor_method)
  d <- ncol(R)
  if (d < 3L) stop("Need at least three items to assess redundancy.", call. = FALSE)

  # proportion[i,j]: share of items k whose correlation with i differs
  # significantly from its correlation with j.
  prop <- matrix(NA_real_, d, d, dimnames = list(labs, labs))
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      if (i == j) next
      pv <- .psn_hittner2003(R[, i], R[, j], R[i, j], n)
      pv <- pv[is.finite(pv)]
      prop[i, j] <- if (length(pv)) mean(pv < p) else NA_real_
    }
  }

  ut <- upper.tri(prop)
  ij <- which(ut, arr.ind = TRUE)
  pr <- prop[ut]
  co <- R[ut]
  keep <- which(pr < threshold & co > cor_min)
  ord <- keep[order(pr[keep])]
  out <- data.frame(item1 = labs[ij[ord, 1L]], item2 = labs[ij[ord, 2L]],
                    proportion = pr[ord], correlation = co[ord],
                    stringsAsFactors = FALSE, row.names = NULL)
  attr(out, "proportion_matrix") <- prop
  attr(out, "p") <- p
  attr(out, "threshold") <- threshold
  attr(out, "cor_min") <- cor_min
  class(out) <- c("psychnet_redundancy", "data.frame")
  out
}

#' Print a redundancy result
#'
#' @param x A `psychnet_redundancy` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_redundancy <- function(x, ...) {
  cat(sprintf("# redundant pairs (proportion < %.2f, r > %.2f): %d found\n",
              attr(x, "threshold"), attr(x, "cor_min"), nrow(x)))
  if (nrow(x) == 0L) cat("  none\n") else print(`class<-`(x, "data.frame"))
  invisible(x)
}
