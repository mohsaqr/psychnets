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
#' **Data handling.** Non-numeric columns are dropped, and so is any column
#' that is constant once incomplete rows are removed. Rows with a missing value
#' in any numeric column are dropped (listwise deletion, as `goldbricker`
#' does), with a warning of class `psychnet_rows_dropped` giving the count; the
#' remaining number of rows is the `n` of every correlation-difference test.
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
#'   is stored with the result (attribute `proportion_matrix`). `print()`
#'   shows the flagged pairs, or a note that none were flagged, and returns `x`
#'   invisibly.
#' @param x A `psychnet_redundancy` object.
#' @param ... Unused.
#' @references Hallquist, M. N., Wright, A. G. C., & Molenaar, P. C. M. (2021).
#'   Problems with centrality measures in psychopathology networks.
#'   *Multivariate Behavioral Research*, 56(2), 199-223.
#' @examples
#' redundancy(SRL_Claude)
#' @export
redundancy <- function(data, p = 0.05, threshold = 0.25, cor_min = 0.5,
                       cor_method = c("auto", "pearson", "spearman", "kendall")) {
  cor_method <- .psn_match_arg(cor_method)
  .psn_check_number(p, "p", lower = 0, upper = 1, open = TRUE)
  .psn_check_number(threshold, "threshold", lower = 0, upper = 1)
  .psn_check_number(cor_min, "cor_min", lower = -1, upper = 1)
  mat <- .psn_item_data(data, NULL, "redundancy")$mat
  n <- nrow(mat)
  if (n < 4L)
    .psn_stop(sprintf(
      "redundancy() needs at least 4 complete rows (the test uses n - 3); got %d.", n),
      class = c("psychnet_degenerate_sample", "psychnet_bad_input"))
  labs <- colnames(mat)
  R <- if (cor_method == "auto") cor_auto(mat) else stats::cor(mat, method = cor_method)
  d <- ncol(R)
  if (d < 3L) .psn_stop("Need at least three items to assess redundancy.", call. = FALSE)

  # proportion[i,j]: share of items k whose correlation with i differs
  # significantly from its correlation with j.
  # (The comparisons of i and j with themselves give |r| = 1 -> NA and are
  # excluded, exactly as goldbricker excludes them.)
  prop <- matrix(NA_real_, d, d, dimnames = list(labs, labs))
  ij_all <- which(row(prop) != col(prop), arr.ind = TRUE)
  prop[ij_all] <- vapply(seq_len(nrow(ij_all)), function(r) {
    i <- ij_all[r, 1L]; j <- ij_all[r, 2L]
    pv <- .psn_hittner2003(R[, i], R[, j], R[i, j], n)
    pv <- pv[is.finite(pv)]
    if (length(pv)) mean(pv < p) else NA_real_
  }, numeric(1))

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

#' @describeIn redundancy The flagged pairs with the thresholds used.
#' @export
print.psychnet_redundancy <- function(x, ...) {
  cat(sprintf("# redundant pairs (proportion < %.2f, r > %.2f): %d found\n",
              attr(x, "threshold"), attr(x, "cor_min"), nrow(x)))
  if (nrow(x) == 0L) cat("  none\n") else print(`class<-`(x, "data.frame"))
  invisible(x)
}
