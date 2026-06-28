# Centrality stability under case-dropping (Epskamp, Borsboom & Fried 2018),
# adapted from Nestimate. For each drop proportion the network is re-estimated
# on random case-dropped subsets and the subset centralities are correlated
# (Spearman) with the full-sample centralities. The CS-coefficient is the
# largest drop proportion at which that correlation stays >= `threshold` with
# probability >= `certainty`.

#' Centrality-stability coefficient (case-dropping subset bootstrap)
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param measures Centrality measures to assess. Default both
#'   `c("strength", "expected_influence")`.
#' @param drop_prop Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.
#' @param iter Subsets per proportion. Default 100.
#' @param threshold Minimum acceptable rank correlation. Default 0.7.
#' @param certainty Probability the correlation must exceed `threshold`.
#'   Default 0.95.
#' @param labels Optional node labels.
#' @param ... Passed to the estimator.
#' @return An object of class `psychnet_stability` with `$cs` (CS-coefficient
#'   per measure) and a tidy `$table` of mean correlations by drop proportion.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' cs <- net_stability(x, drop_prop = c(0.3, 0.5, 0.7), iter = 20)
#' cs$cs
#' @export
net_stability <- function(data, method = "glasso",
                                 measures = c("strength", "expected_influence"),
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 iter = 100L, threshold = 0.7, certainty = 0.95,
                                 labels = NULL, ...) {
  measures <- match.arg(measures, c("strength", "expected_influence"),
                        several.ok = TRUE)
  stopifnot(length(drop_prop) >= 1L, all(drop_prop > 0), all(drop_prop < 1),
            is.numeric(iter), length(iter) == 1L, is.finite(iter), iter >= 1,
            threshold > 0, threshold <= 1, certainty > 0, certainty <= 1)
  iter <- as.integer(iter)   # a fractional count corrupts the stored %d field
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  if (is.null(labels)) labels <- colnames(mat)

  full_cent <- net_centralities(psychnet(mat, method = method,
                                           labels = labels, ...))

  # corr_storage[[measure]]: iter x length(drop_prop) Spearman correlations.
  corr_storage <- lapply(measures, function(m)
    matrix(NA_real_, iter, length(drop_prop)))
  names(corr_storage) <- measures

  for (pj in seq_along(drop_prop)) {
    keep_n <- max(2L, round(n * (1 - drop_prop[pj])))
    for (it in seq_len(iter)) {
      idx <- sample.int(n, keep_n, replace = FALSE)
      fit <- tryCatch(
        psychnet(mat[idx, , drop = FALSE], method = method,
                         labels = labels, ...),
        error = function(e) NULL)
      if (is.null(fit)) next
      ct <- net_centralities(fit)
      for (m in measures) {
        corr_storage[[m]][it, pj] <- suppressWarnings(
          stats::cor(full_cent[[m]], ct[[m]], method = "spearman"))
      }
    }
  }

  cs <- vapply(measures, function(m) {
    prop_above <- colMeans(corr_storage[[m]] >= threshold, na.rm = TRUE)
    valid <- which(prop_above >= certainty)
    # CS = the largest drop PROPORTION that stays stable, robust to the order in
    # which drop_prop was supplied (not the largest index).
    if (length(valid) == 0L) 0 else max(drop_prop[valid])
  }, numeric(1))

  tab <- do.call(rbind, lapply(measures, function(m) {
    cm <- corr_storage[[m]]
    data.frame(measure = m, drop_prop = drop_prop,
               mean_cor = colMeans(cm, na.rm = TRUE),
               prop_above = colMeans(cm >= threshold, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))

  structure(list(cs = cs, table = tab, drop_prop = drop_prop,
                 threshold = threshold, certainty = certainty,
                 iter = iter, method = method),
            class = "psychnet_stability")
}

#' Print a centrality-stability result
#'
#' @param x A `psychnet_stability` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_stability <- function(x, ...) {
  cat(sprintf("<psychnet_stability> %s, %d subsets/proportion\n",
              x$method, x$iter))
  cat(sprintf("  CS-coefficient (cor >= %.2f with %.0f%% certainty):\n",
              x$threshold, 100 * x$certainty))
  for (m in names(x$cs)) cat(sprintf("    %-20s %.2f\n", m, x$cs[[m]]))
  invisible(x)
}
