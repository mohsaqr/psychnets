# Centrality stability under case-dropping (Epskamp, Borsboom & Fried 2018),
# adapted from Nestimate. For each drop proportion the network is re-estimated
# on random case-dropped subsets and the subset centralities are correlated
# (Spearman) with the full-sample centralities. The CS-coefficient is the
# largest drop proportion at which that correlation stays >= `threshold` with
# probability >= `certainty`.

#' Centrality-stability coefficient (case-dropping subset bootstrap)
#'
#' @param data Data frame or matrix (rows = observations), resampled exactly as
#'   given (see [net_boot()]).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param measures Centrality measures to assess. Defaults to the two
#'   recommended for psychometric networks (`c("strength",
#'   "expected_influence")`); `"betweenness"`/`"closeness"` and custom measures
#'   (via `centrality_fn`) are also accepted. See [net_centralities()].
#' @param centrality_fn Optional function supplying any non-built-in `measures`
#'   (see [net_centralities()]).
#' @param drop_prop Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.
#' @param iter Subsets per proportion. Default 100.
#' @param threshold Minimum acceptable rank correlation. Default 0.7.
#' @param certainty Probability the correlation must exceed `threshold`.
#'   Default 0.95.
#' @param labels Optional node labels.
#' @param estimator_args Named list of estimator arguments. Use this for names
#'   consumed by the stability diagnostic itself, such as estimator `threshold`.
#' @param ... Passed to the estimator.
#' @return An object of class `psychnet_stability` with `$cs` (CS-coefficient
#'   per measure) and a tidy `$table` (columns `measure`, `drop_prop`,
#'   `mean_cor`, `sd_cor`, `prop_above`) of the case-dropping correlations by
#'   drop proportion. Visualise it with [plot.psychnet_stability()].
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' cs <- net_stability(x, drop_prop = c(0.3, 0.6), iter = 10)
#' cs$cs
#' @export
net_stability <- function(data, method = "glasso",
                                 measures = c("strength", "expected_influence"),
                                 centrality_fn = NULL,
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 iter = 100L, threshold = 0.7, certainty = 0.95,
                                 labels = NULL, estimator_args = list(), ...) {
  # Group object -> case-drop each level from its stored cross-sectional data.
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_stability, "net_stability",
      "psychnet_stability_group",
      list(measures = measures, centrality_fn = centrality_fn,
           drop_prop = drop_prop, iter = iter, threshold = threshold,
           certainty = certainty, estimator_args = estimator_args)))
  }
  stopifnot(length(drop_prop) >= 1L, all(drop_prop > 0), all(drop_prop < 1),
            is.numeric(iter), length(iter) == 1L, is.finite(iter), iter >= 1,
            threshold > 0, threshold <= 1, certainty > 0, certainty <= 1)
  iter <- as.integer(iter)   # a fractional count corrupts the stored %d field
  inp <- .resample_input(data, method,
                         .resample_dots(list(...), estimator_args))
  mat <- inp$data; dots <- inp$dots
  n <- nrow(mat)
  keep_sizes <- pmax(2L, round(n * (1 - drop_prop)))
  if (n < 3L || any(keep_sizes >= n))
    stop("Each `drop_prop` must remove at least one case and retain at least two cases.",
         call. = FALSE)
  cent_of <- function(fit) net_centralities(fit, measures = measures,
                                             centrality_fn = centrality_fn)
  fit_net <- function(m, align = NULL)
    .psn_refit_network(m, method = method, labels = align, dots = dots)

  full <- fit_net(mat, labels)
  if (is.null(full)) stop("Full-sample estimation failed.", call. = FALSE)
  if (is.null(labels)) labels <- full$nodes$label
  full_cent <- cent_of(full)

  # corr_storage[[measure]]: iter x length(drop_prop) Spearman correlations.
  corr_storage <- lapply(measures, function(m)
    matrix(NA_real_, iter, length(drop_prop)))
  names(corr_storage) <- measures
  fit_ok <- matrix(FALSE, iter, length(drop_prop))

  for (pj in seq_along(drop_prop)) {
    keep_n <- max(2L, round(n * (1 - drop_prop[pj])))
    for (it in seq_len(iter)) {
      idx <- sample.int(n, keep_n, replace = FALSE)
      fit <- fit_net(mat[idx, , drop = FALSE], labels)
      if (is.null(fit)) next
      fit_ok[it, pj] <- TRUE
      ct <- cent_of(fit)
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
               sd_cor = apply(cm, 2L, stats::sd, na.rm = TRUE),
               prop_above = colMeans(cm >= threshold, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))

  n_attempted <- iter * length(drop_prop)
  n_success <- sum(fit_ok)
  n_failed <- n_attempted - n_success
  if (n_failed > 0L)
    warning(sprintf("%d of %d case-drop fits failed.",
                    n_failed, n_attempted), call. = FALSE)

  structure(list(cs = cs, table = tab, drop_prop = drop_prop,
                 threshold = threshold, certainty = certainty,
                 iter = iter, method = method,
                 n_success = n_success, n_failed = n_failed),
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
