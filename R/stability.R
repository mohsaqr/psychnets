# Centrality stability under case-dropping (Epskamp, Borsboom & Fried 2018),
# adapted from Nestimate. For each drop proportion the network is re-estimated
# on random case-dropped subsets and the subset centralities are correlated
# with the full-sample centralities. The CS-coefficient follows
# bootnet::corStability() (1.7.1): cor0() correlations, a subset is stable when
# its correlation is strictly above `threshold`, and CS is the largest drop
# proportion whose share of stable subsets is strictly above `certainty`.

#' Centrality-stability coefficient (case-dropping subset bootstrap)
#'
#' Asks how far the sample can shrink before the ranking of node centralities
#' changes. For each drop proportion the network is re-estimated on `iter`
#' random subsets without replacement, and each subset's centralities are
#' correlated with the full-sample centralities. The CS-coefficient of a
#' measure is the largest drop proportion at which the share of subsets whose
#' correlation is strictly greater than `threshold` is itself strictly greater
#' than `certainty`, the rule of `bootnet::corStability()`; Epskamp, Borsboom
#' and Fried (2018) suggest that it should not fall below 0.25 and preferably
#' exceed 0.5 before centrality differences are interpreted.
#'
#' The correlation follows bootnet's `cor0()`: 0 when either centrality vector
#' is constant across nodes, undefined when fewer than two values are
#' available. Unlike bootnet, which drops undefined correlations from the share
#' (`na.rm = TRUE`), psychnets counts an undefined correlation (and a failed
#' subset fit) as failing, so the CS-coefficient is never raised by subsets
#' that could not be compared. A measure whose full-sample centrality is
#' constant (for example strength in a two-node network) has CS 0 and raises a
#' `psychnet_constant_centrality` warning.
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
#' @param threshold Correlation a subset must exceed (strictly) to count as
#'   stable. A single number in (0, 1). Default 0.7.
#' @param certainty Share of stable subsets that must be exceeded (strictly).
#'   A single number in (0, 1). Default 0.95.
#' @param cor_method Correlation between subset and full-sample centralities:
#'   `"pearson"` (default, as in bootnet) or `"spearman"` (rank stability).
#' @param labels Optional node labels.
#' @param estimator_args Named list of estimator arguments. `threshold`,
#'   `cor_method` and `iter` are arguments of this diagnostic, so an
#'   estimator's own `threshold` or `cor_method` must go here, e.g.
#'   `estimator_args = list(threshold = TRUE)`.
#' @param ... Passed to the estimator by `net_stability()`; ignored by the
#'   methods.
#' @param x A `psychnet_stability` object.
#' @param what For `as.data.frame()`: `"table"` (default) for the
#'   case-dropping table, `"cs"` for the CS-coefficients, or
#'   `"correlations"` for every subset's correlation.
#' @param row.names,optional Ignored (S3 consistency).
#' @return `net_stability()` returns an object of class `psychnet_stability`,
#'   read with `as.data.frame()`, `print()` and `plot()`.
#'
#'   `as.data.frame(x)` (`what = "table"`) returns one row per measure per drop
#'   proportion with columns `measure`, `drop_prop` (nominal), `actual_drop`
#'   (the proportion actually dropped after rounding to whole cases),
#'   `mean_cor`, `sd_cor` (over the defined correlations), `prop_above` (share
#'   of subsets with correlation strictly above `threshold`; undefined
#'   correlations and failed fits count as not above) and `n_valid` (subsets
#'   with a defined correlation). `as.data.frame(x, what = "cs")` returns one
#'   row per measure with columns `measure` and `cs`.
#'   `as.data.frame(x, what = "correlations")` returns one row per measure per
#'   subset with columns `measure`, `drop_prop`, `subset` and `correlation`
#'   (`NA` when undefined or when the subset fit failed), from which the
#'   CS-coefficient can be recomputed.
#'
#'   A failed subset fit raises a `psychnet_resample_failed` warning carrying
#'   the first error; a failed full-sample fit raises a
#'   `psychnet_estimation_failed` error with the estimator's message.
#'   `print()` shows the method, the subsets per drop proportion and each
#'   CS-coefficient; `plot()` draws the mean correlation against the
#'   proportion dropped, one line per measure with a +/- 1 SD band, the
#'   threshold, and the CS-coefficients in the legend. Both return `x`
#'   invisibly.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy: A tutorial paper. *Behavior
#'   Research Methods*, 50(1), 195-212. \doi{10.3758/s13428-017-0862-1}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' # method = "pcor" and small iter keep this example fast; the "glasso"
#' # default and iter = 100 are what a real analysis should use.
#' cs <- net_stability(x, method = "pcor", drop_prop = c(0.3, 0.6), iter = 10)
#' cs                              # the CS-coefficient of each measure
#' as.data.frame(cs, what = "cs")  # the same, as a data frame
#' as.data.frame(cs)               # correlations by drop proportion
#' plot(cs)                        # mean correlation by drop proportion
#' @export
net_stability <- function(data, method = "glasso",
                                 measures = c("strength", "expected_influence"),
                                 centrality_fn = NULL,
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 iter = 100L, threshold = 0.7, certainty = 0.95,
                                 cor_method = c("pearson", "spearman"),
                                 labels = NULL, estimator_args = list(), ...) {
  # Group object -> case-drop each level from its stored cross-sectional data.
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_stability, "net_stability",
      "psychnet_stability_group",
      list(measures = measures, centrality_fn = centrality_fn,
           drop_prop = drop_prop, iter = iter, threshold = threshold,
           certainty = certainty, cor_method = cor_method,
           estimator_args = .resample_dots(list(...), estimator_args))))
  }
  cor_method <- .psn_match_arg(cor_method)
  .psn_check_number(iter, "iter", lower = 1)
  .psn_check_number(threshold, "threshold", lower = 0, upper = 1, open = TRUE)
  .psn_check_number(certainty, "certainty", lower = 0, upper = 1, open = TRUE)
  if (!is.character(measures) || !length(measures) || anyNA(measures) ||
      anyDuplicated(measures))
    .psn_stop("`measures` must be a character vector of distinct measure names.")
  iter <- as.integer(iter)   # a fractional count corrupts the stored %d field
  inp <- .resample_input(data, method,
                         .resample_dots(list(...), estimator_args))
  mat <- inp$data; dots <- inp$dots
  n <- nrow(mat)
  keep_sizes <- .psn_keep_sizes(drop_prop, n)
  cent_of <- function(fit) net_centralities(fit, measures = measures,
                                             centrality_fn = centrality_fn)
  fit_net <- function(m, align = NULL)
    .psn_refit_network(m, method = method, labels = align, dots = dots)

  full <- .psn_fit_or_stop(fit_net(mat, labels), "Full-sample estimation")
  if (is.null(labels)) labels <- full$nodes$label
  full_cent <- cent_of(full)
  constant <- measures[vapply(measures, function(m)
    .psn_is_constant(full_cent[[m]]), logical(1))]
  if (length(constant)) {
    .psn_warn(sprintf(paste0(
      "Full-sample %s centrality is constant across nodes, so every ",
      "case-drop correlation with it is 0 and its CS-coefficient is 0."),
      paste(constant, collapse = ", ")), class = "psychnet_constant_centrality")
  }

  # Subsets drawn in the parent (proportion outer, iteration inner), then each
  # draw is fitted and correlated on its own; a failing draw is recorded.
  idx_list <- .psn_casedrop_indices(n, keep_sizes, iter)
  draws <- lapply(idx_list, function(idx) .psn_try_draw({
    ct <- cent_of(fit_net(mat[idx, , drop = FALSE], labels))
    lapply(stats::setNames(measures, measures), function(m)
      .psn_cor0(full_cent[[m]], ct[[m]], method = cor_method))
  }))
  fit_ok <- .psn_check_draws(draws, "case-drop fits", min_success = 1L)

  # corr_storage[[measure]]: iter x length(drop_prop), iteration fastest.
  corr_storage <- lapply(stats::setNames(measures, measures), function(m)
    matrix(vapply(draws, function(d)
      if (.psn_draw_ok(d)) d[[m]] else NA_real_, numeric(1)),
      iter, length(drop_prop)))
  # bootnet's rule, strict on both sides; an undefined correlation fails.
  cs <- vapply(measures, function(m)
    .psn_cs_coefficient(corr_storage[[m]], drop_prop, threshold, certainty),
    numeric(1))

  tab <- do.call(rbind, lapply(measures, function(m) {
    s <- apply(corr_storage[[m]], 2L, .psn_summarise_draws)
    data.frame(measure = m, drop_prop = drop_prop,
               actual_drop = 1 - keep_sizes / n,
               mean_cor = s["mean", ], sd_cor = s["sd", ],
               prop_above = colMeans(!is.na(corr_storage[[m]]) &
                                       corr_storage[[m]] > threshold),
               n_valid = as.integer(s["n_valid", ]),
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  n_success <- sum(fit_ok)
  structure(list(cs = cs, table = tab, drop_prop = drop_prop,
                 threshold = threshold, certainty = certainty,
                 cor_method = cor_method, correlations = corr_storage,
                 iter = iter, method = method,
                 n_success = n_success,
                 n_failed = length(fit_ok) - n_success),
            class = "psychnet_stability")
}

#' @describeIn net_stability Method, subsets and CS-coefficient per measure.
#' @export
print.psychnet_stability <- function(x, ...) {
  cat(sprintf("<psychnet_stability> %s, %d subsets/proportion\n",
              x$method, x$iter))
  cat(sprintf("  CS-coefficient (%s cor > %.2f in > %.0f%% of subsets):\n",
              x$cor_method %||% "spearman", x$threshold, 100 * x$certainty))
  cat(sprintf("    %-20s %.2f\n", names(x$cs), x$cs), sep = "")
  invisible(x)
}

#' @describeIn net_stability The case-dropping table (`what = "table"`),
#'   the CS-coefficient per measure (`what = "cs"`), or every subset's
#'   correlation (`what = "correlations"`).
#' @export
as.data.frame.psychnet_stability <- function(x, row.names = NULL,
                                             optional = FALSE, ...,
                                             what = c("table", "cs",
                                                      "correlations")) {
  what <- .psn_match_arg(what)
  switch(what,
    table = x$table,
    cs = data.frame(measure = names(x$cs), cs = unname(x$cs),
                    stringsAsFactors = FALSE, row.names = NULL),
    correlations = {
      if (is.null(x$correlations))
        .psn_stop("This psychnet_stability object stores no correlations ",
                  "(it was made by an older psychnets); re-run net_stability().")
      do.call(rbind, lapply(names(x$correlations), function(m) {
        cm <- x$correlations[[m]]
        data.frame(measure = m,
                   drop_prop = rep(x$drop_prop, each = nrow(cm)),
                   subset = rep(seq_len(nrow(cm)), times = ncol(cm)),
                   correlation = as.vector(cm),
                   stringsAsFactors = FALSE, row.names = NULL)
      }))
    })
}
