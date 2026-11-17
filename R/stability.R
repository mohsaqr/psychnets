# Centrality stability under case-dropping (Epskamp, Borsboom & Fried 2018),
# adapted from Nestimate. For each drop proportion the network is re-estimated
# on random case-dropped subsets and the subset centralities are correlated
# with the full-sample centralities. The CS-coefficient follows
# bootnet::corStability() (1.7.1): cor0() correlations, a subset is stable when
# its correlation is strictly above `threshold`, and CS is the largest drop
# proportion whose share of stable subsets is strictly above `certainty`.

#' Assess the stability of centralities when cases are dropped
#'
#' @description
#' `net_stability()` assesses how stable the centralities of a network are when
#' cases are dropped, using the case-dropping bootstrap (Epskamp, Borsboom &
#' Fried, 2018). The network is re-estimated on subsets of the data from which
#' an increasing proportion of cases has been removed, the centralities of each
#' subset are correlated with those of the full sample, and the
#' correlation-stability coefficient (CS-coefficient) summarizes the result:
#' the largest proportion of cases that can be dropped while the subset
#' centralities still correlate highly with the full-sample centralities in
#' nearly all subsets.
#'
#' Centrality indices are often used to rank the nodes of a psychometric
#' network, for example to identify the symptoms most strongly connected to the
#' others, and such a ranking is informative only if a similar ranking would be
#' found in another sample. Epskamp, Borsboom and Fried (2018) suggest that the
#' CS-coefficient should not be below 0.25, and preferably above 0.5, before
#' differences in centrality are interpreted.
#'
#' The function draws `iter` subsets without replacement at every proportion in
#' `drop_prop`, re-estimates the network on each, and correlates each subset's
#' centralities with the full-sample centralities, for every measure in
#' `measures`. The CS-coefficient of a measure is the largest proportion at
#' which the share of subsets whose correlation exceeds `threshold` itself
#' exceeds `certainty`, or 0 when no proportion qualifies.
#'
#' It returns the summary by drop proportion, the CS-coefficient of each
#' measure and the correlation of every subset, read with `as.data.frame()`.
#'
#' For the stability of the edge weights rather than the centralities, use
#' [net_casedrop_reliability()]; for intervals around each edge and
#' centrality, use [net_boot()].
#'
#' @details
#' Both comparisons in the CS rule are strict, as in `bootnet::corStability()`.
#' The correlation of a subset is 0 when either centrality vector is constant
#' across nodes, and undefined when fewer than two values are available. An
#' undefined correlation, or a subset on which the network could not be
#' estimated, counts as not exceeding `threshold`, so such subsets never raise
#' the CS-coefficient (`bootnet` instead leaves undefined correlations out). A
#' measure whose full-sample centrality is constant has a CS-coefficient of 0.
#'
#' Each subset is passed to the estimator as the full data would be (see
#' [net_boot()]). The subsets are drawn before any network is fitted, so the
#' result is reproducible with [set.seed()].
#'
#' @param data A data frame or matrix with one row per observation, or a
#'   [psychnet_group] estimated from a data set (assessed per group).
#' @param method The estimator, as in [psychnet()]. Default `"glasso"`.
#' @param measures Centrality measures to assess. Default
#'   `c("strength", "expected_influence")`; any measure accepted by
#'   [net_centralities()] can be added.
#' @param centrality_fn Optional function supplying measures that are not
#'   built in (see [net_centralities()]).
#' @param drop_prop Proportions of cases to drop, each in (0, 1). Default
#'   `seq(0.1, 0.9, by = 0.1)`. Each must remove at least one case and keep at
#'   least two.
#' @param iter Number of subsets per proportion. Default `100`.
#' @param threshold Correlation a subset must exceed to count as stable, a
#'   number in (0, 1). Default `0.7`.
#' @param certainty Share of stable subsets that must be exceeded, a number in
#'   (0, 1). Default `0.95`.
#' @param cor_method Correlation between subset and full-sample centralities:
#'   `"pearson"` (default, as in `bootnet`) or `"spearman"` (stability of the
#'   ranking).
#' @param labels Optional node labels, one per node.
#' @param estimator_args Named list of arguments for the estimator. Use it for
#'   an estimator argument that shares a name with an argument of
#'   `net_stability()`, e.g. `estimator_args = list(threshold = TRUE)`.
#' @param ... Further arguments for the estimator (`net_stability()`); ignored
#'   by the methods.
#' @param x A `psychnet_stability` object.
#' @param what For `as.data.frame()`: `"table"` (default), `"cs"` or
#'   `"correlations"`.
#' @param row.names,optional Ignored.
#' @return `net_stability()` returns an object of class `psychnet_stability`.
#'   Its tables are read with `as.data.frame()`:
#'
#'   * `as.data.frame(x)` (`what = "table"`): one row per measure and drop
#'     proportion, with columns `measure`, `drop_prop` (nominal),
#'     `actual_drop` (the proportion dropped after rounding to whole cases),
#'     `mean_cor` and `sd_cor` (over the defined correlations), `prop_above`
#'     (share of subsets with a correlation above `threshold`) and `n_valid`
#'     (subsets with a defined correlation).
#'   * `what = "cs"`: one row per measure, with columns `measure` and `cs`.
#'   * `what = "correlations"`: one row per measure and subset, with columns
#'     `measure`, `drop_prop`, `subset` and `correlation` (`NA` when undefined
#'     or when the subset could not be estimated), from which the
#'     CS-coefficient can be recomputed.
#'
#'   `print()` shows the estimator, the number of subsets per proportion and
#'   the CS-coefficient of each measure; `print()` and `plot()` return `x`
#'   invisibly. For a [psychnet_group], `net_stability()` returns one result
#'   per group (see [psychnet_result_group]).
#' @section Conditions:
#'   `psychnet_bad_input` for an invalid argument, including drop proportions
#'   that remove no case or keep fewer than two, and for a group that was not
#'   estimated from a data set; `psychnet_estimation_failed` when the network
#'   cannot be estimated on the full data; `psychnet_resample_failed` (warning)
#'   when some subsets cannot be estimated, naming the first error, and as an
#'   error when none can; `psychnet_constant_centrality` (warning) when a
#'   full-sample centrality is constant across nodes;
#'   `psychnet_constant_variable` (warning) when a constant column is dropped.
#' @seealso [net_casedrop_reliability()] for the edge weights; [net_boot()]
#'   for bootstrap intervals; [net_centralities()] for the measures.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy: A tutorial paper. *Behavior
#'   Research Methods*, 50(1), 195-212. \doi{10.3758/s13428-017-0862-1}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' # method = "pcor" and a small iter keep this example fast; a real analysis
#' # uses the "glasso" default and iter = 100.
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

#' @describeIn net_stability The estimator, the subsets per proportion and the
#'   CS-coefficient of each measure.
#' @export
print.psychnet_stability <- function(x, ...) {
  cat(sprintf("<psychnet_stability> %s, %d subsets/proportion\n",
              x$method, x$iter))
  cat(sprintf("  CS-coefficient (%s cor > %.2f in > %.0f%% of subsets):\n",
              x$cor_method %||% "spearman", x$threshold, 100 * x$certainty))
  cat(sprintf("    %-20s %.2f\n", names(x$cs), x$cs), sep = "")
  invisible(x)
}

#' @describeIn net_stability The summary by drop proportion
#'   (`what = "table"`), the CS-coefficient of each measure (`what = "cs"`),
#'   or the correlation of every subset (`what = "correlations"`).
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
