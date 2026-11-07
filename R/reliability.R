# Two robustness diagnostics adapted from Nestimate, in the psychnets idiom
# (raw `data` in, re-estimate via psychnet(method=), base R only):
#   * net_casedrop_reliability() - edge-weight case-dropping stability (the
#     edge-vector complement to net_stability()'s centrality CS-coefficient).
#   * net_split_reliability()  - split-half reliability of the edge structure.
# Both are estimator-agnostic: they route every refit through psychnet(method=).

# Off-diagonal edge vector of a weight matrix (upper triangle when undirected,
# every off-diagonal cell when directed). Internal.
.psn_edge_vector <- function(weights, directed) {
  mask <- if (isTRUE(directed)) row(weights) != col(weights) else
          upper.tri(weights)
  as.vector(weights[mask])
}

# ---- shared resampling helpers ----------------------------------------------
#
# Used by net_boot(), net_stability(), net_casedrop_reliability() and
# net_split_reliability().

# TRUE when the finite values of `v` are all equal, to a tolerance relative to
# the vector's own magnitude (so float noise counts as constant). Scale-free on
# purpose: a correlation is unit-invariant, so constancy must be too -- an
# absolute floor made a small but varying centrality (1e-15 * 1:3) "constant".
#' @noRd
.psn_is_constant <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(FALSE)
  (max(v) - min(v)) <= 100 * .Machine$double.eps * max(abs(v))
}

# Correlation under bootnet's `cor0()` rule (bootnet 1.7.1, used by
# corStability()): NA when either side has fewer than two observed values or a
# non-finite one; 0 when either side is constant; otherwise stats::cor().
# `constant` = "na" (split-half reliability) is undefined only when BOTH sides
# are constant -- two empty halves agree perfectly, so 0 would be wrong -- and
# 0 when exactly one is (an empty half against a non-empty one is a real
# disagreement; dropping it would inflate the mean correlation).
#' @noRd
.psn_cor0 <- function(x, y, method = "pearson", constant = c("zero", "na")) {
  constant <- .psn_match_arg(constant)
  xo <- x[!is.na(x)]
  yo <- y[!is.na(y)]
  if (length(xo) < 2L || length(yo) < 2L ||
      !all(is.finite(xo)) || !all(is.finite(yo))) return(NA_real_)
  cx <- .psn_is_constant(x); cy <- .psn_is_constant(y)
  if (cx && cy) return(if (identical(constant, "zero")) 0 else NA_real_)
  if (cx || cy) return(0)
  stats::cor(x, y, method = method)
}

# The four bootnet-style edge-vector similarity metrics. Internal.
.psn_edge_metrics <- function(a, b, cor_method, constant = "zero") {
  d <- abs(a - b)
  c(mean_abs_dev   = mean(d),
    median_abs_dev = stats::median(d),
    correlation    = .psn_cor0(a, b, cor_method, constant = constant),
    max_abs_dev    = max(d))
}

# Mean / sd / percentile summary of one statistic across draws, over its
# defined (finite) values only; NA -- never NaN -- when too few are defined.
#' @noRd
.psn_summarise_draws <- function(v) {
  v <- v[is.finite(v)]
  k <- length(v)
  c(mean = if (k) mean(v) else NA_real_,
    sd = if (k >= 2L) stats::sd(v) else NA_real_,
    lower = if (k) stats::quantile(v, 0.025, names = FALSE) else NA_real_,
    upper = if (k) stats::quantile(v, 0.975, names = FALSE) else NA_real_,
    n_valid = k)
}

# The observed / full-sample fit of a resampling verb. An estimator error is
# re-raised with its own message (and its psychnet classes, if any), so
# `net_boot(x, foo = 1)` says "unused argument (foo = 1)" instead of a bare
# "estimation failed".
#' @noRd
.psn_fit_or_stop <- function(expr, what) {
  tryCatch(expr, error = function(e) {
    inherited <- setdiff(class(e), c("simpleError", "error", "condition"))
    .psn_stop(sprintf("%s failed: %s", what, conditionMessage(e)),
              class = unique(c("psychnet_estimation_failed", inherited)))
  })
}

# One resample draw. An error inside the draw (estimator, centrality,
# predictability) is caught HERE, in the worker, and recorded with its message,
# so a failing draw costs only itself: with forked workers an uncaught error
# would otherwise poison every draw prescheduled to the same worker, and the
# result would depend on the number of cores. The column-drop notices of the
# estimator input check are muted inside a draw (.psn_with_quiet_columns()):
# the observed / full-sample fit raised them once already, and a column that
# is constant in one resample is re-expanded as an isolated node.
#' @noRd
.psn_try_draw <- function(expr) {
  tryCatch(.psn_with_quiet_columns(expr), error = function(e)
    structure(list(message = conditionMessage(e)),
              class = "psychnet_draw_failure"))
}

# Did a draw succeed? A failure marker, a crashed fork ("try-error") or a NULL
# returned by mclapply for a killed worker all count as failed.
#' @noRd
.psn_draw_ok <- function(d) {
  is.list(d) && !inherits(d, "psychnet_draw_failure")
}

# The error message a failed draw carries.
#' @noRd
.psn_draw_message <- function(d) {
  if (inherits(d, "psychnet_draw_failure")) return(d$message)
  if (inherits(d, "try-error"))
    return(conditionMessage(attr(d, "condition")) %||% as.character(d))
  "the worker process returned no result"
}

# Count the failed draws of a resampling verb, stop when fewer than
# `min_success` succeeded, and otherwise warn once (class
# `psychnet_resample_failed`) with the first error message. Returns the
# logical success vector invisibly.
#' @noRd
.psn_check_draws <- function(draws, what, min_success = 2L) {
  ok <- vapply(draws, .psn_draw_ok, logical(1))
  n_total <- length(draws)
  n_failed <- sum(!ok)
  first <- if (n_failed) .psn_draw_message(draws[[which(!ok)[1L]]]) else NULL
  if (sum(ok) < min_success) {
    .psn_stop(sprintf("Fewer than %d of %d %s could be estimated%s.",
                      min_success, n_total, what,
                      if (is.null(first)) "" else
                        sprintf(" (first error: %s)", first)),
              class = c("psychnet_resample_failed", "psychnet_estimation_failed"))
  }
  if (n_failed > 0L) {
    .psn_warn(sprintf(paste0("%d of %d %s failed; summaries exclude those ",
                             "draws (first error: %s)."),
                      n_failed, n_total, what, first),
              class = "psychnet_resample_failed")
  }
  invisible(ok)
}

# The CS-coefficient rule of bootnet::corStability() (1.7.1) on a matrix of
# case-drop correlations (rows = subsets, columns = drop proportions): a subset
# is stable when its correlation is strictly above `threshold`, and CS is the
# largest drop proportion whose share of stable subsets is strictly above
# `certainty`, 0 when none is. Unlike bootnet (which uses na.rm = TRUE), an
# undefined correlation counts as not stable.
#' @noRd
.psn_cs_coefficient <- function(corr, drop_prop, threshold, certainty) {
  stable <- !is.na(corr) & corr > threshold
  valid <- which(colMeans(stable) > certainty)
  # The largest drop PROPORTION, robust to the order drop_prop was given in.
  if (length(valid)) max(drop_prop[valid]) else 0
}

# Validate the drop proportions of a case-dropping verb and return the
# retained sample size per proportion.
#' @noRd
.psn_keep_sizes <- function(drop_prop, n) {
  if (!is.numeric(drop_prop) || !length(drop_prop) || anyNA(drop_prop) ||
      any(drop_prop <= 0) || any(drop_prop >= 1)) {
    .psn_stop("`drop_prop` must be a numeric vector of proportions in (0, 1).")
  }
  keep_sizes <- pmax(2L, as.integer(round(n * (1 - drop_prop))))
  if (n < 3L || any(keep_sizes >= n))
    .psn_stop(paste0("Each `drop_prop` must remove at least one case and ",
                     "retain at least two cases."))
  keep_sizes
}

# Draw the case-dropping subsets in the parent, in the order the former nested
# loop drew them (proportion outer, iteration inner), so a seeded call
# reproduces the same subsets. Returns one row index vector per
# (iteration, proportion), iteration varying fastest.
#' @noRd
.psn_casedrop_indices <- function(n, keep_sizes, iter) {
  pj <- rep(seq_along(keep_sizes), each = iter)
  lapply(pj, function(j) sample.int(n, keep_sizes[j], replace = FALSE))
}

# Re-estimate a complete psychnet object, aligned to the full label set. A node
# that becomes constant in a resample is represented as an isolated node rather
# than changing the dimension or causing the draw to be discarded. Estimator
# errors propagate: the caller decides whether that is fatal (observed fit, via
# .psn_fit_or_stop()) or a recorded draw failure (via .psn_try_draw()).
#' @noRd
.psn_refit_network <- function(mat, method, labels = NULL, dots = list()) {
  source_by_label <- NULL
  if (!is.null(labels)) {
    # Labels refer to modeled nodes, which can be fewer than the raw columns
    # (e.g. a character id is ignored by a GGM). Rename those columns before
    # estimation so a node dropped in one draw can still be aligned by name.
    cols <- if (!is.data.frame(mat) || identical(.resolve_method(method), "mgm"))
      seq_len(ncol(mat)) else which(vapply(mat, is.numeric, logical(1)))
    if (length(labels) != length(cols)) {
      if (all(labels %in% colnames(mat)[cols])) {
        cols <- match(labels, colnames(mat))
      } else {
        # Listwise deletion or an explicit vars selection can reduce the node
        # set further. Resolve that set using the same estimator rather than
        # mistaking custom labels for absent variables and filling zeros.
        probe <- do.call(psychnet, c(list(mat, method = method), dots))
        cols <- match(probe$nodes$name, colnames(mat))
      }
    }
    if (length(labels) != length(cols) || anyNA(cols) || anyNA(labels) ||
        anyDuplicated(labels))
      .psn_stop("`labels` must uniquely name every modeled node.", call. = FALSE)
    source_by_label <- stats::setNames(colnames(mat)[cols] %||% labels, labels)
    if (!is.null(dots$vars) && is.character(dots$vars)) {
      selected <- match(dots$vars, unname(source_by_label))
      dots$vars[!is.na(selected)] <- labels[selected[!is.na(selected)]]
    }
    colnames(mat)[cols] <- labels
  }
  fit <- do.call(psychnet, c(list(mat, method = method), dots))
  if (is.null(labels)) return(fit)
  fit$nodes$name <- unname(source_by_label[fit$nodes$label])
  if (!is.null(fit$nodewise))
    fit$nodewise$source_names <- fit$nodes$name
  W <- fit$weights
  if (!identical(rownames(W), labels)) {
    old <- rownames(W)
    full <- matrix(0, length(labels), length(labels),
                   dimnames = list(labels, labels))
    keep <- intersect(old, labels)
    full[keep, keep] <- W[keep, keep]
    fit$weights <- full
    fit$nodes <- data.frame(id = seq_along(labels), label = labels,
                            name = unname(source_by_label[labels]),
                            x = NA_real_, y = NA_real_, stringsAsFactors = FALSE)
    fit$edges <- .edges_index_df(full, fit$directed)
    # Precision-based predictability is undefined for a node absent from the
    # fitted draw. Preserve the available block and mark the rest missing.
    if (!is.null(fit$precision)) {
      K <- matrix(NA_real_, length(labels), length(labels),
                  dimnames = list(labels, labels))
      K[keep, keep] <- fit$precision[keep, keep, drop = FALSE]
      fit$precision <- K
    }
  }
  fit
}

# Lightweight compatibility wrapper used by the reliability diagnostics.
.psn_refit <- function(mat, method, labels, dots = list()) {
  fit <- .psn_refit_network(mat, method, labels, dots)
  W <- fit$weights
  list(weights = W, directed = isTRUE(fit$directed))
}

# ---- edge-weight case-dropping stability -----------------------------------

#' Edge-weight stability coefficient (case-dropping subset bootstrap)
#'
#' The edge-vector complement to [net_stability()]. For each drop proportion the
#' network is re-estimated on random case-dropped subsets and the subset
#' edge-weight vector is compared with the full-sample one. The edge-weight
#' CS-coefficient follows the rule of `bootnet::corStability()` (Epskamp,
#' Borsboom & Fried 2018): a subset is stable when its edge-vector correlation
#' is strictly greater than `threshold`, and the CS-coefficient is the largest
#' drop proportion at which the share of stable subsets is strictly greater
#' than `certainty` (0 when no proportion qualifies).
#'
#' The correlation follows bootnet's `cor0()`: 0 when either edge vector is
#' constant (for example an empty subset network), undefined when fewer than
#' two edges are available. Unlike bootnet, which drops undefined correlations
#' from the share (`na.rm = TRUE`), psychnets counts an undefined correlation
#' (and a failed subset fit) as failing, so the CS-coefficient is never raised
#' by subsets that could not be compared.
#'
#' @param data Data frame or matrix (rows = observations), resampled exactly as
#'   given (see [net_boot()]), or a `psychnet_group` (case-dropped per level).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param drop_prop Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.
#' @param iter Subsets per proportion. Default 100.
#' @param threshold Correlation a subset must exceed (strictly) to count as
#'   stable. A single number in (0, 1). Default 0.7.
#' @param certainty Share of stable subsets that must be exceeded (strictly).
#'   A single number in (0, 1). Default 0.95.
#' @param cor_method Correlation method for the edge-vector comparison:
#'   `"pearson"` (default, as in bootnet), `"spearman"`, or `"kendall"`.
#' @param labels Optional node labels.
#' @param estimator_args Named list of estimator arguments. `cor_method`,
#'   `iter`, `threshold` and `drop_prop` are arguments of this diagnostic (the
#'   edge comparison), so an estimator's own `cor_method` or `threshold` must
#'   go here, e.g. `estimator_args = list(cor_method = "spearman")`; values in
#'   `...` with those names never reach the estimator.
#' @param ... Passed to the estimator by `net_casedrop_reliability()`;
#'   ignored by the methods.
#' @param x A `psychnet_casedrop` object.
#' @return A tidy `data.frame` (class `psychnet_casedrop`), one row per metric
#'   per drop proportion, with columns `metric`, `drop_prop` (nominal),
#'   `actual_drop` (the proportion actually dropped after rounding to whole
#'   cases), `mean`, `sd` and `n_valid` (subsets in which the metric was
#'   defined; `mean`/`sd` are `NA` when too few were). The edge-weight
#'   CS-coefficient is stored with it and shown when the result is printed.
#'   Undefined correlations raise a `psychnet_undefined_statistic` warning; a
#'   constant full-sample edge vector (e.g. an empty network) raises a
#'   `psychnet_constant_centrality` warning (its CS-coefficient is 0).
#'   A failed subset fit raises a `psychnet_resample_failed` warning carrying
#'   the first error; a failed full-sample fit raises a
#'   `psychnet_estimation_failed` error with the estimator's message. `plot()`
#'   draws the four similarity metrics (correlation, mean/median/max absolute
#'   deviation) against the proportion dropped, each with a +/- 1 SD band, the
#'   correlation panel carrying the threshold and the CS-coefficient. `print()`
#'   and `plot()` return `x` invisibly.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy. *Behavior Research Methods*,
#'   50(1), 195-212.
#' @examples
#' # `method`, `iter`, and `drop_prop` are chosen so the example runs quickly;
#' # the defaults (method = "glasso", iter = 100, drop_prop = seq(0.1, 0.9,
#' # 0.1)) are what a real reliability assessment should use.
#' rel <- net_casedrop_reliability(SRL_Claude, method = "pcor", iter = 5,
#'                                 drop_prop = c(0.25, 0.5))
#' rel
#' plot(rel)
#' @export
net_casedrop_reliability <- function(data, method = "glasso",
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 iter = 100L, threshold = 0.7, certainty = 0.95,
                                 cor_method = c("pearson", "spearman", "kendall"),
                                 labels = NULL, estimator_args = list(), ...) {
  # Group object -> case-drop each level from its stored cross-sectional data.
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_casedrop_reliability, "net_casedrop_reliability",
      "psychnet_casedrop_group",
      list(drop_prop = drop_prop, iter = iter, threshold = threshold,
           certainty = certainty, cor_method = cor_method,
           estimator_args = .resample_dots(list(...), estimator_args))))
  }
  cor_method <- .psn_match_arg(cor_method)
  .psn_check_number(iter, "iter", lower = 1)
  .psn_check_number(threshold, "threshold", lower = 0, upper = 1, open = TRUE)
  .psn_check_number(certainty, "certainty", lower = 0, upper = 1, open = TRUE)
  iter <- as.integer(iter)
  inp <- .resample_input(data, method,
                         .resample_dots(list(...), estimator_args))
  mat <- inp$data; dots <- inp$dots
  n <- nrow(mat)
  keep_sizes <- .psn_keep_sizes(drop_prop, n)

  full <- .psn_fit_or_stop(.psn_refit(mat, method, labels, dots),
                           "Full-sample estimation")
  # Labels come from the observed fit: the estimator decides which columns are
  # nodes, and every draw is re-expanded to that node set.
  if (is.null(labels)) labels <- rownames(full$weights)
  orig <- .psn_edge_vector(full$weights, full$directed)
  if (length(orig) >= 2L && .psn_is_constant(orig)) {
    .psn_warn(paste0("The full-sample network has constant edge weights (e.g. ",
                     "no edges): every edge-vector correlation is 0, so the ",
                     "edge-weight CS-coefficient is 0."),
              class = "psychnet_constant_centrality")
  }

  metric_names <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")
  idx_list <- .psn_casedrop_indices(n, keep_sizes, iter)
  draws <- lapply(idx_list, function(idx) .psn_try_draw({
    sub <- .psn_refit(mat[idx, , drop = FALSE], method, labels, dots)
    as.list(.psn_edge_metrics(orig, .psn_edge_vector(sub$weights, sub$directed),
                              cor_method))
  }))
  fit_ok <- .psn_check_draws(draws, "case-drop fits", min_success = 1L)
  # store[[metric]]: iter x length(drop_prop), iteration varying fastest.
  store <- stats::setNames(lapply(metric_names, function(nm)
    matrix(vapply(draws, function(d)
      if (.psn_draw_ok(d)) d[[nm]] else NA_real_, numeric(1)),
      iter, length(drop_prop))), metric_names)

  tab <- do.call(rbind, lapply(metric_names, function(nm) {
    s <- apply(store[[nm]], 2L, .psn_summarise_draws)
    data.frame(metric = nm, drop_prop = drop_prop,
               actual_drop = 1 - keep_sizes / n,
               mean = s["mean", ], sd = s["sd", ],
               n_valid = as.integer(s["n_valid", ]),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  n_undefined <- sum(fit_ok) - sum(is.finite(store$correlation))
  if (n_undefined > 0L) {
    .psn_warn(sprintf(paste0(
      "The edge-vector correlation is undefined in %d of %d estimated subsets ",
      "(fewer than two edges, or a non-finite weight); those subsets count as ",
      "failing the CS threshold."), n_undefined, sum(fit_ok)),
      class = "psychnet_undefined_statistic")
  }
  # bootnet's rule, strict on both sides; an undefined correlation fails.
  prop_above <- colMeans(!is.na(store$correlation) &
                           store$correlation > threshold)
  cs <- .psn_cs_coefficient(store$correlation, drop_prop, threshold, certainty)
  n_success <- sum(fit_ok); n_failed <- length(fit_ok) - n_success

  # The result IS the tidy table (one row per metric per drop proportion);
  # the CS-coefficient and raw draws ride along as attributes for plot()/print().
  attr(tab, "cs") <- cs
  attr(tab, "correlations") <- store$correlation
  attr(tab, "prop_above") <- prop_above
  attr(tab, "threshold") <- threshold
  attr(tab, "certainty") <- certainty
  attr(tab, "iter") <- iter
  attr(tab, "cor_method") <- cor_method
  attr(tab, "method") <- method
  attr(tab, "n_cases") <- n
  attr(tab, "n_edges") <- length(orig)
  attr(tab, "n_success") <- n_success
  attr(tab, "n_failed") <- n_failed
  class(tab) <- c("psychnet_casedrop", "data.frame")
  tab
}

#' @describeIn net_casedrop_reliability The metric table and the edge-weight
#'   CS-coefficient.
#' @export
print.psychnet_casedrop <- function(x, ...) {
  cat(sprintf("# edge-weight stability: %s | CS = %.2f (%s cor > %.2f in > %.0f%% of subsets)\n",
              attr(x, "method"), attr(x, "cs"), attr(x, "cor_method"),
              attr(x, "threshold"), 100 * attr(x, "certainty")))
  print(`class<-`(x, "data.frame"))
  invisible(x)
}

#' @describeIn net_casedrop_reliability Similarity metrics by drop proportion.
#' @export
plot.psychnet_casedrop <- function(x, ...) {
  panels <- c("correlation", "mean_abs_dev", "median_abs_dev", "max_abs_dev")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  tab <- `class<-`(x, "data.frame")
  graphics::par(mfrow = c(2, 2), mar = c(4, 4.2, 3, 1))
  for (nm in panels) {
    d <- tab[tab$metric == nm, , drop = FALSE]
    d <- d[order(d$drop_prop), , drop = FALSE]
    yl <- if (nm == "correlation") c(0, 1)
          else c(0, max(d$mean + d$sd, na.rm = TRUE))
    graphics::plot(d$drop_prop, d$mean, type = "n", ylim = yl,
                   xlab = "proportion of cases dropped", ylab = nm,
                   main = nm, bty = "n")
    graphics::polygon(c(d$drop_prop, rev(d$drop_prop)),
                      c(pmax(yl[1], d$mean - d$sd), rev(d$mean + d$sd)),
                      col = grDevices::adjustcolor(.psn_pal$sig, 0.18),
                      border = NA)
    graphics::lines(d$drop_prop, d$mean, col = .psn_pal$sig, lwd = 2)
    graphics::points(d$drop_prop, d$mean, col = .psn_pal$sig, pch = 19, cex = 0.8)
    if (nm == "correlation") {
      graphics::abline(h = attr(x, "threshold"), lty = 2, col = .psn_pal$ref)
      graphics::mtext(sprintf("CS = %.2f", attr(x, "cs")), side = 3,
                      line = -1.3, cex = 0.8)
    }
  }
  invisible(x)
}

# ---- split-half reliability ------------------------------------------------

#' Split-half reliability of the network edge structure
#'
#' Repeatedly splits the sample into two halves, estimates a network on each, and
#' compares their edge-weight vectors. Reports, across splits, the edge-weight
#' correlation between halves plus the mean/median/maximum absolute edge
#' deviation - a psychometric reliability view of the estimated structure.
#'
#' @param data Data frame or matrix (rows = observations), resampled exactly as
#'   given (see [net_boot()]), or a `psychnet_group` (split-half per level).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param iter Number of split-half iterations. Default 100.
#' @param split Fraction of rows in the first half. Default 0.5.
#' @param cor_method Correlation method for the between-halves edge comparison:
#'   `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param labels Optional node labels.
#' @param estimator_args Named list of estimator arguments. `cor_method`,
#'   `iter` and `split` are arguments of this diagnostic (the between-halves
#'   comparison), so an estimator's own `cor_method` or `threshold` must go
#'   here, e.g. `estimator_args = list(cor_method = "spearman")`; values in
#'   `...` with those names never reach the estimator.
#' @param ... Passed to the estimator by `net_split_reliability()`; ignored by
#'   the methods.
#' @param x A `psychnet_reliability` object.
#' @return A tidy `data.frame` (class `psychnet_reliability`), one row per metric
#'   with columns `metric`, `mean`, `sd`, `lower`, `upper` (2.5% and 97.5%
#'   percentiles across splits) and `n_valid` (splits in which the metric was
#'   defined). The edge correlation is 0 for a split in which exactly one half
#'   has a constant edge vector (for example an empty network against a
#'   non-empty one, as in bootnet's `cor0()`), and undefined when both halves
#'   are constant (they agree, but no correlation exists) or fewer than two
#'   edges are available; undefined splits are excluded from the summary, counted out of
#'   `n_valid`, and announced by a `psychnet_undefined_statistic` warning. A
#'   metric defined in no split has `NA` summaries. A failed split fit raises a
#'   `psychnet_resample_failed` warning carrying the first error; a failed
#'   full-sample fit raises a `psychnet_estimation_failed` error with the
#'   estimator's message. The per-split draws are stored with the result for
#'   `plot()`, which draws a histogram of each metric across splits with its
#'   mean marked. `print()` and `plot()` return `x` invisibly.
#' @examples
#' # `method` and `iter` are chosen so the example runs quickly; the defaults
#' # (method = "glasso", iter = 100) are what a real assessment should use.
#' sh <- net_split_reliability(SRL_Claude, method = "pcor", iter = 10)
#' sh
#' plot(sh)
#' @export
net_split_reliability <- function(data, method = "glasso", iter = 100L,
                                split = 0.5,
                                cor_method = c("pearson", "spearman", "kendall"),
                                labels = NULL, estimator_args = list(), ...) {
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_split_reliability, "net_split_reliability",
      "psychnet_reliability_group",
      list(iter = iter, split = split, cor_method = cor_method,
           estimator_args = .resample_dots(list(...), estimator_args))))
  }
  cor_method <- .psn_match_arg(cor_method)
  .psn_check_number(iter, "iter", lower = 1)
  .psn_check_number(split, "split", lower = 0, upper = 1, open = TRUE)
  iter <- as.integer(iter)
  inp <- .resample_input(data, method,
                         .resample_dots(list(...), estimator_args))
  mat <- inp$data; dots <- inp$dots
  n <- nrow(mat); n_half <- max(2L, round(n * split))
  if (n_half < 2L || n - n_half < 2L) {
    .psn_stop("`split` must leave at least two observations in each half.",
         call. = FALSE)
  }
  # The node set both halves are re-expanded to is whatever the estimator makes
  # of the full data (a factor column is a node for mgm, dropped by glasso), so
  # it is read off one full-sample fit rather than off the input columns.
  if (is.null(labels)) {
    full <- .psn_fit_or_stop(.psn_refit(mat, method, NULL, dots),
                             "Full-sample estimation")
    labels <- rownames(full$weights)
  }

  metric_names <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")
  # Halves are drawn in the parent, one per iteration, before any fitting.
  idx_list <- lapply(seq_len(iter), function(i)
    sample.int(n, n_half, replace = FALSE))
  draws <- lapply(idx_list, function(idx_a) .psn_try_draw({
    a <- .psn_refit(mat[idx_a, , drop = FALSE], method, labels, dots)
    b <- .psn_refit(mat[-idx_a, , drop = FALSE], method, labels, dots)
    # Two empty halves agree perfectly, so a constant edge vector leaves the
    # correlation undefined here rather than 0 (the CS rule's convention).
    as.list(.psn_edge_metrics(.psn_edge_vector(a$weights, a$directed),
                              .psn_edge_vector(b$weights, b$directed),
                              cor_method, constant = "na"))
  }))
  fit_ok <- .psn_check_draws(draws, "split-half iterations")
  iters <- as.data.frame(stats::setNames(lapply(metric_names, function(nm)
    vapply(draws, function(d) if (.psn_draw_ok(d)) d[[nm]] else NA_real_,
           numeric(1))), metric_names))
  n_failed <- sum(!fit_ok)

  summ <- do.call(rbind, lapply(metric_names, function(m) {
    s <- .psn_summarise_draws(iters[[m]])
    data.frame(metric = m, mean = s[["mean"]], sd = s[["sd"]],
               lower = s[["lower"]], upper = s[["upper"]],
               n_valid = as.integer(s[["n_valid"]]),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  undefined <- summ$n_valid < sum(fit_ok)
  if (any(undefined)) {
    .psn_warn(sprintf(paste0(
      "%s undefined in %s of %d estimated splits (a half with a constant edge ",
      "vector, e.g. an empty network, or fewer than two edges); the summary ",
      "uses the defined splits only (see `n_valid`)."),
      paste(summ$metric[undefined], collapse = ", "),
      paste(sum(fit_ok) - summ$n_valid[undefined], collapse = ", "),
      sum(fit_ok)), class = "psychnet_undefined_statistic")
  }

  # The result IS the tidy per-metric summary; the raw per-split draws ride
  # along in an attribute for plot().
  attr(summ, "iterations") <- iters
  attr(summ, "iter") <- iter
  attr(summ, "split") <- split
  attr(summ, "cor_method") <- cor_method
  attr(summ, "method") <- method
  attr(summ, "n") <- n
  attr(summ, "n_success") <- sum(fit_ok)
  attr(summ, "n_failed") <- n_failed
  class(summ) <- c("psychnet_reliability", "data.frame")
  summ
}

#' @describeIn net_split_reliability The between-halves metric table.
#' @export
print.psychnet_reliability <- function(x, ...) {
  cat(sprintf("# split-half reliability: %s | %d iterations (%.0f/%.0f split)\n",
              attr(x, "method"), attr(x, "iter"),
              100 * attr(x, "split"), 100 * (1 - attr(x, "split"))))
  print(`class<-`(x, "data.frame"))
  invisible(x)
}

#' @describeIn net_split_reliability Histograms of the metrics across splits.
#' @export
plot.psychnet_reliability <- function(x, ...) {
  it <- attr(x, "iterations")
  labs <- c(correlation = "edge correlation", mean_abs_dev = "mean |edge diff|",
            median_abs_dev = "median |edge diff|", max_abs_dev = "max |edge diff|")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (m in c("correlation", "mean_abs_dev", "median_abs_dev", "max_abs_dev")) {
    v <- it[[m]]
    v <- v[is.finite(v)]
    if (!length(v)) {
      graphics::plot.new()
      graphics::title(main = labs[[m]])
      graphics::text(0.5, 0.5, "unavailable")
      next
    }
    graphics::hist(v, breaks = 25, col = .psn_pal$nonsig, border = "white",
                   main = labs[[m]], xlab = labs[[m]])
    graphics::abline(v = mean(v, na.rm = TRUE), col = .psn_pal$pos, lwd = 2)
  }
  invisible(x)
}
