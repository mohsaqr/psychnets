# Two robustness diagnostics adapted from Nestimate, in the psychnets idiom
# (raw `data` in, re-estimate via psychnet(method=), base R only):
#   * casedrop_reliability() - edge-weight case-dropping stability (the
#     edge-vector complement to net_stability()'s centrality CS-coefficient).
#   * network_reliability()  - split-half reliability of the edge structure.
# Both are estimator-agnostic: they route every refit through psychnet(method=).

# Off-diagonal edge vector of a weight matrix (upper triangle when undirected,
# every off-diagonal cell when directed). Internal.
.psn_edge_vector <- function(weights, directed) {
  mask <- if (isTRUE(directed)) row(weights) != col(weights) else
          upper.tri(weights)
  as.vector(weights[mask])
}

# The four bootnet-style edge-vector similarity metrics. Internal.
.psn_edge_metrics <- function(a, b, cor_method) {
  d <- abs(a - b)
  cor_ok <- stats::sd(a, na.rm = TRUE) > 0 && stats::sd(b, na.rm = TRUE) > 0
  c(mean_abs_dev   = mean(d, na.rm = TRUE),
    median_abs_dev = stats::median(d, na.rm = TRUE),
    correlation    = if (cor_ok) stats::cor(a, b, method = cor_method) else NA_real_,
    max_abs_dev    = max(d, na.rm = TRUE))
}

# Re-estimate and return the weight matrix + directedness (NULL on failure),
# always aligned to the full `labels` set. If a subsample makes a column
# constant, psychnet drops it; we re-expand the result to the full node set with
# the dropped node isolated (zero edges) so the edge vector stays comparable to
# the full-sample one. Without this, such draws would error and be silently
# skipped, biasing the stability/reliability estimate upward.
.psn_refit <- function(mat, method, labels, ...) {
  if (!is.null(labels) && length(labels) == ncol(mat)) colnames(mat) <- labels
  fit <- tryCatch(psychnet(mat, method = method, ...), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  W <- fit$weights
  if (!is.null(labels) && !identical(rownames(W), labels)) {
    full <- matrix(0, length(labels), length(labels),
                   dimnames = list(labels, labels))
    keep <- intersect(rownames(W), labels)
    full[keep, keep] <- W[keep, keep]
    W <- full
  }
  list(weights = W, directed = isTRUE(fit$directed))
}

# ---- edge-weight case-dropping stability -----------------------------------

#' Edge-weight stability coefficient (case-dropping subset bootstrap)
#'
#' The edge-vector complement to [net_stability()]. For each drop proportion the
#' network is re-estimated on random case-dropped subsets and the subset
#' edge-weight vector is compared with the full-sample one. The edge-weight
#' CS-coefficient is the largest drop proportion at which the edge-vector
#' correlation stays `>= threshold` with probability `>= certainty` (Epskamp,
#' Borsboom & Fried 2018).
#'
#' @param data Numeric data frame or matrix (rows = observations), or a
#'   `psychnet_group` (case-dropped per level).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param drop_prop Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.
#' @param iter Subsets per proportion. Default 100.
#' @param threshold Minimum acceptable edge-vector correlation. Default 0.7.
#' @param certainty Probability the correlation must exceed `threshold`.
#'   Default 0.95.
#' @param cor_method Correlation method for the edge-vector comparison:
#'   `"spearman"` (default, robust to the wide range of edge weights),
#'   `"pearson"`, or `"kendall"`.
#' @param labels Optional node labels.
#' @param ... Passed to the estimator.
#' @return A tidy `data.frame` (class `psychnet_casedrop`), one row per metric
#'   per drop proportion, with columns `metric`, `drop_prop`, `mean`, `sd`. The
#'   edge-weight CS-coefficient is carried in `attr(x, "cs")` and shown when the
#'   result is printed. Visualise it with [plot.psychnet_casedrop()].
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy. *Behavior Research Methods*,
#'   50(1), 195-212.
#' @examples
#' # `iter` and `drop_prop` are kept small here so the example runs quickly;
#' # the defaults (iter = 100, drop_prop = seq(0.1, 0.9, 0.1)) are what a real
#' # reliability assessment should use.
#' casedrop_reliability(SRL_Claude, iter = 5, drop_prop = c(0.25, 0.5))
#' @export
casedrop_reliability <- function(data, method = "glasso",
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 iter = 100L, threshold = 0.7, certainty = 0.95,
                                 cor_method = c("spearman", "pearson", "kendall"),
                                 labels = NULL, ...) {
  # Group object -> case-drop each level from its stored cross-sectional data.
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, casedrop_reliability, "casedrop_reliability",
      "psychnet_casedrop_group",
      list(drop_prop = drop_prop, iter = iter, threshold = threshold,
           certainty = certainty, cor_method = cor_method)))
  }
  cor_method <- match.arg(cor_method)
  stopifnot(length(drop_prop) >= 1L, all(drop_prop > 0), all(drop_prop < 1),
            is.numeric(iter), length(iter) == 1L, is.finite(iter), iter >= 1,
            threshold > 0, threshold <= 1, certainty > 0, certainty <= 1)
  iter <- as.integer(iter)
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  if (is.null(labels)) labels <- colnames(mat)

  full <- .psn_refit(mat, method, labels, ...)
  if (is.null(full)) stop("Full-sample estimation failed.", call. = FALSE)
  orig <- .psn_edge_vector(full$weights, full$directed)

  metric_names <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")
  store <- stats::setNames(
    lapply(metric_names, function(.) matrix(NA_real_, iter, length(drop_prop))),
    metric_names)

  for (pj in seq_along(drop_prop)) {
    keep_n <- max(2L, round(n * (1 - drop_prop[pj])))
    if (keep_n >= n) next
    for (it in seq_len(iter)) {
      idx <- sample.int(n, keep_n, replace = FALSE)
      sub <- .psn_refit(mat[idx, , drop = FALSE], method, labels, ...)
      if (is.null(sub)) next
      m <- .psn_edge_metrics(orig, .psn_edge_vector(sub$weights, sub$directed),
                             cor_method)
      for (nm in metric_names) store[[nm]][it, pj] <- m[[nm]]
    }
  }

  tab <- do.call(rbind, lapply(metric_names, function(nm) {
    M <- store[[nm]]
    data.frame(metric = nm, drop_prop = drop_prop,
               mean = colMeans(M, na.rm = TRUE),
               sd = apply(M, 2L, stats::sd, na.rm = TRUE),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  prop_above <- colMeans(store$correlation >= threshold, na.rm = TRUE)
  ok <- which(prop_above >= certainty)
  cs <- if (length(ok)) max(drop_prop[ok]) else 0

  # The result IS the tidy table (one row per metric per drop proportion);
  # the CS-coefficient and raw draws ride along as attributes for plot()/print().
  attr(tab, "cs") <- cs
  attr(tab, "correlations") <- store$correlation
  attr(tab, "threshold") <- threshold
  attr(tab, "certainty") <- certainty
  attr(tab, "iter") <- iter
  attr(tab, "cor_method") <- cor_method
  attr(tab, "method") <- method
  attr(tab, "n_cases") <- n
  attr(tab, "n_edges") <- length(orig)
  class(tab) <- c("psychnet_casedrop", "data.frame")
  tab
}

#' Print an edge-weight stability result
#'
#' @param x A `psychnet_casedrop` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_casedrop <- function(x, ...) {
  cat(sprintf("# edge-weight stability: %s | CS = %.2f (%s cor >= %.2f at %.0f%%)\n",
              attr(x, "method"), attr(x, "cs"), attr(x, "cor_method"),
              attr(x, "threshold"), 100 * attr(x, "certainty")))
  print(`class<-`(x, "data.frame"))
  invisible(x)
}

#' Plot edge-weight case-dropping stability
#'
#' Draws the four similarity metrics (correlation, mean/median/max absolute
#' deviation) against the proportion of cases dropped, each with a +/- 1 SD band;
#' the correlation panel carries the acceptance threshold and CS-coefficient.
#'
#' @param x A `psychnet_casedrop` object.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' # Small `iter` / `drop_prop` for a fast example; see [casedrop_reliability()]
#' # for the defaults a real assessment should use.
#' plot(casedrop_reliability(SRL_Claude, iter = 5, drop_prop = c(0.25, 0.5)))
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
#' @param data Numeric data frame or matrix (rows = observations), or a
#'   `psychnet_group` (split-half per level).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param iter Number of split-half iterations. Default 100.
#' @param split Fraction of rows in the first half. Default 0.5.
#' @param cor_method Correlation method for the between-halves edge comparison:
#'   `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param labels Optional node labels.
#' @param ... Passed to the estimator.
#' @return A tidy `data.frame` (class `psychnet_reliability`), one row per metric
#'   with columns `metric`, `mean`, `sd`, `lower`, `upper`. The per-split draws
#'   are carried in `attr(x, "iterations")` for [plot.psychnet_reliability()].
#' @examples
#' # `iter` is kept small here so the example runs quickly; the default
#' # (iter = 100) is what a real reliability assessment should use.
#' network_reliability(SRL_Claude, iter = 10)
#' @export
network_reliability <- function(data, method = "glasso", iter = 100L,
                                split = 0.5,
                                cor_method = c("pearson", "spearman", "kendall"),
                                labels = NULL, ...) {
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, network_reliability, "network_reliability",
      "psychnet_reliability_group",
      list(iter = iter, split = split, cor_method = cor_method)))
  }
  cor_method <- match.arg(cor_method)
  stopifnot(is.numeric(iter), length(iter) == 1L, is.finite(iter), iter >= 1,
            split > 0, split < 1)
  iter <- as.integer(iter)
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat); n_half <- max(2L, round(n * split))
  if (is.null(labels)) labels <- colnames(mat)

  res <- t(vapply(seq_len(iter), function(i) {
    idx_a <- sample.int(n, n_half, replace = FALSE)
    a <- .psn_refit(mat[idx_a, , drop = FALSE], method, labels, ...)
    b <- .psn_refit(mat[-idx_a, , drop = FALSE], method, labels, ...)
    if (is.null(a) || is.null(b)) return(rep(NA_real_, 4L))
    .psn_edge_metrics(.psn_edge_vector(a$weights, a$directed),
                      .psn_edge_vector(b$weights, b$directed), cor_method)
  }, numeric(4L)))
  colnames(res) <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")
  iters <- as.data.frame(res, stringsAsFactors = FALSE)

  summ <- do.call(rbind, lapply(colnames(res), function(m) {
    v <- iters[[m]]
    data.frame(metric = m, mean = mean(v, na.rm = TRUE),
               sd = stats::sd(v, na.rm = TRUE),
               lower = stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE),
               upper = stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE),
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  # The result IS the tidy per-metric summary; the raw per-split draws ride
  # along in an attribute for plot().
  attr(summ, "iterations") <- iters
  attr(summ, "iter") <- iter
  attr(summ, "split") <- split
  attr(summ, "cor_method") <- cor_method
  attr(summ, "method") <- method
  attr(summ, "n") <- n
  class(summ) <- c("psychnet_reliability", "data.frame")
  summ
}

#' Print a split-half reliability result
#'
#' @param x A `psychnet_reliability` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_reliability <- function(x, ...) {
  cat(sprintf("# split-half reliability: %s | %d iterations (%.0f/%.0f split)\n",
              attr(x, "method"), attr(x, "iter"),
              100 * attr(x, "split"), 100 * (1 - attr(x, "split"))))
  print(`class<-`(x, "data.frame"))
  invisible(x)
}

#' Plot split-half reliability
#'
#' Histograms of the four between-halves edge metrics across split-half
#' iterations, with each observed mean marked.
#'
#' @param x A `psychnet_reliability` object.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' # Small `iter` for a fast example; see [network_reliability()] for the
#' # default a real assessment should use.
#' plot(network_reliability(SRL_Claude, iter = 10))
#' @export
plot.psychnet_reliability <- function(x, ...) {
  it <- attr(x, "iterations")
  labs <- c(correlation = "edge correlation", mean_abs_dev = "mean |edge diff|",
            median_abs_dev = "median |edge diff|", max_abs_dev = "max |edge diff|")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (m in c("correlation", "mean_abs_dev", "median_abs_dev", "max_abs_dev")) {
    v <- it[[m]]
    graphics::hist(v, breaks = 25, col = .psn_pal$nonsig, border = "white",
                   main = labs[[m]], xlab = labs[[m]])
    graphics::abline(v = mean(v, na.rm = TRUE), col = .psn_pal$pos, lwd = 2)
  }
  invisible(x)
}
