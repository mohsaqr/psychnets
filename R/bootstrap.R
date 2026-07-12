# Nonparametric bootstrap for edge-weight and centrality accuracy, adapted from
# Nestimate's boot_glasso into a general resample-and-re-estimate loop over any
# psychnet estimator (Epskamp, Borsboom & Fried 2018). The raw resample draws
# are retained on the object so the bootstrapped difference tests of
# difference_test() can read the paired draws.

# Resolve the worker count: NULL -> two thirds of the detected cores, capped at
# the number available; forced to 1 on Windows, where mclapply cannot fork.
.resolve_cores <- function(cores) {
  if (.Platform$OS.type == "windows") return(1L)
  nc <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  if (is.na(nc) || nc < 1L) nc <- 1L
  # CRAN's check harness forbids using more than two cores.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) nc <- min(nc, 2L)
  if (is.null(cores)) return(max(1L, as.integer(floor(nc * 2 / 3))))
  cores <- as.integer(cores)
  stopifnot(length(cores) == 1L, !is.na(cores), cores >= 1L)
  min(cores, nc)
}

# Closed-form GGM node predictability from a precision matrix: R^2 = 1 - 1/Kii,
# clamped to [0, 1] (Haslbeck & Waldorp 2018). NULL when no precision is stored.
#' @noRd
.psn_predictability <- function(fit) {
  K <- fit$precision
  if (is.null(K)) return(NULL)
  pmin(pmax(1 - 1 / diag(as.matrix(K)), 0), 1)
}

# Pairwise difference p-value matrix from stored draws (columns = items): the
# two-sided bootstrap p, 2 * min(P(diff > 0), P(diff < 0)), optionally adjusted
# for multiplicity over the unique pairs. Returns NULL past `max_items`.
#' @noRd
.psn_diff_pmat <- function(draws, labels, p_adjust = "none", max_items = 500L) {
  m <- ncol(draws)
  if (m > max_items) return(NULL)
  p_mat <- matrix(1, m, m, dimnames = list(labels, labels))
  bm <- draws[stats::complete.cases(draws), , drop = FALSE]
  if (nrow(bm) < 2L) { diag(p_mat) <- 0; return(p_mat) }
  for (i in seq_len(m - 1L)) {
    js <- seq.int(i + 1L, m)
    d <- bm[, i] - bm[, js, drop = FALSE]
    pv <- 2 * pmin(colMeans(d > 0), colMeans(d < 0))
    p_mat[i, js] <- pv; p_mat[js, i] <- pv
  }
  if (p_adjust != "none") {
    ut <- upper.tri(p_mat)
    adj <- stats::p.adjust(p_mat[ut], method = p_adjust)
    p_mat[ut] <- adj
    p_mat[lower.tri(p_mat)] <- t(p_mat)[lower.tri(p_mat)]
  }
  diag(p_mat) <- 0
  p_mat
}

#' Bootstrap a psychometric network
#'
#' Resamples observations with replacement, re-estimates the network on each
#' resample, and summarizes the sampling distribution of every edge weight and
#' node centrality (mean, percentile confidence interval, and edge inclusion
#' proportion). An edge is flagged `significant` when its percentile interval
#' excludes zero. The raw per-resample draws are stored on the returned object
#' for use by [difference_test()].
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator (see [psychnet()]). Default `"glasso"`.
#' @param n_boot Number of bootstrap resamples. Default 1000.
#' @param ci Confidence level for percentile intervals. Default 0.95.
#' @param measures Centrality measures to bootstrap. Defaults to the two
#'   recommended for psychometric networks (`"strength"`,
#'   `"expected_influence"`); `"betweenness"`/`"closeness"` and custom measures
#'   (via `centrality_fn`) are also accepted. See [net_centralities()].
#' @param centrality_fn Optional function supplying any non-built-in `measures`
#'   (see [net_centralities()]).
#' @param predictability Logical; if `TRUE` and the estimator returns a
#'   precision matrix (GGM family), bootstrap node predictability (R^2) and
#'   report its interval. Default `FALSE`.
#' @param threshold Logical; if `TRUE`, also return the observed network with
#'   every edge whose bootstrap interval includes zero set to zero
#'   (`$thresholded`). Default `FALSE`.
#' @param diff_test Logical; if `TRUE`, also return two-sided bootstrap
#'   difference p-value matrices for edges (`$edge_diff_p`, `NULL` past 500
#'   edges) and for each centrality measure (`$centrality_diff_p`). Default
#'   `FALSE`.
#' @param p_adjust Multiple-comparison adjustment applied to the difference
#'   p-value matrices (any [stats::p.adjust] method). Default `"none"`.
#' @param labels Optional node labels.
#' @param cores Number of CPU cores for the resample loop. `NULL` (default) uses
#'   two thirds of the detected cores; `1` forces a serial run. Parallelism uses
#'   forking (`parallel::mclapply`) and falls back to serial on Windows. Because
#'   every resample index is drawn in the parent process before any fitting, the
#'   result is identical for any number of cores and reproducible from
#'   `set.seed()`.
#' @param engine Optional estimator engine forwarded to each resample fit
#'   (e.g. `"base"`/`"glasso"` for glasso, `"base"`/`"glmnet"` for ising/mgm).
#'   `NULL` (default) uses the estimator's own default.
#' @param ... Passed to the estimator.
#' @return An object of class `psychnet_bootstrap`: tidy `$edges` (with a
#'   `significant` flag) and `$centrality` data frames, the observed network in
#'   `$observed`, raw resample draws in `$edge_boot`, `$str_boot`, `$ei_boot`,
#'   and the general `$centrality_boot` (named list, one matrix per measure).
#'   Optional `$predictability`, `$thresholded`, `$edge_diff_p`,
#'   `$centrality_diff_p`, plus `$lambda_path`/`$lambda_selected` when the
#'   estimator reports them.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- net_boot(x, n_boot = 50, cores = 1)
#' as.data.frame(bs)
#' @export
net_boot <- function(data, method = "glasso", n_boot = 1000L,
                     ci = 0.95, measures = c("strength", "expected_influence"),
                     centrality_fn = NULL, predictability = FALSE,
                     threshold = FALSE, diff_test = FALSE, p_adjust = "none",
                     labels = NULL, cores = NULL, engine = NULL, ...) {
  # Group object -> bootstrap each level from its stored cross-sectional data,
  # reproducing the group's estimator configuration (incl. labels).
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_boot, "net_boot",
      "psychnet_bootstrap_group",
      list(n_boot = n_boot, ci = ci, measures = measures,
           centrality_fn = centrality_fn, predictability = predictability,
           threshold = threshold, diff_test = diff_test, p_adjust = p_adjust,
           cores = cores, engine = engine)))
  }
  stopifnot(is.numeric(n_boot), length(n_boot) == 1L, is.finite(n_boot),
            n_boot >= 1, ci > 0, ci < 1)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  n_boot <- as.integer(n_boot)   # a fractional count corrupts the stored %d field
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  if (is.null(labels)) labels <- colnames(mat)

  # Engine is forwarded only when set, so estimators without an engine argument
  # (cor/pcor) are unaffected unless the caller explicitly asks for one.
  dots <- list(...)
  if (!is.null(engine)) dots$engine <- engine
  fit_net <- function(m)
    do.call(psychnet, c(list(m, method = method, labels = labels), dots))

  obs <- fit_net(mat)
  p <- nrow(obs$nodes)
  if (is.null(labels)) labels <- obs$nodes$label
  # Directed estimators (e.g. relimp) have an asymmetric network: take every
  # off-diagonal cell, not just the upper triangle.
  ut <- if (isTRUE(obs$directed)) row(obs$weights) != col(obs$weights)
        else upper.tri(obs$weights)
  obs_edges <- obs$weights[ut]
  obs_cent  <- net_centralities(obs, measures = measures,
                                centrality_fn = centrality_fn)
  obs_pred  <- if (predictability) .psn_predictability(obs) else NULL
  if (predictability && is.null(obs_pred)) {
    warning("`predictability = TRUE` ignored: estimator '", obs$method,
            "' returns no precision matrix.", call. = FALSE)
  }
  cores <- .resolve_cores(cores)

  # Draw every resample index in the parent, so the fits are pure functions of
  # the data subset; the result is then independent of the number of cores.
  idx_list <- lapply(seq_len(n_boot),
                     function(b) sample.int(n, n, replace = TRUE))

  one <- function(idx) {
    fit <- tryCatch(fit_net(mat[idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ct <- net_centralities(fit, measures = measures,
                           centrality_fn = centrality_fn)
    list(edge = fit$weights[ut],
         cent = ct[measures],
         pred = if (!is.null(obs_pred)) .psn_predictability(fit) else NULL)
  }

  draws <- if (cores > 1L)
    parallel::mclapply(idx_list, one, mc.cores = cores)
  else lapply(idx_list, one)

  alpha <- (1 - ci) / 2
  edge_boot <- matrix(NA_real_, n_boot, length(obs_edges))
  cent_boot <- stats::setNames(
    replicate(length(measures), matrix(NA_real_, n_boot, p), simplify = FALSE),
    measures)
  pred_boot <- if (!is.null(obs_pred)) matrix(NA_real_, n_boot, p) else NULL
  for (b in seq_len(n_boot)) {
    d <- draws[[b]]
    if (is.null(d) || !is.list(d)) next   # NULL fit or a crashed worker
    edge_boot[b, ] <- d$edge
    for (m in measures) cent_boot[[m]][b, ] <- d$cent[[m]]
    if (!is.null(pred_boot)) pred_boot[b, ] <- d$pred
  }

  qci <- function(v) stats::quantile(v, c(alpha, 1 - alpha), na.rm = TRUE,
                                     names = FALSE)
  ij <- which(ut, arr.ind = TRUE)
  edge_ci <- t(apply(edge_boot, 2L, qci))
  edges <- data.frame(
    from = labels[ij[, 1L]], to = labels[ij[, 2L]],
    observed = obs_edges,
    mean = colMeans(edge_boot, na.rm = TRUE),
    lower = edge_ci[, 1L], upper = edge_ci[, 2L],
    prop_nonzero = colMeans(abs(edge_boot) > 1e-12, na.rm = TRUE),
    significant = edge_ci[, 1L] > 0 | edge_ci[, 2L] < 0,
    stringsAsFactors = FALSE, row.names = NULL)

  # Centrality table: node + observed/lower/upper per measure, uniform naming.
  cent <- data.frame(node = labels, stringsAsFactors = FALSE, row.names = NULL)
  for (m in measures) {
    mci <- t(apply(cent_boot[[m]], 2L, qci))
    cent[[m]] <- obs_cent[[m]]
    cent[[paste0(m, "_lower")]] <- mci[, 1L]
    cent[[paste0(m, "_upper")]] <- mci[, 2L]
  }

  # Backward-compatible raw-draw aliases for the two default measures.
  str_boot <- if ("strength" %in% measures) cent_boot[["strength"]] else NULL
  ei_boot  <- if ("expected_influence" %in% measures)
    cent_boot[["expected_influence"]] else NULL

  out <- list(observed = obs, edges = edges, centrality = cent,
              edge_boot = edge_boot, str_boot = str_boot, ei_boot = ei_boot,
              centrality_boot = cent_boot,
              edge_labels = paste(edges$from, edges$to, sep = "--"),
              node_labels = labels, measures = measures,
              n_boot = n_boot, ci = ci, method = obs$method,
              lambda_path = obs$lambda_path, lambda_selected = obs$lambda)

  if (!is.null(pred_boot)) {
    pci <- t(apply(pred_boot, 2L, qci))
    out$predictability <- data.frame(
      node = labels, value = obs_pred,
      lower = pci[, 1L], upper = pci[, 2L],
      stringsAsFactors = FALSE, row.names = NULL)
    out$predictability_boot <- pred_boot
  }

  if (threshold) {
    th <- obs$weights
    not_sig <- !edges$significant
    bad <- ij[not_sig, , drop = FALSE]
    th[bad] <- 0
    if (!isTRUE(obs$directed)) th[bad[, c(2L, 1L), drop = FALSE]] <- 0
    out$thresholded <- th
  }

  if (diff_test) {
    out$edge_diff_p <- .psn_diff_pmat(edge_boot, out$edge_labels, p_adjust)
    out$centrality_diff_p <- lapply(measures, function(m)
      .psn_diff_pmat(cent_boot[[m]], labels, p_adjust))
    names(out$centrality_diff_p) <- measures
  }

  structure(out, class = "psychnet_bootstrap")
}

#' Bootstrapped difference test for edges or centralities
#'
#' Tests, within a single network, whether two edge weights or two node
#' centralities differ. For every pair it forms the per-resample difference from
#' the stored bootstrap draws, takes the percentile interval of that difference,
#' and flags the pair `significant` when the interval excludes zero; it also
#' reports the two-sided bootstrap p-value (Epskamp, Borsboom & Fried 2018).
#' This is the within-network counterpart to the edge accuracy intervals
#' reported by [net_boot()].
#'
#' @param boot A `psychnet_bootstrap` object from [net_boot()].
#' @param type Quantity to compare: `"edge"` (default), or any centrality
#'   measure bootstrapped by [net_boot()] (e.g. `"strength"`,
#'   `"expected_influence"`).
#' @param ci Confidence level for the difference interval. Defaults to the level
#'   used by the bootstrap object.
#' @param p_adjust Multiple-comparison adjustment for the pairwise p-values (any
#'   [stats::p.adjust] method). Default `"none"`.
#' @return A tidy data frame, one row per pair, with `item1`, `item2`, the two
#'   observed values, their observed difference, the percentile interval of the
#'   bootstrap difference (`lower`, `upper`), the two-sided `p_value`, and a
#'   logical `significant`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
#' difference_test(bs, type = "strength")
#' @export
difference_test <- function(boot, type = "edge", ci = NULL,
                            p_adjust = "none") {
  stopifnot(inherits(boot, "psychnet_bootstrap"))
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  draws <- if (type == "edge") boot$edge_boot else boot$centrality_boot[[type]]
  if (is.null(draws)) {
    stop("No stored draws for type = '", type, "'. Available: edge, ",
         paste(boot$measures, collapse = ", "), ".", call. = FALSE)
  }
  labs <- if (type == "edge") boot$edge_labels else boot$node_labels
  obs <- if (type == "edge") boot$edges$observed else boot$centrality[[type]]
  if (is.null(ci)) ci <- boot$ci
  stopifnot(ci > 0, ci < 1)
  alpha <- (1 - ci) / 2

  m <- ncol(draws)
  if (m < 2L) stop("Need at least two items to compare.")
  pairs <- which(upper.tri(matrix(0, m, m)), arr.ind = TRUE)
  i <- pairs[, 1L]; j <- pairs[, 2L]
  stats_k <- t(vapply(seq_len(nrow(pairs)), function(k) {
    d <- draws[, i[k]] - draws[, j[k]]
    c(stats::quantile(d, c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE),
      2 * min(mean(d > 0, na.rm = TRUE), mean(d < 0, na.rm = TRUE)))
  }, numeric(3L)))
  pval <- stats_k[, 3L]
  if (p_adjust != "none") pval <- stats::p.adjust(pval, method = p_adjust)

  out <- data.frame(
    item1 = labs[i], item2 = labs[j],
    value1 = obs[i], value2 = obs[j],
    obs_diff = obs[i] - obs[j],
    lower = stats_k[, 1L], upper = stats_k[, 2L],
    p_value = pval,
    significant = stats_k[, 1L] > 0 | stats_k[, 2L] < 0,
    stringsAsFactors = FALSE, row.names = NULL)
  # `type` (and the observed per-item values) let plot.psychnet_difference draw
  # the bootnet-style significance-box matrix without re-deriving anything.
  attr(out, "diff_type") <- type
  attr(out, "observed") <- stats::setNames(obs, labs)
  class(out) <- c("psychnet_difference", "data.frame")
  out
}

#' Tidy a network bootstrap
#'
#' @param x A `psychnet_bootstrap` object.
#' @param row.names,optional Ignored (S3 consistency).
#' @param ... Unused.
#' @param significant If `TRUE`, return only the edges whose confidence interval
#'   excludes zero. Default `FALSE` (all edges).
#' @return The tidy `$edges` data frame (one row per edge, with its percentile
#'   interval, inclusion proportion, and `significant` flag).
#' @export
as.data.frame.psychnet_bootstrap <- function(x, row.names = NULL,
                                             optional = FALSE, ...,
                                             significant = FALSE) {
  ed <- x$edges
  if (isTRUE(significant)) ed <- ed[ed$significant %in% TRUE, , drop = FALSE]
  ed
}

#' Print a network bootstrap
#'
#' @param x A `psychnet_bootstrap` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_bootstrap <- function(x, ...) {
  cat(sprintf("<psychnet_bootstrap> %s, %d resamples, %.0f%% CI\n",
              x$method, x$n_boot, 100 * x$ci))
  cat(sprintf("  %d edges (%d significant), %d nodes, measures: %s\n",
              nrow(x$edges), sum(x$edges$significant, na.rm = TRUE),
              nrow(x$centrality), paste(x$measures, collapse = ", ")))
  invisible(x)
}
