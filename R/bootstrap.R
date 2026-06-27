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
#' @param method Estimator (see [psychnet()]). Default `"EBICglasso"`.
#' @param n_boot Number of bootstrap resamples. Default 1000.
#' @param ci Confidence level for percentile intervals. Default 0.95.
#' @param labels Optional node labels.
#' @param cores Number of CPU cores for the resample loop. `NULL` (default) uses
#'   two thirds of the detected cores; `1` forces a serial run. Parallelism uses
#'   forking (`parallel::mclapply`) and falls back to serial on Windows. Because
#'   every resample index is drawn in the parent process before any fitting, the
#'   result is identical for any number of cores and reproducible from
#'   `set.seed()`.
#' @param ... Passed to the estimator.
#' @return An object of class `psychnet_bootstrap`: tidy `$edges` (with a
#'   `significant` flag) and `$centrality` data frames, the observed network in
#'   `$observed`, and the raw resample draws in `$edge_boot`, `$str_boot`,
#'   `$ei_boot` (one row per resample) with their labels in `$edge_labels` and
#'   `$node_labels`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- bootstrap_network(x, n_boot = 50, cores = 1)
#' as.data.frame(bs)
#' @export
bootstrap_network <- function(data, method = "EBICglasso", n_boot = 1000L,
                              ci = 0.95, labels = NULL, cores = NULL, ...) {
  stopifnot(is.numeric(n_boot), length(n_boot) == 1L, is.finite(n_boot),
            n_boot >= 1, ci > 0, ci < 1)
  n_boot <- as.integer(n_boot)   # a fractional count corrupts the stored %d field
  mat <- .as_numeric_matrix(data)
  n <- nrow(mat)
  if (is.null(labels)) labels <- colnames(mat)

  obs <- psychnet(mat, method = method, labels = labels, ...)
  p <- nrow(obs$nodes)
  # Directed estimators (e.g. relimp) have an asymmetric network: take every
  # off-diagonal cell, not just the upper triangle.
  ut <- if (isTRUE(obs$directed)) row(obs$weights) != col(obs$weights)
        else upper.tri(obs$weights)
  obs_edges <- obs$weights[ut]
  obs_cent  <- centrality(obs)
  cores <- .resolve_cores(cores)

  # Draw every resample index in the parent, so the fits are pure functions of
  # the data subset; the result is then independent of the number of cores.
  idx_list <- lapply(seq_len(n_boot),
                     function(b) sample.int(n, n, replace = TRUE))

  one <- function(idx) {
    fit <- tryCatch(
      psychnet(mat[idx, , drop = FALSE], method = method,
                       labels = labels, ...),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ct <- centrality(fit)
    list(edge = fit$weights[ut], strength = ct$strength,
         ei = ct$expected_influence)
  }

  draws <- if (cores > 1L)
    parallel::mclapply(idx_list, one, mc.cores = cores)
  else lapply(idx_list, one)

  alpha <- (1 - ci) / 2
  edge_boot <- matrix(NA_real_, n_boot, length(obs_edges))
  str_boot  <- matrix(NA_real_, n_boot, p)
  ei_boot   <- matrix(NA_real_, n_boot, p)
  for (b in seq_len(n_boot)) {
    d <- draws[[b]]
    if (is.null(d) || !is.list(d)) next   # NULL fit or a crashed worker
    edge_boot[b, ] <- d$edge
    str_boot[b, ]  <- d$strength
    ei_boot[b, ]   <- d$ei
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

  str_ci <- t(apply(str_boot, 2L, qci))
  ei_ci  <- t(apply(ei_boot, 2L, qci))
  cent <- data.frame(
    node = labels,
    strength = obs_cent$strength,
    strength_lower = str_ci[, 1L], strength_upper = str_ci[, 2L],
    expected_influence = obs_cent$expected_influence,
    ei_lower = ei_ci[, 1L], ei_upper = ei_ci[, 2L],
    stringsAsFactors = FALSE, row.names = NULL)

  structure(list(observed = obs, edges = edges, centrality = cent,
                 edge_boot = edge_boot, str_boot = str_boot, ei_boot = ei_boot,
                 edge_labels = paste(edges$from, edges$to, sep = "--"),
                 node_labels = labels,
                 n_boot = n_boot, ci = ci, method = obs$method),
            class = "psychnet_bootstrap")
}

#' Bootstrapped difference test for edges or centralities
#'
#' Tests, within a single network, whether two edge weights or two node
#' centralities differ. For every pair it forms the per-resample difference from
#' the stored bootstrap draws, takes the percentile interval of that difference,
#' and flags the pair `significant` when the interval excludes zero (Epskamp,
#' Borsboom & Fried 2018). This is the within-network counterpart to the edge
#' accuracy intervals reported by [bootstrap_network()].
#'
#' @param boot A `psychnet_bootstrap` object from [bootstrap_network()].
#' @param type Quantity to compare: `"edge"` (default), `"strength"`, or
#'   `"expected_influence"`.
#' @param ci Confidence level for the difference interval. Defaults to the level
#'   used by the bootstrap object.
#' @return A tidy data frame, one row per pair, with `item1`, `item2`, the two
#'   observed values, their observed difference, the percentile interval of the
#'   bootstrap difference (`lower`, `upper`), and a logical `significant`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' bs <- bootstrap_network(x, n_boot = 100)
#' difference_test(bs, type = "strength")
#' @export
difference_test <- function(boot, type = c("edge", "strength",
                                           "expected_influence"), ci = NULL) {
  stopifnot(inherits(boot, "psychnet_bootstrap"))
  type <- match.arg(type)
  draws <- switch(type, edge = boot$edge_boot, strength = boot$str_boot,
                  expected_influence = boot$ei_boot)
  if (is.null(draws))
    stop("This bootstrap object has no stored draws; re-run bootstrap_network().")
  labs <- if (type == "edge") boot$edge_labels else boot$node_labels
  obs <- switch(type, edge = boot$edges$observed,
                strength = boot$centrality$strength,
                expected_influence = boot$centrality$expected_influence)
  if (is.null(ci)) ci <- boot$ci
  stopifnot(ci > 0, ci < 1)
  alpha <- (1 - ci) / 2

  m <- ncol(draws)
  if (m < 2L) stop("Need at least two items to compare.")
  pairs <- which(upper.tri(matrix(0, m, m)), arr.ind = TRUE)
  i <- pairs[, 1L]; j <- pairs[, 2L]
  diff_ci <- t(vapply(seq_len(nrow(pairs)), function(k) {
    d <- draws[, i[k]] - draws[, j[k]]
    stats::quantile(d, c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  }, numeric(2L)))

  data.frame(
    item1 = labs[i], item2 = labs[j],
    value1 = obs[i], value2 = obs[j],
    obs_diff = obs[i] - obs[j],
    lower = diff_ci[, 1L], upper = diff_ci[, 2L],
    significant = diff_ci[, 1L] > 0 | diff_ci[, 2L] < 0,
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Tidy a network bootstrap
#'
#' @param x A `psychnet_bootstrap` object.
#' @param ... Unused.
#' @return The tidy `$edges` data frame (one row per edge, with its percentile
#'   interval, inclusion proportion, and `significant` flag).
#' @export
as.data.frame.psychnet_bootstrap <- function(x, ...) x$edges

#' Print a network bootstrap
#'
#' @param x A `psychnet_bootstrap` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_bootstrap <- function(x, ...) {
  cat(sprintf("<psychnet_bootstrap> %s, %d resamples, %.0f%% CI\n",
              x$method, x$n_boot, 100 * x$ci))
  cat(sprintf("  %d edges (%d significant), %d nodes\n",
              nrow(x$edges), sum(x$edges$significant, na.rm = TRUE),
              nrow(x$centrality)))
  invisible(x)
}
