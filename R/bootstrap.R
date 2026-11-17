# Nonparametric bootstrap for edge-weight and centrality accuracy, adapted from
# Nestimate's boot_glasso into a general resample-and-re-estimate loop over any
# psychnet estimator (Epskamp, Borsboom & Fried 2018). The raw resample draws
# are retained on the object so the bootstrapped difference tests of
# difference_test() can read the paired draws.

# Resolve the worker count: NULL -> two thirds of the detected cores, capped at
# the number available; forced to 1 on Windows, where mclapply cannot fork.
.resolve_cores <- function(cores) {
  .psn_check_number(cores, "cores", lower = 1, whole = TRUE, null_ok = TRUE)
  if (.Platform$OS.type == "windows") return(1L)
  # detectCores() returns NA (it does not error) when it cannot tell.
  nc <- parallel::detectCores()
  if (is.na(nc) || nc < 1L) nc <- 1L
  # CRAN's check harness forbids using more than two cores.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) nc <- min(nc, 2L)
  if (is.null(cores)) return(max(1L, as.integer(floor(nc * 2 / 3))))
  min(as.integer(cores), as.integer(nc))
}

# Estimators whose edge weights are non-negative by construction (relimp's
# R^2 shares): a percentile interval of such a weight excludes zero almost
# surely, so "significant" (interval excludes 0) carries no information.
#' @noRd
.psn_nonnegative_method <- function(method) {
  identical(method, "relimp")
}

#' @noRd
.psn_warn_uninformative <- function(method, what) {
  .psn_warn(sprintf(paste0(
    "%s is uninformative for method '%s': its edge weights are non-negative ",
    "by construction, so every bootstrap interval excludes 0 and every edge ",
    "is flagged significant."), what, method),
    class = "psychnet_uninformative_test")
}

# Closed-form GGM node predictability from a precision matrix. With a penalized
# diagonal the fitted marginal variance is diag(solve(K)), not the input S.
#' @noRd
.psn_predictability <- function(fit) {
  K <- fit$precision
  if (is.null(K)) return(NULL)
  labs <- rownames(K)
  good <- is.finite(diag(K))
  out <- stats::setNames(rep(NA_real_, nrow(K)), labs)
  if (!any(good)) return(out)
  Ks <- K[good, good, drop = FALSE]
  glabs <- labs[good]
  sdiag <- if (isTRUE(fit$penalize_diagonal)) {
    diag(solve(Ks))
  } else if (!is.null(fit$cor_matrix)) {
    S <- fit$cor_matrix
    if (!is.null(rownames(S)) && all(glabs %in% rownames(S)))
      diag(S[glabs, glabs, drop = FALSE]) else diag(S)
  } else rep(1, nrow(Ks))
  out[good] <- pmin(pmax(1 - 1 / (diag(Ks) * sdiag), 0), 1)
  out
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
    # Split ties evenly between the two tails. This gives p = 1 when every
    # paired difference is exactly zero and avoids anti-conservative tie bias.
    z <- colMeans(d == 0)
    pv <- pmin(1, 2 * pmin(colMeans(d > 0) + z / 2,
                           colMeans(d < 0) + z / 2))
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

#' Bootstrap the edges and centralities of a network
#'
#' @description
#' `net_boot()` bootstraps a network to show how accurately its edge weights and
#' node centralities are estimated. The nonparametric bootstrap approximates the
#' sampling variation of the network, whose edge weights and centralities would
#' differ in another sample from the same population, by re-estimating it on
#' samples drawn with replacement from the observations; the spread of the
#' bootstrapped values shows which edges, and which differences between them,
#' can be interpreted (Epskamp, Borsboom & Fried, 2018).
#'
#' The function resamples the rows of `data`, re-estimates the network with the
#' chosen estimator on every resample, and summarizes each edge weight and each
#' node centrality by its bootstrap mean and percentile confidence interval. For
#' each edge it also reports the proportion of resamples in which the edge was
#' estimated as nonzero, and flags it as `significant` when its interval
#' excludes zero. Options add intervals for node predictability
#' (`predictability`) and a version of the network that keeps only the edges
#' whose interval excludes zero (`threshold`).
#'
#' It returns these intervals as tables read with `as.data.frame()`, and keeps
#' the bootstrapped values on the object, so that [difference_test()] can test
#' differences between edges, or between node centralities, without resampling
#' again.
#'
#' The bootstrap describes the accuracy of estimates within one sample. For the
#' stability of the order of centralities when cases are dropped, use
#' [net_stability()]; for the stability of the edge weights, use
#' [net_casedrop_reliability()] and [net_split_reliability()]; to test whether
#' two groups have different networks, use [net_compare()].
#'
#' @details
#' Each resample is passed to the estimator exactly as the full data would be,
#' including incomplete rows and factor columns, so the estimator's own
#' handling of missing values and variable types applies to every resample.
#' For `method = "mgm"` the node types are decided once on the full data and
#' kept for every resample. A resample on which a variable has no variation
#' keeps that variable as an isolated node.
#'
#' All resample indices are drawn before any network is fitted, so the result
#' is reproducible with [set.seed()] and does not depend on `cores`.
#'
#' For an estimator whose edge weights are non-negative by construction
#' (`method = "relimp"`, whose weights are shares of explained variance), a
#' percentile interval excludes zero almost surely: every edge is then flagged
#' `significant` and `threshold = TRUE` removes nothing. Pairwise difference
#' tests with [difference_test()] remain informative.
#'
#' @param data A data frame or matrix with one row per observation, or a
#'   [psychnet_group] estimated from a data set, in which case each group is
#'   bootstrapped with the group's own estimator settings.
#' @param method The estimator, as in [psychnet()]. Default `"glasso"`.
#' @param n_boot Number of bootstrap resamples. Default `1000`.
#' @param ci Confidence level of the percentile intervals, a number in
#'   (0, 1). Default `0.95`.
#' @param measures Centrality measures to bootstrap. Default
#'   `c("strength", "expected_influence")`; any measure accepted by
#'   [net_centralities()] can be added.
#' @param centrality_fn Optional function supplying measures that are not
#'   built in (see [net_centralities()]).
#' @param predictability If `TRUE`, also bootstrap node predictability (R^2)
#'   when the estimator returns a precision matrix (the Gaussian graphical
#'   models). Default `FALSE`.
#' @param threshold If `TRUE`, also keep the observed network with every edge
#'   whose interval includes zero set to zero. Default `FALSE`.
#' @param diff_test If `TRUE`, also store matrices of pairwise difference
#'   p-values on the object: for the edges (up to 500 edges) and for each
#'   centrality measure. Default `FALSE`. The tidy version of these tests is
#'   [difference_test()].
#' @param p_adjust Multiplicity adjustment for the matrices stored by
#'   `diff_test`, any method of [stats::p.adjust()]. Default `"none"`.
#' @param labels Optional node labels, one per node.
#' @param cores Number of CPU cores. `NULL` (default) uses two thirds of the
#'   detected cores; `1` runs serially. Runs are serial on Windows.
#' @param engine Optional solver for every refit: `"base"` or `"glasso"` for
#'   `method = "glasso"` or `"huge"`, `"base"` or `"glmnet"` for `"ising"` or
#'   `"mgm"`. `NULL` (default) keeps the estimator's default.
#' @param estimator_args Named list of arguments for the estimator. Use it for
#'   an estimator argument that shares a name with an argument of `net_boot()`,
#'   e.g. `estimator_args = list(threshold = TRUE)`. Values here override those
#'   in `...`.
#' @param ... Further arguments for the estimator (`net_boot()`); ignored by
#'   the methods.
#' @param x A `psychnet_bootstrap` object.
#' @param what For `as.data.frame()`: `"edges"` (default), `"centrality"`,
#'   `"predictability"` (needs `predictability = TRUE`) or `"thresholded"`
#'   (needs `threshold = TRUE`).
#' @param significant For `as.data.frame()` with `what = "edges"`: if `TRUE`,
#'   return only the edges whose interval excludes zero. Default `FALSE`.
#' @param row.names,optional Ignored.
#' @param type For `plot()`: `"edges"` (default; edge weights with their
#'   intervals, sorted by observed weight), `"centrality"` (one panel per
#'   measure), `"edge_diff"` or `"centrality_diff"` (the matrix of pairwise
#'   difference tests for the edges or for one measure), or `"predictability"`
#'   (needs `predictability = TRUE`).
#' @param measure For `plot()` with `type = "centrality"` or
#'   `"centrality_diff"`: the measure(s) to draw. Default: all bootstrapped
#'   measures, or the first one for `"centrality_diff"`.
#' @return `net_boot()` returns an object of class `psychnet_bootstrap`. Its
#'   tables are read with `as.data.frame()`:
#'
#'   * `as.data.frame(x)` (`what = "edges"`): one row per edge, with columns
#'     `from`, `to`, `observed`, `mean` (bootstrap mean), `lower`, `upper`
#'     (percentile interval), `prop_nonzero` (share of resamples in which the
#'     edge is nonzero) and `significant` (the interval excludes zero). For a
#'     directed estimator each ordered pair is a row.
#'   * `what = "centrality"`: one row per node, with `node` and, for each
#'     measure, its observed value and the columns `<measure>_lower` and
#'     `<measure>_upper`.
#'   * `what = "predictability"`: one row per node, with `node`, `value` (the
#'     observed R^2), `lower` and `upper`.
#'   * `what = "thresholded"`: one row per edge, with `from`, `to`,
#'     `observed`, `weight` (the observed weight, or 0 when its interval
#'     includes zero) and `retained`.
#'
#'   The object also stores the observed network and the bootstrapped value of
#'   every edge and centrality in every resample, which [difference_test()]
#'   and the difference plots use. `print()` shows the estimator, the number of
#'   resamples, the interval level and the numbers of edges, significant edges
#'   and nodes, and returns `x` invisibly; `plot()` returns `x` invisibly. For
#'   a [psychnet_group], `net_boot()` returns one result per group (see
#'   [psychnet_result_group]).
#' @section Conditions:
#'   `psychnet_bad_input` for an invalid argument, for a group that was not
#'   estimated from a data set, and from `as.data.frame()` or `plot()` when the
#'   requested table was not computed; `psychnet_estimation_failed` when the
#'   network cannot be estimated on the full data (with the estimator's
#'   message); `psychnet_resample_failed` (warning) when some resamples fail,
#'   naming the first error, and as an error when fewer than two succeed;
#'   `psychnet_constant_variable` (warning) when a constant column is dropped;
#'   `psychnet_ignored_argument` (warning) when `predictability = TRUE` but the
#'   estimator returns no precision matrix; `psychnet_uninformative_test`
#'   (warning) for `threshold = TRUE` or `as.data.frame(significant = TRUE)`
#'   with `method = "relimp"`.
#' @seealso [difference_test()] for pairwise differences;
#'   [net_centralities()] for the measures; [net_stability()],
#'   [net_casedrop_reliability()] and [net_split_reliability()] for
#'   case-dropping and split-half stability; [psychnet()] for the estimators.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy: A tutorial paper. *Behavior
#'   Research Methods*, 50(1), 195-212. \doi{10.3758/s13428-017-0862-1}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' # method = "pcor" and n_boot = 50 keep this example fast; a real analysis
#' # uses the "glasso" default and n_boot >= 1000.
#' bs <- net_boot(x, method = "pcor", n_boot = 50, cores = 1)
#' bs
#' as.data.frame(bs)
#' as.data.frame(bs, significant = TRUE)
#' as.data.frame(bs, what = "centrality")
#' plot(bs)                       # edge weights with their intervals
#' plot(bs, type = "centrality")  # centralities with their intervals
#' @export
net_boot <- function(data, method = "glasso", n_boot = 1000L,
                     ci = 0.95, measures = c("strength", "expected_influence"),
                     centrality_fn = NULL, predictability = FALSE,
                     threshold = FALSE, diff_test = FALSE, p_adjust = "none",
                     labels = NULL, cores = NULL, engine = NULL,
                     estimator_args = list(), ...) {
  # Group object -> bootstrap each level from its stored cross-sectional data,
  # reproducing the group's estimator configuration (incl. labels).
  if (inherits(data, "psychnet_group")) {
    return(.group_data_apply(data, net_boot, "net_boot",
      "psychnet_bootstrap_group",
      list(n_boot = n_boot, ci = ci, measures = measures,
           centrality_fn = centrality_fn, predictability = predictability,
           threshold = threshold, diff_test = diff_test, p_adjust = p_adjust,
           cores = cores, engine = engine,
           estimator_args = .resample_dots(list(...), estimator_args))))
  }
  .psn_check_number(n_boot, "n_boot", lower = 1)
  .psn_check_number(ci, "ci", lower = 0, upper = 1, open = TRUE)
  .psn_check_flag(predictability, "predictability")
  .psn_check_flag(threshold, "threshold")
  .psn_check_flag(diff_test, "diff_test")
  if (!is.character(measures) || !length(measures) || anyNA(measures) ||
      anyDuplicated(measures))
    .psn_stop("`measures` must be a character vector of distinct measure names.")
  p_adjust <- .psn_match_arg(p_adjust, stats::p.adjust.methods)
  n_boot <- as.integer(n_boot)   # a fractional count corrupts the stored %d field
  cores <- .resolve_cores(cores)
  # Engine is forwarded only when set, so estimators without an engine argument
  # (cor/pcor) are unaffected unless the caller explicitly asks for one.
  dots <- .resample_dots(list(...), estimator_args)
  native <- .resample_engine(engine, method)
  if (!is.null(native)) dots$native <- native
  inp <- .resample_input(data, method, dots)
  mat <- inp$data; dots <- inp$dots
  n <- nrow(mat)
  fit_net <- function(m, align = NULL)
    .psn_refit_network(m, method = method, labels = align, dots = dots)

  # Node labels come from the observed fit, not the input columns: the
  # estimator decides which columns are nodes (glasso drops a factor column,
  # mgm keeps it), and the draws must align to that.
  obs <- .psn_fit_or_stop(fit_net(mat, labels), "Observed network estimation")
  p <- nrow(obs$nodes)
  if (is.null(labels)) labels <- obs$nodes$label
  # Directed estimators (e.g. relimp) have an asymmetric network: take every
  # off-diagonal cell, not just the upper triangle.
  directed <- isTRUE(obs$directed)
  ut <- if (directed) row(obs$weights) != col(obs$weights)
        else upper.tri(obs$weights)
  obs_edges <- obs$weights[ut]
  obs_cent  <- net_centralities(obs, measures = measures,
                                centrality_fn = centrality_fn)
  obs_pred  <- if (predictability) .psn_predictability(obs) else NULL
  if (predictability && is.null(obs_pred)) {
    .psn_warn("`predictability = TRUE` ignored: estimator '", obs$method,
              "' returns no precision matrix.",
              class = "psychnet_ignored_argument")
  }
  if (threshold && .psn_nonnegative_method(obs$method))
    .psn_warn_uninformative(obs$method, "`threshold = TRUE`")

  # Draw every resample index in the parent, so the fits are pure functions of
  # the data subset; the result is then independent of the number of cores.
  idx_list <- lapply(seq_len(n_boot),
                     function(b) sample.int(n, n, replace = TRUE))

  # The whole draw runs under .psn_try_draw() inside the worker: an error in
  # any step fails this draw only, identically for any number of cores.
  one <- function(idx) .psn_try_draw({
    fit <- fit_net(mat[idx, , drop = FALSE], labels)
    ct <- net_centralities(fit, measures = measures,
                           centrality_fn = centrality_fn)
    list(edge = fit$weights[ut],
         cent = ct[measures],
         pred = if (!is.null(obs_pred)) .psn_predictability(fit) else NULL)
  })

  draws <- if (cores > 1L)
    parallel::mclapply(idx_list, one, mc.cores = cores)
  else lapply(idx_list, one)

  ok <- .psn_check_draws(draws, "bootstrap resamples")
  n_success <- sum(ok)
  n_failed <- n_boot - n_success
  good <- draws[ok]

  # Raw draws: one row per resample, NA rows for failed resamples.
  fill_rows <- function(field, ncols) {
    out <- matrix(NA_real_, n_boot, ncols)
    out[ok, ] <- do.call(rbind, lapply(good, function(d) field(d)))
    out
  }
  edge_boot <- fill_rows(function(d) d$edge, length(obs_edges))
  cent_boot <- lapply(stats::setNames(measures, measures), function(m)
    fill_rows(function(d) d$cent[[m]], p))
  pred_boot <- if (!is.null(obs_pred)) fill_rows(function(d) d$pred, p) else NULL

  alpha <- (1 - ci) / 2
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
  cent <- do.call(cbind, c(
    list(data.frame(node = labels, stringsAsFactors = FALSE)),
    lapply(measures, function(m) {
      mci <- t(apply(cent_boot[[m]], 2L, qci))
      stats::setNames(data.frame(obs_cent[[m]], mci[, 1L], mci[, 2L]),
                      c(m, paste0(m, c("_lower", "_upper"))))
    })))
  rownames(cent) <- NULL

  # Backward-compatible raw-draw aliases for the two default measures.
  str_boot <- if ("strength" %in% measures) cent_boot[["strength"]] else NULL
  ei_boot  <- if ("expected_influence" %in% measures)
    cent_boot[["expected_influence"]] else NULL

  out <- list(observed = obs, edges = edges, centrality = cent,
              edge_boot = edge_boot, str_boot = str_boot, ei_boot = ei_boot,
              centrality_boot = cent_boot,
              edge_labels = paste(edges$from, edges$to,
                                  sep = if (directed) "->" else "--"),
              node_labels = labels, measures = measures,
              n_boot = n_boot, ci = ci, method = obs$method,
              n_success = n_success, n_failed = n_failed,
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
    if (!directed) th[bad[, c(2L, 1L), drop = FALSE]] <- 0
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

#' Test differences between edges or between centralities within a network
#'
#' @description
#' `difference_test()` tests whether pairs of edges, or pairs of nodes on one
#' centrality measure, differ within a single network. It uses the bootstrapped
#' difference test (Epskamp, Borsboom & Fried, 2018), which computes the
#' difference of a pair in every bootstrap resample and treats the pair as
#' different when the resulting interval excludes zero, that is, when the
#' difference exceeds what sampling variation alone would produce. Interpreting
#' one edge as stronger than another, or one node as more central than another,
#' requires such a difference.
#'
#' The function reads the resampled values stored by [net_boot()] and tests
#' every pair of edges (`type = "edge"`) or every pair of nodes on one
#' bootstrapped centrality measure.
#'
#' For each pair it returns the observed difference, the percentile interval of
#' the bootstrapped difference, a two-sided bootstrap p-value and whether the
#' interval excludes zero.
#'
#' The test compares estimates within one network. To compare the same edge
#' or the global strength between two groups, use [net_compare()]; for the
#' interval of each edge or centrality on its own, use [net_boot()].
#'
#' @details
#' The p-value of a pair is twice the smaller of the shares of resamples in
#' which the difference is positive or negative, with resamples in which it is
#' exactly zero split evenly between the two, and capped at 1. `significant`
#' is decided by the interval at level `ci` and is not changed by `p_adjust`,
#' which adjusts only `p_value`.
#'
#' @param boot A `psychnet_bootstrap` object from [net_boot()].
#' @param type What to compare: `"edge"` (default) or one centrality measure
#'   bootstrapped by [net_boot()], such as `"strength"` or
#'   `"expected_influence"`.
#' @param ci Confidence level of the difference interval, a number in (0, 1).
#'   Default: the level used by `boot`.
#' @param p_adjust Multiplicity adjustment for the p-values, any method of
#'   [stats::p.adjust()]. Default `"none"`.
#' @return A data frame of class `psychnet_difference`, one row per pair, with
#'   columns `item1`, `item2`, `value1`, `value2` (the two observed values),
#'   `obs_diff` (`value1 - value2`), `lower`, `upper` (percentile interval of
#'   the bootstrapped difference), `p_value` and `significant` (the interval
#'   excludes zero). Edges are labelled `A--B`, or `A->B` for a directed
#'   estimator. `plot()` draws the table and returns `x` invisibly.
#' @section Conditions:
#'   `psychnet_bad_input` when `boot` is not a `psychnet_bootstrap` object,
#'   `type` is not a bootstrapped quantity, `ci` or `p_adjust` is invalid, or
#'   fewer than two items can be compared.
#' @seealso [net_boot()] for the resamples; [net_compare()] for differences
#'   between groups.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy: A tutorial paper. *Behavior
#'   Research Methods*, 50(1), 195-212. \doi{10.3758/s13428-017-0862-1}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
#' colnames(x) <- paste0("V", 1:5)
#' # method = "pcor" and n_boot = 50 keep this example fast; a real analysis
#' # uses the "glasso" default and n_boot >= 1000.
#' bs <- net_boot(x, method = "pcor", n_boot = 50, cores = 1)
#' dt <- difference_test(bs, type = "strength")
#' dt
#' plot(dt)                   # matrix of pairwise differences
#' plot(dt, style = "forest") # one interval per pair
#' @export
difference_test <- function(boot, type = "edge", ci = NULL,
                            p_adjust = "none") {
  if (!inherits(boot, "psychnet_bootstrap"))
    .psn_stop("`boot` must be a psychnet_bootstrap object from net_boot().")
  available <- c("edge", boot$measures)
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% available) {
    .psn_stop("`type` must be one of: ", paste(available, collapse = ", "),
              " (the bootstrapped measures).")
  }
  p_adjust <- .psn_match_arg(p_adjust, stats::p.adjust.methods)
  draws <- if (identical(type, "edge")) boot$edge_boot
           else boot$centrality_boot[[type]]
  if (is.null(draws)) {
    .psn_stop("No stored draws for type = '", type, "'. Available: edge, ",
         paste(boot$measures, collapse = ", "), ".", call. = FALSE)
  }
  labs <- if (identical(type, "edge")) boot$edge_labels else boot$node_labels
  obs <- if (identical(type, "edge")) boot$edges$observed
         else boot$centrality[[type]]
  ci <- ci %||% boot$ci
  .psn_check_number(ci, "ci", lower = 0, upper = 1, open = TRUE)
  alpha <- (1 - ci) / 2

  m <- ncol(draws)
  if (m < 2L) .psn_stop("Need at least two items to compare.")
  pairs <- which(upper.tri(matrix(0, m, m)), arr.ind = TRUE)
  i <- pairs[, 1L]; j <- pairs[, 2L]
  stats_k <- t(vapply(seq_len(nrow(pairs)), function(k) {
    d <- draws[, i[k]] - draws[, j[k]]
    c(stats::quantile(d, c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE),
      {
        z <- mean(d == 0, na.rm = TRUE)
        min(1, 2 * min(mean(d > 0, na.rm = TRUE) + z / 2,
                       mean(d < 0, na.rm = TRUE) + z / 2))
      })
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

#' @describeIn net_boot The edge intervals (default), the centrality or
#'   predictability intervals, or the thresholded network, chosen by `what`.
#' @export
as.data.frame.psychnet_bootstrap <- function(x, row.names = NULL,
                                             optional = FALSE, ...,
                                             what = c("edges", "centrality",
                                                      "predictability",
                                                      "thresholded"),
                                             significant = FALSE) {
  what <- .psn_match_arg(what)
  .psn_check_flag(significant, "significant")
  if (significant && !identical(what, "edges"))
    .psn_stop("`significant = TRUE` applies to `what = \"edges\"` only.")
  switch(what,
    edges = {
      ed <- x$edges
      if (significant) {
        if (.psn_nonnegative_method(x$method))
          .psn_warn_uninformative(x$method, "`significant = TRUE`")
        ed <- ed[ed$significant %in% TRUE, , drop = FALSE]
        rownames(ed) <- NULL
      }
      ed
    },
    centrality = x$centrality,
    predictability = {
      if (is.null(x$predictability))
        .psn_stop("No predictability stored: re-run net_boot(predictability = ",
                  "TRUE) with an estimator that returns a precision matrix.")
      x$predictability
    },
    thresholded = {
      th <- x$thresholded
      if (is.null(th))
        .psn_stop("No thresholded network stored: re-run net_boot(threshold = TRUE).")
      ut <- if (isTRUE(x$observed$directed)) row(th) != col(th) else upper.tri(th)
      w <- th[ut]
      data.frame(from = x$edges$from, to = x$edges$to,
                 observed = x$edges$observed, weight = w,
                 retained = x$edges$significant %in% TRUE,
                 stringsAsFactors = FALSE, row.names = NULL)
    })
}

#' @describeIn net_boot The estimator, the number of resamples, the interval
#'   level and the edge counts.
#' @export
print.psychnet_bootstrap <- function(x, ...) {
  cat(sprintf("<psychnet_bootstrap> %s, %d resamples, %.0f%% CI\n",
              x$method, x$n_boot, 100 * x$ci))
  cat(sprintf("  %d edges (%d significant), %d nodes, measures: %s\n",
              nrow(x$edges), sum(x$edges$significant, na.rm = TRUE),
              nrow(x$centrality), paste(x$measures, collapse = ", ")))
  invisible(x)
}
