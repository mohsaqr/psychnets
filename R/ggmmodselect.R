# Stepwise unregularized Gaussian graphical model selection (Epskamp 2020,
# qgraph::ggmModSelect), clean-room base R. The L1 glasso is used only to
# generate candidate edge sets; the reported precision is the UNREGULARIZED
# maximum-likelihood fit constrained to the selected graph, so there is no
# shrinkage of the retained edges. EBIC chooses among graphs, and an optional
# stepwise add/drop search refines the winner. The certificate is the
# constrained-MLE stationarity residual (see [ggm_support_kkt()]).

# Unregularized Gaussian-MRF MLE constrained to a fixed graph `support`
# (symmetric logical). Block-coordinate descent where each column solve is an
# unpenalized least squares restricted to that node's graph-neighbors.
#' @noRd
.ggm_fit_support <- function(S, support, max_outer = 1e4, tol_outer = 1e-8) {
  p <- ncol(S)
  W <- S; diag(W) <- diag(S)
  Beta <- matrix(0, p, p)
  for (outer in seq_len(max_outer)) {
    max_diff <- 0
    for (j in seq_len(p)) {
      idx <- seq_len(p)[-j]
      nbr <- which(support[idx, j])
      beta <- numeric(length(idx))
      if (length(nbr)) {
        W11 <- W[idx, idx, drop = FALSE]
        beta[nbr] <- solve(W11[nbr, nbr, drop = FALSE], S[idx, j][nbr])
      }
      w12 <- as.numeric(W[idx, idx, drop = FALSE] %*% beta)
      d <- max(abs(w12 - W[idx, j]))
      if (d > max_diff) max_diff <- d
      W[idx, j] <- w12; W[j, idx] <- w12
      Beta[j, idx] <- beta
    }
    if (max_diff < tol_outer) break
  }
  Theta <- matrix(0, p, p)
  for (j in seq_len(p)) {
    idx   <- seq_len(p)[-j]
    beta  <- Beta[j, idx]
    denom <- W[j, j] - sum(W[idx, j] * beta)
    tjj   <- if (abs(denom) > 1e-12) 1 / denom else 1e6
    Theta[j, j]   <- tjj
    Theta[idx, j] <- -beta * tjj
  }
  (Theta + t(Theta)) / 2
}

# A small recorder for fits that error during the model search: `try(expr)`
# returns the value, or NULL after recording the error; `count()` and
# `first()` report what was skipped so the caller can surface it.
#' @noRd
.ggm_failure_log <- function() {
  n_failed <- 0L
  first_msg <- NULL
  list(
    try = function(expr) {
      tryCatch(expr, error = function(e) {
        n_failed <<- n_failed + 1L
        if (is.null(first_msg)) first_msg <<- conditionMessage(e)
        NULL
      })
    },
    count = function() n_failed,
    first = function() first_msg)
}

# EBIC of an unregularized constrained fit. The model dimension is the number of
# edges in the SUPPORT (free off-diagonal parameters), not the count of nonzero
# MLE entries: a structurally present edge can estimate to exactly zero on a
# degenerate S, which would otherwise undercount the penalty.
#' @noRd
.support_ebic <- function(Theta, S, n, gamma, p, support) {
  ld <- determinant(Theta, logarithm = TRUE)
  if (ld$sign <= 0) return(Inf)
  loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% Theta)))
  npar   <- sum(support[upper.tri(support)])
  -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)
}

#' Stepwise Gaussian graphical model selection (ggmModSelect)
#'
#' Selects a GGM by extended-BIC model search over edge sets generated from the
#' glasso path, refitting the *unregularized* maximum-likelihood precision on
#' each candidate graph, with an optional stepwise add/drop search. Unlike the
#' graphical lasso, retained edges are not shrunk. Equivalent in purpose to
#' `qgraph::ggmModSelect()`, but pure base R and self-certified via
#' [ggm_support_kkt()].
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix; if given, `n` is required.
#' @param n Sample size (required when `cor_matrix` is supplied): a single
#'   finite number >= 1.
#' @param gamma EBIC hyperparameter. Default `0`, matching
#'   `qgraph::ggmModSelect()`: because the selected graph is refit with the
#'   *unregularized* MLE (no edge shrinkage), the extra EBIC penalty is not
#'   needed, so gamma 0 (plain BIC) is the method's intended setting. (The
#'   regularized GGMs [ebic_glasso()] and [huge_network()] keep gamma 0.5.)
#' @param stepwise If `TRUE` (default), refine the best glasso-path graph by a
#'   greedy single-edge add/drop search.
#' @param nlambda Number of glasso penalties scanned for candidate graphs (a
#'   whole number >= 2).
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest, in
#'   (0, 1).
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed: a single number >= 0, default 0. With a positive value `$kkt` is `NA` for the returned
#'   graph and the pre-threshold residual is stored as `$fit_kkt`.
#' @param cor_method Correlation used when `data` is supplied: `"pearson"`
#'   (default), `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial,
#'   the `qgraph`/`bootnet` default for ordinal items). See [cor_auto()].
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param native Solver switch for generating candidate supports: `TRUE`
#'   (default) uses the pure-R solver; `FALSE` delegates to the `glasso` Fortran
#'   package (in `Suggests`). The reported precision is the unregularized refit
#'   either way. See [ebic_glasso()].
#' @param labels Optional node labels, one unique label per network node.
#' @return A `psychnet` object whose `$weights` is the partial-correlation matrix,
#'   with `$precision`, `$support` (the selected graph), `$gamma`, `$ebic`,
#'   `$cor_matrix`, and `$kkt`. Invalid arguments raise `psychnet_bad_input`.
#'   A candidate or stepwise fit that errors is skipped and reported by one
#'   `psychnet_fit_failed` warning; a selected fit whose certificate exceeds
#'   `1e-6` raises a `psychnet_no_converge` warning; if no candidate can be
#'   fitted at all the call errors (`psychnet_no_converge`).
#' @references Williams, D. R., Rhemtulla, M., Wysocki, A. C., & Rast, P.
#'   (2019). On nonregularized estimation of psychological networks.
#'   *Multivariate Behavioral Research*, 54(5), 719-750.
#'   \doi{10.1080/00273171.2019.1575716}
#'
#'   Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for
#'   Gaussian graphical models. *Advances in Neural Information Processing
#'   Systems*, 23, 604-612.
#' @examples
#' S <- 0.5^abs(outer(1:6, 1:6, "-"))
#' ggm_modselect(cor_matrix = S, n = 250)
#' @export
ggm_modselect <- function(data = NULL, cor_matrix = NULL, n = NULL,
                          gamma = 0, stepwise = TRUE, nlambda = 100L,
                          lambda_min_ratio = 0.01, threshold = 0,
                          cor_method = c("pearson", "spearman", "kendall", "auto"),
                          na_method = c("pairwise", "listwise"),
                          native = TRUE, labels = NULL) {
  na_method <- .psn_match_arg(na_method)
  cor_method <- .psn_match_arg(cor_method)
  engine <- .resolve_native(native, "glasso")
  .check_ebic_args(gamma, nlambda, lambda_min_ratio, threshold)
  .psn_check_flag(stepwise, "stepwise")
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "required")
  S <- inp$S; n <- inp$n; labels <- inp$labels
  p <- ncol(S)
  # Candidate and stepwise fits that error (a singular block on a degenerate
  # S) are skipped so one bad graph cannot abort the search, but every skip is
  # counted and reported below as a `psychnet_fit_failed` warning.
  failures <- .ggm_failure_log()

  # No off-diagonal association: the empty graph is the optimum, skip the search.
  if (max(abs(S[upper.tri(S)])) <= 1e-12) {
    best_support <- matrix(FALSE, p, p)
    best_theta <- .ggm_fit_support(S, best_support)
    best_ebic <- .support_ebic(best_theta, S, n, gamma, p, best_support)
    stepwise <- FALSE
  } else {

  # Candidate supports from the glasso path (deduplicated by edge pattern).
  lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
  seen <- character(0)
  best_ebic <- Inf; best_support <- NULL; best_theta <- NULL
  for (lam in lambda_path) {
    wi <- failures$try(.glasso_solve(S, lam, engine, 1e-4, 1e-4)$wi)
    if (is.null(wi)) next
    sup <- (abs(wi) > 1e-6) & upper.tri(wi)
    sup <- sup | t(sup)
    key <- paste(which(sup[upper.tri(sup)]), collapse = ",")
    if (key %in% seen) next
    seen <- c(seen, key)
    theta <- failures$try(.ggm_fit_support(S, sup))
    if (is.null(theta)) next
    eb <- .support_ebic(theta, S, n, gamma, p, sup)
    if (eb < best_ebic) { best_ebic <- eb; best_support <- sup; best_theta <- theta }
  }
  if (is.null(best_support)) {
    .psn_stop(sprintf(paste0(
      "Model selection failed: none of the candidate graphs could be fitted ",
      "(%d failed fits; first error: %s)."), failures$count(),
      failures$first() %||% "none"),
      class = c("psychnet_no_converge", "psychnet_bad_input"))
  }
  }

  # Greedy stepwise single-edge add/drop refinement.
  if (stepwise) {
    ut <- which(upper.tri(S), arr.ind = TRUE)
    improved <- TRUE; guard <- 0L
    while (improved && guard < 200L) {
      improved <- FALSE; guard <- guard + 1L
      for (e in seq_len(nrow(ut))) {
        i <- ut[e, 1L]; j <- ut[e, 2L]
        cand <- best_support
        cand[i, j] <- cand[j, i] <- !cand[i, j]
        theta <- failures$try(.ggm_fit_support(S, cand))
        if (is.null(theta)) next
        eb <- .support_ebic(theta, S, n, gamma, p, cand)
        if (eb < best_ebic - 1e-10) {
          best_ebic <- eb; best_support <- cand; best_theta <- theta
          improved <- TRUE
        }
      }
    }
  }

  # Defensive idempotent refit: guarantees best_theta == fit(best_support) even
  # if no candidate improved the initial support (same tolerance throughout).
  best_theta <- .ggm_fit_support(S, best_support)
  dimnames(best_theta) <- dimnames(best_support) <- dimnames(S) <- list(labels, labels)
  pcor <- .precision_to_pcor(best_theta)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)
  fit_kkt <- ggm_support_kkt(best_theta, S, best_support)
  if (failures$count() > 0L) {
    .psn_warn(sprintf(paste0(
      "%d candidate graph fit(s) failed and were skipped during model ",
      "selection (first error: %s); the selected graph is the best of the ",
      "fits that succeeded."), failures$count(), failures$first()),
      class = "psychnet_fit_failed")
  }
  if (!is.finite(fit_kkt) || fit_kkt > 1e-6) {
    .psn_warn(sprintf(paste0(
      "The constrained MLE on the selected graph did not reach its optimum ",
      "(certificate %.2e > 1e-6); the returned partial correlations are not ",
      "certified."), fit_kkt), class = "psychnet_no_converge")
  }

  .new_psychnet(
    graph = pcor, labels = labels, method = "ggm",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = best_theta, support = best_support, gamma = gamma,
      cor_matrix = S, ebic = best_ebic,
      kkt = if (threshold > 0) NA_real_ else fit_kkt,
      fit_kkt = fit_kkt,
      certificate_target = if (threshold > 0) "pre_threshold_fit" else
        "returned_network"
    )
  )
}
