# Clean-room graphical lasso (Friedman, Hastie & Tibshirani 2008, Biostatistics)
# and EBIC model selection (Foygel & Drton 2010). Pure base R, no compiled
# dependency. The objective
#     min_{Theta > 0}  -log det Theta + tr(S Theta) + rho * sum_{i != j}|Theta_ij|
# is strictly convex, so its minimiser is unique; glasso_kkt() grades any
# candidate against that unique optimum with no reference solver.
#
# The covariance (W) block-coordinate descent solves, for each column, an
# ordinary lasso by Gauss-Seidel soft-thresholding, then reconstructs the
# precision matrix from the converged W and the lasso coefficients. The outer
# column sweep and inner coordinate descent are Gauss-Seidel: each update
# consumes the latest values, a genuine sequential dependency (the inner dot
# products are vectorised).

# --- soft-threshold ----------------------------------------------------------
.soft <- function(z, g) sign(z) * pmax(abs(z) - g, 0)

# --- single-column lasso (coordinate descent) --------------------------------
# Pathwise coordinate descent with covariance/residual updates (Friedman,
# Hastie & Tibshirani 2010): keep a running r = W11 %*% beta so each coordinate
# reads r[k] in O(1), and apply an O(p) rank-1 update to r only when a
# coordinate actually moves -- in a sparse network most stay zero, so the inner
# sweeps are far cheaper than recomputing a full dot product per coordinate. The
# update is algebraically identical to the naive sweep, so it converges to the
# same unique optimum of the strictly convex objective.
.glasso_lasso_column <- function(W11, s12, beta, rho, max_inner, tol_inner, dgg) {
  pp  <- length(s12)
  r   <- as.numeric(W11 %*% beta)           # running W11 %*% beta
  for (inner in seq_len(max_inner)) {
    max_diff <- 0
    for (k in seq_len(pp)) {
      wkk <- dgg[k]
      if (wkk < 1e-12) {
        new_k <- 0
      } else {
        partial <- s12[k] - (r[k] - wkk * beta[k])
        thr <- abs(partial) - rho           # soft-threshold, inlined (scalar)
        new_k <- if (thr > 0) (if (partial > 0) thr else -thr) / wkk else 0
      }
      delta <- new_k - beta[k]
      if (delta != 0) {
        r <- r + W11[, k] * delta           # rank-1 residual update
        beta[k] <- new_k
        ad <- abs(delta); if (ad > max_diff) max_diff <- ad
      }
    }
    if (max_diff < tol_inner) break
  }
  beta
}

# --- single graphical-lasso fit at fixed penalty -----------------------------
.glasso_fit <- function(S, rho,
                        max_outer = 1e4, tol_outer = 1e-8,
                        max_inner = 1e4, tol_inner = 1e-10,
                        w_init = NULL, beta_init = NULL,
                        penalize_diagonal = FALSE) {
  p <- ncol(S)
  W <- if (is.null(w_init)) S else w_init
  # The working covariance diagonal is diag(S) when the diagonal is unpenalized
  # and diag(S) + rho when it is. `dW` (not diag(S)) is what the column solver
  # needs: it is the actual W11[k, k] every inner update divides by, so the two
  # must be derived from the same quantity or the solver converges to a point
  # that is not the optimum of either objective.
  dW <- diag(S) + if (penalize_diagonal) rho else 0
  diag(W) <- dW
  Beta <- if (is.null(beta_init)) matrix(0, p, p) else beta_init

  # diag(W) is held fixed (only off-diagonals update), so each column's
  # W11[k, k] is just dW[idx] -- precompute it instead of diag(W11).
  for (outer in seq_len(max_outer)) {
    max_diff <- 0
    for (j in seq_len(p)) {
      idx  <- seq_len(p)[-j]
      W11  <- W[idx, idx, drop = FALSE]
      s12  <- S[idx, j]
      beta <- .glasso_lasso_column(W11, s12, Beta[j, idx], rho,
                                   max_inner, tol_inner, dW[idx])
      w12  <- as.numeric(W11 %*% beta)
      d <- max(abs(w12 - W[idx, j]))
      if (d > max_diff) max_diff <- d
      W[idx, j] <- w12
      W[j, idx] <- w12
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
  Theta <- (Theta + t(Theta)) / 2
  dimnames(Theta) <- dimnames(W) <- dimnames(S)
  list(wi = Theta, w = W, beta = Beta)
}

# Validate a caller-supplied penalty grid. Returns NULL (compute the path) or
# the vector unchanged -- the order is meaningful, so it is never sorted.
#' @noRd
.check_lambda_path <- function(lambda_path) {
  if (is.null(lambda_path)) return(NULL)
  if (!is.numeric(lambda_path) || length(lambda_path) < 1L) {
    .psn_stop("`lambda_path` must be a numeric vector of length >= 1.", call. = FALSE)
  }
  if (anyNA(lambda_path) || any(!is.finite(lambda_path)) || any(lambda_path <= 0)) {
    .psn_stop("`lambda_path` values must be finite and strictly positive.",
         call. = FALSE)
  }
  as.numeric(lambda_path)
}

# The EBIC-path arguments shared by ebic_glasso(), huge_network() and
# ggm_modselect(), plus the absolute-weight `threshold`.
#' @noRd
.check_ebic_args <- function(gamma, nlambda, lambda_min_ratio, threshold) {
  .psn_check_number(gamma, "gamma", lower = 0)
  .psn_check_number(nlambda, "nlambda", lower = 2, whole = TRUE)
  .psn_check_number(lambda_min_ratio, "lambda_min_ratio", lower = 0, upper = 1,
                    open = TRUE)
  .check_threshold(threshold)
  invisible(NULL)
}

# `refit` is TRUE / FALSE / "unregularized".
#' @noRd
.check_refit <- function(refit) {
  if (is.character(refit) && length(refit) == 1L &&
      identical(refit, "unregularized")) {
    return("unregularized")
  }
  if (is.logical(refit) && length(refit) == 1L && !is.na(refit)) return(refit)
  .psn_stop('`refit` must be TRUE, FALSE, or "unregularized".', call. = FALSE)
}

# --- log-spaced lambda path --------------------------------------------------
.compute_lambda_path <- function(S, nlambda, lambda_min_ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))
  if (lambda_max <= 0) {
    .psn_stop("All off-diagonal correlations are zero; nothing to regularise.",
         call. = FALSE)
  }
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
          length.out = nlambda))
}

# Empty-graph fit. When no off-diagonal association exists the optimum is the
# diagonal precision diag(1/S_ii) at lambda 0 (FHT 2008 KKT: W = S, every
# inactive residual 0 <= rho). Returns the same shape as .select_ebic().
.empty_glasso <- function(S, n) {
  wi <- diag(1 / diag(S), nrow = ncol(S))
  dimnames(wi) <- dimnames(S)
  loglik <- (n / 2) * (sum(log(1 / diag(S))) - ncol(S))   # tr(S wi) = p
  ebic <- -2 * loglik                                      # npar = 0
  # The empty graph's support is the diagonal, so refit = "unregularized" can
  # certify it with ggm_support_kkt() like any other selected support.
  support <- diag(TRUE, ncol(S)); dimnames(support) <- dimnames(S)
  list(wi = wi, lambda = 0, ebic = ebic, ebic_path = ebic, support = support)
}

# Translate the public `native` switch into the internal solver name and check
# the optional package is present. `native = TRUE` (default) is the pure-R,
# dependency-free, self-certified solver; `native = FALSE` delegates to the
# established external package -- `ext` is "glasso" (Fortran) for the GGMs and
# "glmnet" for the nodewise estimators -- for byte-identical reference parity,
# at the cost of that solver's looser convergence.
#' @noRd
.resolve_native <- function(native, ext) {
  .psn_check_flag(native, "native")
  if (native) return("base")
  if (!requireNamespace(ext, quietly = TRUE)) {
    .psn_stop(sprintf(paste0("native = FALSE needs the '%s' package; install it ",
                        "or keep native = TRUE (the default)."), ext),
         call. = FALSE)
  }
  ext
}

# One fixed-penalty graphical-lasso solve, dispatched by engine. Returns
# list(wi, w, beta) (beta is NULL for the glasso engine). The base path is
# unchanged, so engine = "base" is byte-identical to the previous behaviour.
#' @noRd
.glasso_solve <- function(S, rho, engine, tol_outer, tol_inner, warm = NULL,
                          penalize_diagonal = FALSE) {
  if (engine == "base") {
    return(.glasso_fit(S, rho, tol_outer = tol_outer, tol_inner = tol_inner,
                       w_init = warm$w, beta_init = warm$beta,
                       penalize_diagonal = penalize_diagonal))
  }
  g <- glasso::glasso(S, rho, penalize.diagonal = penalize_diagonal,
                      thr = tol_outer)
  wi <- (g$wi + t(g$wi)) / 2
  dimnames(wi) <- dimnames(S)
  list(wi = wi, w = g$w, beta = NULL)
}

# --- two-tier EBIC selection: fast scan, then tight refit at the winner ------
# `refit` controls what is returned at the selected penalty:
#   TRUE            re-solve to machine precision (the certified optimum; default)
#   FALSE           the scan fit exactly as computed at scan_tol
#   "unregularized" the unpenalized constrained MLE on the selected support
# A lambda whose fit is non-PD is skipped rather than raising, because
# resampling callers hit that routinely on near-singular draws. The skip is
# not silent in effect: the lambda scores EBIC = Inf, so it is visible in the
# returned `ebic_path`, and if every lambda fails the call errors. (The
# tryCatch below therefore deliberately catches any solver error.)
.select_ebic <- function(S, lambda_path, n, gamma, engine = "base",
                         scan_tol = 1e-4, refit_tol = 1e-8, refit = TRUE,
                         penalize_diagonal = FALSE) {
  p <- ncol(S)
  ebic_vals <- numeric(length(lambda_path))
  w_prev <- NULL; beta_prev <- NULL
  best_idx <- 1L; best_ebic <- Inf; have_best <- FALSE
  best_scan <- NULL

  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(
      .glasso_solve(S, lambda_path[i], engine, scan_tol, scan_tol,
                    warm = list(w = w_prev, beta = beta_prev),
                    penalize_diagonal = penalize_diagonal),
      error = function(e) NULL
    )
    if (is.null(fit)) { ebic_vals[i] <- Inf; next }
    w_prev <- fit$w; beta_prev <- fit$beta

    ld <- determinant(fit$wi, logarithm = TRUE)
    if (ld$sign <= 0) { ebic_vals[i] <- Inf; next }
    loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% fit$wi)))
    npar   <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic_vals[i] <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)

    # Strict `<` scanning in path order: ties go to the first (largest-lambda)
    # minimiser, i.e. the sparsest graph wins.
    if (ebic_vals[i] < best_ebic) {
      best_ebic <- ebic_vals[i]; best_idx <- i; have_best <- TRUE
      best_scan <- fit$wi
    }
  }
  if (!have_best) .psn_stop("All glasso fits failed; check the input data.",
                       call. = FALSE)

  best_lambda <- lambda_path[best_idx]
  support <- NULL
  confirmed_wi <- NULL
  if (engine == "base" && !isFALSE(refit)) {
    # Confirm the selection at the refit tolerance. The scan solves each lambda
    # only to scan_tol, which can leave an entry that is zero at the optimum
    # just above the 1e-10 edge-count threshold; one miscounted edge raises that
    # lambda's EBIC by exactly log(n) + 4 gamma log(p) and can hide the true
    # minimum (the NCT parity audit found this in 4 of 29 divergent fits). The
    # true minimum therefore lies within that margin of the scan minimum, so
    # every lambda within two edge penalties is refit tightly and rescored, and
    # the lowest exact EBIC wins (ties to the larger lambda, as in the scan).
    pen <- log(n) + 4 * gamma * log(p)
    cand <- which(is.finite(ebic_vals) & ebic_vals <= best_ebic + 2 * pen)
    rescore <- function(i) {
      wi <- tryCatch(
        .glasso_solve(S, lambda_path[i], engine, refit_tol, refit_tol * 1e-2,
                      penalize_diagonal = penalize_diagonal)$wi,
        error = function(e) NULL)
      if (is.null(wi)) return(list(wi = NULL, ebic = Inf))
      wi[abs(wi) <= 1e-10 & row(wi) != col(wi)] <- 0
      ld <- determinant(wi, logarithm = TRUE)
      if (ld$sign <= 0) return(list(wi = NULL, ebic = Inf))
      ll <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% wi)))
      k <- sum(wi[upper.tri(wi)] != 0)
      list(wi = wi, ebic = -2 * ll + k * log(n) + 4 * k * gamma * log(p))
    }
    confirmed <- vector("list", length(lambda_path))
    repeat {
      rescored <- lapply(cand, rescore)
      confirmed[cand] <- rescored
      ebic_vals[cand] <- vapply(rescored, function(r) r$ebic, numeric(1))
      if (!any(is.finite(ebic_vals)))
        .psn_stop("All glasso refits failed; check the input data.", call. = FALSE)
      best_idx <- which.min(ebic_vals)
      # A correction can move the minimum to an unconfirmed scan row. Confirm
      # that row too, so the returned path and selected fit agree.
      if (!is.null(confirmed[[best_idx]])) break
      cand <- best_idx
    }
    best_lambda <- lambda_path[best_idx]
    best_ebic <- ebic_vals[best_idx]
    confirmed_wi <- confirmed[[best_idx]]$wi
    best_scan <- confirmed_wi
  }
  if (identical(refit, "unregularized")) {
    # Zero-constrained unregularized refit: keep the EBIC-selected support and
    # re-solve the graph-restricted MLE, certified by ggm_support_kkt() instead
    # of glasso_kkt() (the returned network is no longer a penalized optimum).
    support <- abs(best_scan) > 1e-10
    diag(support) <- TRUE
    best_wi <- .ggm_fit_support(S, support)
  } else if (!is.null(confirmed_wi)) {
    best_wi <- confirmed_wi
  } else if (isFALSE(refit)) {
    best_wi <- best_scan
  } else {
    # base refits to the certified optimum; glasso stays at its own tolerance so
    # the result is exactly what the glasso package returns.
    refit_outer <- if (engine == "glasso") scan_tol else refit_tol
    refit_inner <- if (engine == "glasso") scan_tol else refit_tol * 1e-2
    best_wi <- .glasso_solve(S, best_lambda, engine, refit_outer, refit_inner,
                             penalize_diagonal = penalize_diagonal)$wi
  }
  # Return exactly the support the EBIC counted: an off-diagonal entry it
  # treated as absent (|.| <= 1e-10, as in `npar`) becomes an exact zero. At
  # lambda = max|S_ij| (the empty graph) the soft-threshold otherwise leaves
  # ~1e-17 residues, which read as edges downstream -- a permutation test
  # compares statistics with `>=`, so a residue against an exact zero is a
  # false difference.
  absent <- abs(best_wi) <= 1e-10 & row(best_wi) != col(best_wi)
  best_wi[absent] <- 0
  dimnames(best_wi) <- dimnames(S)
  list(wi = best_wi, lambda = best_lambda, ebic = best_ebic,
       ebic_path = ebic_vals, support = support)
}

# Shared input contract of the two certificates: a finite square `theta` and a
# finite `cor_matrix` of the same size.
#' @noRd
.check_kkt_matrices <- function(theta, cor_matrix) {
  if (!is.numeric(theta) || nrow(theta) != ncol(theta) || !all(is.finite(theta)))
    .psn_stop("`theta` must be a finite numeric square matrix.")
  if (!is.numeric(cor_matrix) || !identical(dim(cor_matrix), dim(theta)) ||
      !all(is.finite(cor_matrix)))
    .psn_stop("`cor_matrix` must be a finite numeric matrix the size of `theta`.")
  invisible(NULL)
}

# --- precision -> partial correlation ----------------------------------------
.precision_to_pcor <- function(wi) {
  x <- -stats::cov2cor(wi)
  diag(x) <- 0
  # `+ 0` turns the -0 that negating an exact-zero precision entry produces
  # into 0, so an absent edge prints as 0.000, not -0.000 (values unchanged).
  (x + t(x)) / 2 + 0
}

#' Check that a precision matrix solves the graphical lasso
#'
#' @description
#' `glasso_kkt()` checks whether a precision matrix solves the graphical lasso
#' for a given correlation or covariance matrix and penalty. The graphical
#' lasso estimates a Gaussian graphical model as the precision (inverse
#' covariance) matrix that maximizes the Gaussian likelihood minus a penalty
#' on the absolute partial associations (Friedman et al., 2008), and a matrix
#' is its solution if and only if it satisfies the optimality
#' (Karush-Kuhn-Tucker) conditions of the problem.
#'
#' The objective is strictly convex, so it has exactly one optimum, and
#' checking the conditions shows whether a network is the graphical lasso
#' solution without solving the problem a second time. In words, the
#' covariance implied by the precision matrix must equal the observed
#' covariance on the diagonal, differ from it by exactly the penalty for every
#' pair joined by an edge, and by at most the penalty for every pair that is
#' not.
#'
#' `glasso_kkt()` returns the largest violation of the conditions; a value
#' near zero shows that the matrix is the optimum, up to numerical precision.
#'
#' Networks from [ebic_glasso()] and [huge_network()] store this value, and
#' [certificate()] reports it. Unpenalized networks restricted to a graph are
#' checked by [ggm_support_kkt()], and nodewise Ising and mixed graphical
#' models by [glm_lasso_kkt()].
#'
#' @details
#' The objective, with the penalty on the off-diagonal entries, is
#' \deqn{\min_{\Theta \succ 0} -\log\det\Theta + \mathrm{tr}(S\Theta) +
#'       \rho \sum_{i \neq j} |\Theta_{ij}|.}
#' With \eqn{W = \Theta^{-1}}, the optimality conditions are
#' \eqn{W_{ii} = S_{ii}};
#' \eqn{W_{ij} - S_{ij} = \rho\,\mathrm{sign}(\Theta_{ij})} where
#' \eqn{\Theta_{ij} \neq 0}; and \eqn{|W_{ij} - S_{ij}| \le \rho} otherwise.
#' The function returns the largest absolute deviation from these conditions.
#' Because the objective is strictly convex, a precision matrix with zero
#' violation is its unique global optimum.
#'
#' When the penalty also applies to the diagonal (`penalize_diagonal = TRUE`,
#' the objective adds \eqn{\rho \sum_i |\Theta_{ii}|}), the diagonal condition
#' becomes \eqn{W_{ii} = S_{ii} + \rho}. Checking a diagonally penalized fit
#' with the default would report a violation of exactly \eqn{\rho}, so this
#' argument must match the fit being checked.
#'
#' The objective is defined only for symmetric positive-definite matrices, so
#' a `theta` that is asymmetric (by more than `1e-8`) or not positive definite
#' is not a solution and returns `Inf`.
#'
#' @param theta The precision matrix to check: a finite, square numeric
#'   matrix.
#' @param cor_matrix The correlation or covariance matrix the model was fitted
#'   to, of the same size as `theta`.
#' @param rho The penalty: a single finite number >= 0.
#' @param active_tol Absolute value above which an off-diagonal entry of
#'   `theta` counts as an edge. Default `1e-8`.
#' @param penalize_diagonal Whether the fit being checked also penalized the
#'   diagonal. Default `FALSE`, the default of [ebic_glasso()].
#' @return A single number, the largest absolute violation of the optimality
#'   conditions: 0 at the exact optimum, `Inf` for an asymmetric or
#'   non-positive-definite `theta`.
#' @section Conditions:
#'   `psychnet_bad_input` when `theta` is not a finite square numeric matrix,
#'   `cor_matrix` is not a finite matrix of the same size, `rho` or
#'   `active_tol` is not a single finite number >= 0, or `penalize_diagonal`
#'   is not `TRUE` or `FALSE`.
#' @references
#'   Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#'   covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
#'   432-441. \doi{10.1093/biostatistics/kxm045}
#' @seealso [certificate()], which reports the stored value of a fitted
#'   network; [ebic_glasso()] and [huge_network()]; [ggm_support_kkt()] and
#'   [glm_lasso_kkt()] for the other optimality checks.
#' @examples
#' S <- 0.5^abs(outer(1:5, 1:5, "-"))
#' # At rho = 0 the optimum is the inverse of S, so the violation is ~0 ...
#' glasso_kkt(solve(S), S, rho = 0)
#' # ... and a perturbed precision matrix violates the conditions.
#' glasso_kkt(solve(S) + diag(0.1, 5), S, rho = 0)
#' # A fitted network reports the same value through certificate().
#' certificate(ebic_glasso(cor_matrix = S, n = 200))
#' @export
glasso_kkt <- function(theta, cor_matrix, rho, active_tol = 1e-8,
                       penalize_diagonal = FALSE) {
  theta <- as.matrix(theta); cor_matrix <- as.matrix(cor_matrix)
  .check_kkt_matrices(theta, cor_matrix)
  .psn_check_number(rho, "rho", lower = 0)
  .psn_check_number(active_tol, "active_tol", lower = 0)
  .psn_check_flag(penalize_diagonal, "penalize_diagonal")
  # Feasibility is part of the optimality conditions: the objective is defined
  # only on symmetric Theta > 0, so an asymmetric, indefinite, or numerically
  # non-invertible matrix is infeasible, not optimal.
  if (max(abs(theta - t(theta))) > 1e-8) return(Inf)
  if (min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    return(Inf)
  W <- tryCatch(solve(theta), error = function(e) NULL)
  if (is.null(W)) return(Inf)
  # W_ii = S_ii, or S_ii + rho when the diagonal carries the penalty too.
  diag_target <- diag(cor_matrix) + if (penalize_diagonal) rho else 0
  diag_v <- max(abs(diag(W) - diag_target))
  off <- upper.tri(theta)
  r  <- (W - cor_matrix)[off]
  th <- theta[off]
  active <- abs(th) > active_tol
  v_a <- if (any(active))  max(abs(r[active] - rho * sign(th[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(r[!active]) - rho, 0)) else 0
  max(diag_v, v_a, v_i)
}

#' Check that a precision matrix is the maximum-likelihood network on a graph
#'
#' @description
#' `ggm_support_kkt()` checks whether a precision matrix is the
#' maximum-likelihood Gaussian graphical model on a given graph, for the
#' correlation or covariance matrix it was fitted to. Some Gaussian graphical
#' models are estimated without a penalty: a graph is chosen first, and the
#' precision matrix is then the maximum-likelihood estimate under the
#' constraint that every pair of variables not joined by an edge has a
#' partial correlation of zero.
#'
#' This constrained estimate is unique and is characterized by two
#' conditions: the covariance implied by the model reproduces the observed
#' covariance for every pair joined by an edge and for every variance, and the
#' precision matrix is zero for every pair that is not joined. A matrix
#' meeting both conditions is the maximum-likelihood network on that graph.
#'
#' `ggm_support_kkt()` returns the largest violation of these conditions; a
#' value near zero shows that the matrix is the constrained optimum, up to
#' numerical precision.
#'
#' Networks from [ggm_modselect()], [logo_network()] and
#' `ebic_glasso(refit = "unregularized")` store this value, and
#' [certificate()] reports it. Penalized graphical lasso networks are checked
#' by [glasso_kkt()].
#'
#' @details
#' For a Gaussian graphical model on a graph \eqn{G}, with
#' \eqn{W = \Theta^{-1}}, the maximum-likelihood conditions are
#' \eqn{W_{ij} = S_{ij}} for every \eqn{(i,j)} in \eqn{G} and on the diagonal,
#' and \eqn{\Theta_{ij} = 0} for every \eqn{(i,j)} not in \eqn{G}. The function
#' returns the largest absolute deviation from these conditions; a `theta`
#' that is asymmetric (by more than `1e-8`) or not positive definite returns
#' `Inf`.
#'
#' @param theta The precision matrix to check: a finite, square numeric
#'   matrix.
#' @param cor_matrix The correlation or covariance matrix the model was fitted
#'   to, of the same size as `theta`.
#' @param support A logical matrix of the same size, `TRUE` where an edge is
#'   allowed, with no missing values. It is made symmetric, so an upper- or
#'   lower-triangular matrix describes the same graph; the diagonal is
#'   ignored.
#' @param active_tol A single finite number >= 0 (default `1e-8`). It is
#'   checked but does not affect the result: every entry of `theta` outside
#'   the graph counts as a violation of its absolute size.
#' @return A single number, the largest absolute violation of the conditions:
#'   0 at the exact constrained optimum, `Inf` for an asymmetric or
#'   non-positive-definite `theta`.
#' @section Conditions:
#'   `psychnet_bad_input` when `theta` is not a finite square numeric matrix,
#'   `cor_matrix` is not a finite matrix of the same size, `support` is not a
#'   logical matrix of that size without missing values, or `active_tol` is
#'   not a single finite number >= 0.
#' @seealso [certificate()], which reports the stored value of a fitted
#'   network; [ggm_modselect()] and [logo_network()]; [glasso_kkt()] for
#'   penalized networks.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' # An AR(1) correlation has a tridiagonal inverse: the chain graph is exact.
#' chain <- abs(row(S) - col(S)) <= 1
#' ggm_support_kkt(solve(S), S, support = chain)
#' # The same precision violates an empty-graph (diagonal) constraint.
#' ggm_support_kkt(solve(S), S, support = diag(6) > 0)
#' # A fitted network reports its value through certificate().
#' certificate(ggm_modselect(cor_matrix = S, n = 250))
#' @export
ggm_support_kkt <- function(theta, cor_matrix, support, active_tol = 1e-8) {
  theta <- as.matrix(theta); cor_matrix <- as.matrix(cor_matrix)
  .check_kkt_matrices(theta, cor_matrix)
  .psn_check_number(active_tol, "active_tol", lower = 0)
  support <- as.matrix(support)
  if (!is.logical(support) || !identical(dim(support), dim(theta)) ||
      anyNA(support)) {
    .psn_stop(sprintf(paste0(
      "`support` must be a logical %d x %d matrix (TRUE where an edge is ",
      "allowed) with no missing values."), nrow(theta), ncol(theta)))
  }
  # An edge is a property of the unordered pair: a triangular support (only
  # upper or only lower filled in) means the symmetric graph.
  support <- support | t(support)
  if (max(abs(theta - t(theta))) > 1e-8) return(Inf)
  if (min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    return(Inf)
  W <- tryCatch(solve(theta), error = function(e) NULL)
  if (is.null(W)) return(Inf)
  diag_v <- max(abs(diag(W) - diag(cor_matrix)))
  off <- upper.tri(theta)
  sup <- support[off]
  r   <- (W - cor_matrix)[off]
  th  <- theta[off]
  v_match <- if (any(sup))  max(abs(r[sup])) else 0       # W = S on the graph
  v_zero  <- if (any(!sup)) max(abs(th[!sup])) else 0     # Theta = 0 off-graph
  max(diag_v, v_match, v_zero)
}

#' Estimate a Gaussian graphical model with the EBIC graphical lasso
#'
#' @description
#' `ebic_glasso()` estimates a Gaussian graphical model: a network in which
#' variables are nodes joined by partial correlations, and a missing edge means
#' two variables are independent given all the others. It fits the graphical
#' lasso, which shrinks weak partial correlations to zero, and chooses the
#' penalty by the extended Bayesian information criterion (Epskamp & Fried,
#' 2018).
#'
#' An edge is the association between two variables after controlling for all
#' other variables in the network. The penalty sets weak partial correlations
#' to exactly zero, so that only the associations supported by the data are
#' retained, and the EBIC balances the fit of the model against the number of
#' edges. `ebic_glasso()` fits the graphical lasso over a sequence of
#' penalties, selects the penalty with the lowest EBIC and returns the
#' partial-correlation network at that penalty. It takes raw data, or a
#' correlation matrix together with its sample size. `gamma` sets how strongly
#' the EBIC favours sparse networks, `refit = "unregularized"` keeps the
#' selected edges but estimates them without shrinkage, and `threshold`
#' removes weak edges after estimation.
#'
#' The returned network is the optimum of the penalized likelihood at the
#' selected penalty, which [certificate()] checks.
#'
#' This regularized partial-correlation network is the standard estimator for
#' continuous and ordinal data in psychological network analysis. For data far
#' from normality, [huge_network()] fits the same model after a
#' rank-based transformation. [ggm_modselect()] selects an unregularized
#' network by stepwise search, [pcor_network()] keeps every partial
#' correlation, and [ising_fit()] and [mgm_fit()] estimate networks of binary
#' and mixed variables. The model is the one fitted by `qgraph::EBICglasso()`.
#'
#' @details
#' For a precision matrix \eqn{\Theta} with \eqn{E} nonzero off-diagonal
#' entries, fitted to a correlation matrix \eqn{S} of \eqn{p} variables and
#' \eqn{n} observations, the criterion is
#' \deqn{\mathrm{EBIC} = -2L + E \log n + 4 \gamma E \log p, \quad
#'       L = \frac{n}{2}\left(\log\det\Theta - \mathrm{tr}(S\Theta)\right).}
#' The penalties form a log-spaced sequence of `nlambda` values from the
#' largest absolute off-diagonal correlation, at which the network is empty,
#' down to `lambda_min_ratio` times that value. When two penalties give the
#' same EBIC, the one earlier in the path is selected, which on the default
#' path is the larger penalty and the sparser network. Unless `lambda_path` is
#' given, a correlation matrix with no off-diagonal association returns the
#' empty network with penalty 0.
#'
#' With `na_method = "pairwise"` and incomplete data, the sample size entering
#' the EBIC is the median number of observations available for a pair of
#' variables.
#'
#' The certificate is computed by [glasso_kkt()]. With
#' `refit = "unregularized"` the network is no longer a penalized optimum and
#' is checked by [ggm_support_kkt()] instead. With `refit = FALSE` or
#' `native = FALSE` the network is solved to a looser tolerance, and its
#' certificate can exceed the default tolerance of [certificate()]. After a
#' positive `threshold` the certificate is missing; `print()` still shows the
#' value for the network before thresholding.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#'   Not needed when `cor_matrix` is given.
#' @param cor_matrix A correlation matrix to use instead of `data`; a
#'   covariance matrix is rescaled to a correlation matrix. It must be
#'   symmetric and positive semi-definite. `n` is then required.
#' @param n Sample size, required with `cor_matrix`: a single finite number
#'   >= 1. Ignored when `data` is given.
#' @param gamma EBIC hyperparameter, a number >= 0; larger values favour
#'   sparser networks and 0 gives the ordinary BIC. Default `0.5`.
#' @param nlambda Number of penalties scanned, a whole number >= 2. Default
#'   `100`.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest, in
#'   (0, 1). Default `0.01`.
#' @param threshold Absolute partial correlation below which an edge is set to
#'   zero after estimation: a single number >= 0. Default `0`. A logical value
#'   (the `threshold = TRUE` rule of `qgraph`) is refused. A positive value
#'   leaves the certificate missing.
#' @param cor_method Correlation computed from `data`: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric and polyserial
#'   correlations for ordinal items; see [cor_auto()]).
#' @param na_method Handling of missing values in `data`: `"pairwise"`
#'   (default) uses pairwise-complete correlations adjusted to the nearest
#'   positive-definite matrix; `"listwise"` drops incomplete rows. The two
#'   agree for complete data.
#' @param native `TRUE` (default) uses the package's own solver. `FALSE`
#'   solves each penalty with the `glasso` package, which must be installed,
#'   to reproduce results of the `glasso` and `qgraph` packages.
#' @param lambda_path Optional numeric vector of positive penalties to scan,
#'   in the given order, replacing `nlambda` and `lambda_min_ratio`. Useful
#'   when resampling: the path of the full data, reused for every resample,
#'   keeps the selected penalties comparable. Default `NULL`.
#' @param penalize_diagonal Whether the penalty also applies to the diagonal
#'   of the precision matrix. Default `FALSE`.
#' @param refit What to return at the selected penalty: `TRUE` (default) the
#'   optimum at that penalty; `FALSE` the fit obtained while scanning the
#'   penalties, less precise and kept to reproduce earlier results;
#'   `"unregularized"` the unpenalized maximum-likelihood network restricted
#'   to the selected edges, without shrinkage.
#' @param labels Optional node labels, one unique label per node of the
#'   network (the columns left after non-numeric and constant columns are
#'   dropped).
#' @return A fitted network of class `psychnet` ([psychnet-class]) whose edge
#'   weights are partial correlations. `as.data.frame()` returns one row per
#'   edge with columns `from`, `to` and `weight`; `print()` shows the selected
#'   penalty and the certificate; [certificate()] reports the check. The
#'   object also stores the precision matrix (`precision`), the selected
#'   penalty (`lambda`), `gamma`, the correlation matrix used (`cor_matrix`),
#'   the EBIC of the selected network (`ebic`), the penalties scanned and
#'   their EBIC (`lambda_path`, `ebic_path`) and, with
#'   `refit = "unregularized"`, the selected graph (`support`). Nodes keep the
#'   column order of the input.
#' @section Conditions:
#'   `psychnet_bad_input` when an argument is out of range, `threshold` is
#'   logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
#'   matrix with a positive diagonal, `n` is missing with `cor_matrix`,
#'   `data` has fewer than two usable variables, a pair of variables is never
#'   observed together, `labels` do not match the nodes, `native = FALSE`
#'   without the `glasso` package, or every penalty fails to fit;
#'   `psychnet_constant_variable` (warning) when a constant or all-missing
#'   column is dropped. A `psychnet_dropped_columns` message names
#'   non-numeric columns that were dropped.
#' @seealso [certificate()] and [glasso_kkt()] for the check;
#'   [huge_network()], [ggm_modselect()] and [pcor_network()] for related
#'   Gaussian networks; [psychnet()], which routes `method = "glasso"` here;
#'   [psychnet-class] for the result.
#' @references
#'   Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#'   covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
#'   432-441. \doi{10.1093/biostatistics/kxm045}
#'
#'   Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for
#'   Gaussian graphical models. *Advances in Neural Information Processing
#'   Systems*, 23, 604-612.
#'
#'   Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
#'   correlation networks. *Psychological Methods*, 23(4), 617-634.
#'   \doi{10.1037/met0000167}
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 250)
#' fit
#' as.data.frame(fit)
#' certificate(fit)
#' @export
ebic_glasso <- function(data = NULL, cor_matrix = NULL, n = NULL,
                        gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,
                        threshold = 0,
                        cor_method = c("pearson", "spearman", "kendall", "auto"),
                        na_method = c("pairwise", "listwise"),
                        native = TRUE, lambda_path = NULL,
                        penalize_diagonal = FALSE, refit = TRUE,
                        labels = NULL) {
  na_method <- .psn_match_arg(na_method)
  cor_method <- .psn_match_arg(cor_method)
  engine <- .resolve_native(native, "glasso")
  .check_ebic_args(gamma, nlambda, lambda_min_ratio, threshold)
  .psn_check_flag(penalize_diagonal, "penalize_diagonal")
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "required")
  S <- inp$S; n <- inp$n; labels <- inp$labels
  refit <- .check_refit(refit)
  lambda_path <- .check_lambda_path(lambda_path)

  dimnames(S) <- list(labels, labels)
  # A user-supplied path is honoured even on a zero-association matrix: a
  # resampling caller scanning a fixed grid must get the same object shape for
  # every draw, including a degenerate one.
  if (is.null(lambda_path) && max(abs(S[upper.tri(S)])) <= 1e-12) {
    sel <- .empty_glasso(S, n); lambda_path <- 0
  } else {
    if (is.null(lambda_path)) {
      lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
    }
    sel <- .select_ebic(S, lambda_path, n, gamma, engine = engine,
                        refit = refit, penalize_diagonal = penalize_diagonal)
  }

  pcor <- .precision_to_pcor(sel$wi)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  # The certificate must match the objective actually solved: an unregularized
  # support refit is a constrained MLE, not a penalized optimum. `$support` is
  # stored only in that case; the default object stays as lean as before.
  unreg <- identical(refit, "unregularized")
  kkt <- if (unreg) {
    ggm_support_kkt(sel$wi, S, sel$support)
  } else {
    glasso_kkt(sel$wi, S, sel$lambda, penalize_diagonal = penalize_diagonal)
  }
  returned_kkt <- if (threshold > 0) NA_real_ else kkt

  .new_psychnet(
    graph = pcor, labels = labels, method = "glasso",
    directed = FALSE, n_obs = n,
    extra = c(
      list(precision = sel$wi, lambda = sel$lambda, gamma = gamma,
           cor_matrix = S, ebic = sel$ebic, ebic_path = sel$ebic_path,
           lambda_path = lambda_path, na_method = na_method, native = native,
           penalize_diagonal = penalize_diagonal, refit = refit),
      if (unreg) list(support = sel$support),
      list(kkt = returned_kkt,
           fit_kkt = kkt,
           certificate_target = if (threshold > 0) "pre_threshold_fit" else
             "returned_network")
    )
  )
}
