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
.glasso_lasso_column <- function(W11, s12, beta, rho, max_inner, tol_inner) {
  pp <- length(s12)
  for (inner in seq_len(max_inner)) {
    max_diff <- 0
    for (k in seq_len(pp)) {
      partial <- s12[k] - (sum(W11[k, ] * beta) - W11[k, k] * beta[k])
      wkk <- W11[k, k]
      new_k <- if (wkk < 1e-12) 0 else .soft(partial, rho) / wkk
      d <- abs(new_k - beta[k])
      if (d > max_diff) max_diff <- d
      beta[k] <- new_k
    }
    if (max_diff < tol_inner) break
  }
  beta
}

# --- single graphical-lasso fit at fixed penalty -----------------------------
.glasso_fit <- function(S, rho,
                        max_outer = 1e4, tol_outer = 1e-8,
                        max_inner = 1e4, tol_inner = 1e-10,
                        w_init = NULL, beta_init = NULL) {
  p <- ncol(S)
  W <- if (is.null(w_init)) S else w_init
  diag(W) <- diag(S)                                   # penalize.diagonal = FALSE
  Beta <- if (is.null(beta_init)) matrix(0, p, p) else beta_init

  for (outer in seq_len(max_outer)) {
    max_diff <- 0
    for (j in seq_len(p)) {
      idx  <- seq_len(p)[-j]
      W11  <- W[idx, idx, drop = FALSE]
      s12  <- S[idx, j]
      beta <- .glasso_lasso_column(W11, s12, Beta[j, idx], rho,
                                   max_inner, tol_inner)
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

# --- log-spaced lambda path --------------------------------------------------
.compute_lambda_path <- function(S, nlambda, lambda_min_ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))
  if (lambda_max <= 0) {
    stop("All off-diagonal correlations are zero; nothing to regularise.",
         call. = FALSE)
  }
  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
          length.out = nlambda))
}

# --- two-tier EBIC selection: fast scan, then tight refit at the winner ------
.select_ebic <- function(S, lambda_path, n, gamma,
                         scan_tol = 1e-4, refit_tol = 1e-8) {
  p <- ncol(S)
  ebic_vals <- numeric(length(lambda_path))
  w_prev <- NULL; beta_prev <- NULL
  best_idx <- 1L; best_ebic <- Inf; have_best <- FALSE

  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(
      .glasso_fit(S, lambda_path[i], tol_outer = scan_tol, tol_inner = scan_tol,
                  w_init = w_prev, beta_init = beta_prev),
      error = function(e) NULL
    )
    if (is.null(fit)) { ebic_vals[i] <- Inf; next }
    w_prev <- fit$w; beta_prev <- fit$beta

    ld <- determinant(fit$wi, logarithm = TRUE)
    if (ld$sign <= 0) { ebic_vals[i] <- Inf; next }
    loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% fit$wi)))
    npar   <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic_vals[i] <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)

    if (ebic_vals[i] < best_ebic) {
      best_ebic <- ebic_vals[i]; best_idx <- i; have_best <- TRUE
    }
  }
  if (!have_best) stop("All glasso fits failed; check the input data.",
                       call. = FALSE)

  best_wi <- .glasso_fit(S, lambda_path[best_idx],
                         tol_outer = refit_tol, tol_inner = refit_tol * 1e-2)$wi
  dimnames(best_wi) <- dimnames(S)
  list(wi = best_wi, lambda = lambda_path[best_idx], ebic = best_ebic,
       ebic_path = ebic_vals)
}

# --- precision -> partial correlation ----------------------------------------
.precision_to_pcor <- function(wi) {
  x <- -stats::cov2cor(wi)
  diag(x) <- 0
  (x + t(x)) / 2
}

#' Graphical-lasso stationarity (KKT) residual
#'
#' A dependency-free correctness certificate for a fitted Gaussian graphical
#' model. For the convex objective
#' \deqn{\min_{\Theta \succ 0} -\log\det\Theta + \mathrm{tr}(S\Theta) +
#'       \rho \sum_{i \neq j} |\Theta_{ij}|}
#' (off-diagonal penalty), let \eqn{W = \Theta^{-1}}. The subgradient
#' optimality conditions are \eqn{W_{ii} = S_{ii}};
#' \eqn{W_{ij} - S_{ij} = \rho\,\mathrm{sign}(\Theta_{ij})} where
#' \eqn{\Theta_{ij} \neq 0}; and \eqn{|W_{ij} - S_{ij}| \le \rho} otherwise. By
#' strict convexity, a precision matrix with zero violation is the unique global
#' optimum, so a near-zero return certifies correctness independently of any
#' reference solver.
#'
#' @param theta Precision matrix to test.
#' @param cor_matrix Correlation / covariance the model was fit to.
#' @param rho Scalar penalty.
#' @param active_tol Magnitude above which an off-diagonal entry is "active".
#' @return Maximum absolute stationarity violation (scalar); 0 = exact optimum.
#' @examples
#' S <- 0.5^abs(outer(1:5, 1:5, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 200)
#' glasso_kkt(fit$precision, S, fit$lambda)
#' @export
glasso_kkt <- function(theta, cor_matrix, rho, active_tol = 1e-8) {
  # Positive-definiteness is part of the optimality conditions: the objective is
  # defined only on Theta > 0, so an indefinite matrix is infeasible, not optimal.
  if (min(eigen((theta + t(theta)) / 2, symmetric = TRUE,
                only.values = TRUE)$values) <= 0) return(Inf)
  W <- solve(theta)
  diag_v <- max(abs(diag(W) - diag(cor_matrix)))
  off <- upper.tri(theta)
  r  <- (W - cor_matrix)[off]
  th <- theta[off]
  active <- abs(th) > active_tol
  v_a <- if (any(active))  max(abs(r[active] - rho * sign(th[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(r[!active]) - rho, 0)) else 0
  max(diag_v, v_a, v_i)
}

#' Constrained Gaussian-MRF (graph-restricted MLE) stationarity residual
#'
#' Certificate for an *unregularized* Gaussian graphical model whose precision
#' is constrained to a fixed graph (the estimator behind [ggm_modselect()] and
#' [logo_network()]). The maximum-likelihood / maximum-entropy conditions for a
#' Gaussian Markov random field on a graph \eqn{G} are exact:
#' \eqn{W_{ij} = S_{ij}} for every \eqn{(i,j)} on the graph and on the diagonal
#' (\eqn{W = \Theta^{-1}}), and \eqn{\Theta_{ij} = 0} for every \eqn{(i,j)} not
#' on the graph. A near-zero return certifies the constrained optimum with no
#' reference solver.
#'
#' @param theta Precision matrix to test.
#' @param cor_matrix Correlation / covariance the model was fit to.
#' @param support Logical p x p matrix; `TRUE` where an edge is allowed.
#' @param active_tol Magnitude above which an off-support entry counts as a
#'   nonzero violation.
#' @return Maximum absolute stationarity violation (scalar); 0 = exact optimum.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ggm_modselect(cor_matrix = S, n = 250)
#' ggm_support_kkt(fit$precision, S, fit$support)
#' @export
ggm_support_kkt <- function(theta, cor_matrix, support, active_tol = 1e-8) {
  if (min(eigen((theta + t(theta)) / 2, symmetric = TRUE,
                only.values = TRUE)$values) <= 0) return(Inf)
  W <- solve(theta)
  diag_v <- max(abs(diag(W) - diag(cor_matrix)))
  off <- upper.tri(theta)
  sup <- support[off]
  r   <- (W - cor_matrix)[off]
  th  <- theta[off]
  v_match <- if (any(sup))  max(abs(r[sup])) else 0       # W = S on the graph
  v_zero  <- if (any(!sup)) max(abs(th[!sup])) else 0     # Theta = 0 off-graph
  max(diag_v, v_match, v_zero)
}

#' EBIC-regularized Gaussian graphical model (graphical lasso)
#'
#' Selects an L1 penalty by the extended BIC (Foygel & Drton 2010) over a
#' log-spaced path, then refits the chosen penalty to machine precision so the
#' returned network is the certified global optimum of the convex objective.
#' Equivalent in purpose to `qgraph::EBICglasso()` / `bootnet`'s `"EBICglasso"`
#' default, but pure base R and self-certified (see [glasso_kkt()]).
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix; if given, `n` is required and
#'   `data` is ignored.
#' @param n Sample size (required when `cor_matrix` is supplied).
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param nlambda Number of penalties on the path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#'   Default 0.01.
#' @param threshold Partial correlations with absolute value below this are set
#'   to zero. Default 0.
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default, pairwise-complete correlations + nearest-PD projection) or
#'   `"listwise"` (drop incomplete rows). Identical for complete data.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the partial-correlation matrix,
#'   with `$precision`, `$lambda`, `$gamma`, `$cor_matrix`, `$ebic`, and `$kkt`
#'   (the stationarity residual of the returned network).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 250)
#' fit
#' as.data.frame(fit)
#' @export
ebic_glasso <- function(data = NULL, cor_matrix = NULL, n = NULL,
                        gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,
                        threshold = 0, na_method = c("pairwise", "listwise"),
                        labels = NULL) {
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, na_method = na_method)
    S <- ci$S; n <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- as.matrix(cor_matrix)
    if (is.null(n)) stop("`n` is required when `cor_matrix` is supplied.",
                         call. = FALSE)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
  }
  stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma >= 0,
            nlambda >= 2L, lambda_min_ratio > 0, lambda_min_ratio < 1)

  lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
  sel <- .select_ebic(S, lambda_path, n, gamma)

  pcor <- .precision_to_pcor(sel$wi)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)
  dimnames(S) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "EBICglasso",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = sel$wi, lambda = sel$lambda, gamma = gamma,
      cor_matrix = S, ebic = sel$ebic, ebic_path = sel$ebic_path,
      lambda_path = lambda_path, na_method = na_method,
      kkt = glasso_kkt(sel$wi, S, sel$lambda)
    )
  )
}
