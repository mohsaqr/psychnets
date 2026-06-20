# Non-convex regularized Gaussian graphical model (Williams 2020, GGMncv),
# clean-room base R. A non-convex penalty (SCAD, MCP, or atan) is solved by a
# single one-step local linear approximation (Zou & Li 2008; Fan, Xue & Zou
# 2014): fit the L1 (lasso) glasso, read off per-edge weights from the penalty
# derivative at that solution, then take ONE weighted-glasso step. The weighted
# problem is convex, so glasso_kkt_weighted() certifies the returned network.

# Per-penalty default shape hyperparameter.
#' @noRd
.ncv_default_hyper <- function(penalty, hyper) {
  if (!is.null(hyper)) return(hyper)
  switch(penalty, scad = 3.7, mcp = 3, atan = 0.5, lasso = NA_real_)
}

# One-step LLA weights w_ij = p'_lambda(|theta_init_ij|) / lambda, as a matrix.
#' @noRd
.ncv_weights <- function(t_abs, lambda, penalty, hyper) {
  if (penalty == "lasso") {
    return(matrix(1, nrow(t_abs), ncol(t_abs)))
  }
  if (penalty == "scad") {
    a <- hyper
    w <- ifelse(t_abs <= lambda, 1,
                pmax(a * lambda - t_abs, 0) / ((a - 1) * lambda))
  } else if (penalty == "mcp") {
    g <- hyper
    w <- pmax(1 - t_abs / (g * lambda), 0)
  } else {                                            # atan (Wang & Zhu 2016)
    g <- hyper
    w <- (g + 2 / pi) * g / (g^2 + t_abs^2)
  }
  w
}

# One non-convex fit at fixed lambda: lasso init -> LLA weights -> weighted fit.
#' @noRd
.ggmncv_fit <- function(S, lambda, penalty, hyper, tol) {
  init <- .glasso_fit(S, lambda, tol_outer = tol, tol_inner = tol)$wi
  W <- .ncv_weights(abs(init), lambda, penalty, hyper)
  Rho <- lambda * W
  diag(Rho) <- 0
  fit <- .glasso_fit(S, Rho, tol_outer = tol, tol_inner = tol * 1e-2)
  list(wi = fit$wi, Rho = Rho)
}

#' Non-convex regularized Gaussian graphical model (GGMncv)
#'
#' Estimates a GGM with a non-convex penalty (SCAD, MCP, or atan) via one-step
#' local linear approximation, selecting the penalty by EBIC. Non-convex
#' penalties reduce the shrinkage bias of the L1 (lasso) glasso on strong edges.
#' Equivalent in purpose to `GGMncv::ggmncv()`, but pure base R and self-
#' certified via [glasso_kkt_weighted()].
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix; if given, `n` is required.
#' @param n Sample size (required when `cor_matrix` is supplied).
#' @param penalty `"atan"` (default), `"scad"`, `"mcp"`, or `"lasso"`.
#' @param hyper Penalty shape hyperparameter (SCAD `a`, MCP `gamma`, atan
#'   `gamma`); per-penalty default if `NULL`.
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param nlambda Number of penalties on the path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed. Default 0.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the partial-correlation matrix,
#'   with `$precision`, `$lambda`, `$gamma`, `$penalty`, `$hyper`,
#'   `$penalty_matrix` (the per-edge penalties), `$cor_matrix`, `$ebic`, and
#'   `$kkt`.
#' @examples
#' S <- 0.5^abs(outer(1:6, 1:6, "-"))
#' ggmncv_network(cor_matrix = S, n = 250, penalty = "scad")
#' @export
ggmncv_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                           penalty = c("atan", "scad", "mcp", "lasso"),
                           hyper = NULL, gamma = 0.5, nlambda = 100L,
                           lambda_min_ratio = 0.01, threshold = 0,
                           labels = NULL) {
  penalty <- match.arg(penalty)
  if (is.null(cor_matrix)) {
    mat <- .as_numeric_matrix(data)
    S   <- stats::cor(mat)
    n   <- nrow(mat)
    if (is.null(labels)) labels <- colnames(mat)
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
  hyper <- .ncv_default_hyper(penalty, hyper)
  p <- ncol(S)

  lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
  best_ebic <- Inf; best_lambda <- lambda_path[1L]
  for (lam in lambda_path) {
    fit <- tryCatch(.ggmncv_fit(S, lam, penalty, hyper, tol = 1e-4),
                    error = function(e) NULL)
    if (is.null(fit)) next
    ld <- determinant(fit$wi, logarithm = TRUE)
    if (ld$sign <= 0) next
    loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% fit$wi)))
    npar   <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic   <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)
    if (ebic < best_ebic) { best_ebic <- ebic; best_lambda <- lam }
  }

  fit <- .ggmncv_fit(S, best_lambda, penalty, hyper, tol = 1e-8)
  dimnames(fit$wi) <- dimnames(fit$Rho) <- dimnames(S) <- list(labels, labels)

  pcor <- .precision_to_pcor(fit$wi)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "GGMncv",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = fit$wi, lambda = best_lambda, gamma = gamma,
      penalty = penalty, hyper = hyper, penalty_matrix = fit$Rho,
      cor_matrix = S, ebic = best_ebic,
      kkt = glasso_kkt_weighted(fit$wi, S, fit$Rho)
    )
  )
}
