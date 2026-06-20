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

# EBIC of an unregularized constrained fit (npar = number of retained edges).
#' @noRd
.support_ebic <- function(Theta, S, n, gamma, p) {
  ld <- determinant(Theta, logarithm = TRUE)
  if (ld$sign <= 0) return(Inf)
  loglik <- (n / 2) * (as.numeric(ld$modulus) - sum(diag(S %*% Theta)))
  npar   <- sum(Theta[upper.tri(Theta)] != 0)
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
#' @param n Sample size (required when `cor_matrix` is supplied).
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param stepwise If `TRUE` (default), refine the best glasso-path graph by a
#'   greedy single-edge add/drop search.
#' @param nlambda Number of glasso penalties scanned for candidate graphs.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed. Default 0.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the partial-correlation matrix,
#'   with `$precision`, `$support` (the selected graph), `$gamma`, `$ebic`,
#'   `$cor_matrix`, and `$kkt`.
#' @examples
#' S <- 0.5^abs(outer(1:6, 1:6, "-"))
#' ggm_modselect(cor_matrix = S, n = 250)
#' @export
ggm_modselect <- function(data = NULL, cor_matrix = NULL, n = NULL,
                          gamma = 0.5, stepwise = TRUE, nlambda = 100L,
                          lambda_min_ratio = 0.01, threshold = 0,
                          labels = NULL) {
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
  p <- ncol(S)

  # Candidate supports from the glasso path (deduplicated by edge pattern).
  lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
  seen <- character(0)
  best_ebic <- Inf; best_support <- NULL; best_theta <- NULL
  for (lam in lambda_path) {
    wi <- tryCatch(.glasso_fit(S, lam, tol_outer = 1e-4, tol_inner = 1e-4)$wi,
                   error = function(e) NULL)
    if (is.null(wi)) next
    sup <- (abs(wi) > 1e-6) & upper.tri(wi)
    sup <- sup | t(sup)
    key <- paste(which(sup[upper.tri(sup)]), collapse = ",")
    if (key %in% seen) next
    seen <- c(seen, key)
    theta <- tryCatch(.ggm_fit_support(S, sup), error = function(e) NULL)
    if (is.null(theta)) next
    eb <- .support_ebic(theta, S, n, gamma, p)
    if (eb < best_ebic) { best_ebic <- eb; best_support <- sup; best_theta <- theta }
  }
  if (is.null(best_support)) stop("Model selection failed.", call. = FALSE)

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
        theta <- tryCatch(.ggm_fit_support(S, cand), error = function(e) NULL)
        if (is.null(theta)) next
        eb <- .support_ebic(theta, S, n, gamma, p)
        if (eb < best_ebic - 1e-10) {
          best_ebic <- eb; best_support <- cand; best_theta <- theta
          improved <- TRUE
        }
      }
    }
  }

  best_theta <- .ggm_fit_support(S, best_support)   # tight final fit
  dimnames(best_theta) <- dimnames(best_support) <- dimnames(S) <- list(labels, labels)
  pcor <- .precision_to_pcor(best_theta)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "ggmModSelect",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = best_theta, support = best_support, gamma = gamma,
      cor_matrix = S, ebic = best_ebic,
      kkt = ggm_support_kkt(best_theta, S, best_support)
    )
  )
}
