# Nonparanormal graphical model (Liu, Lafferty & Wasserman 2009; Liu et al.
# 2012), clean-room base R. The nonparanormal relaxes the Gaussian assumption:
# each variable is assumed Gaussian only after an unknown monotone transform.
# We estimate that transform by normal-scoring the ranks (or via a rank
# correlation for the "skeptic"), then hand the resulting correlation matrix to
# the existing EBIC graphical-lasso path -- so the certificate is the ordinary
# glasso_kkt() on the transformed correlation. Equivalent in purpose to
# huge::huge() with the nonparanormal transform.

# Nonparanormal correlation matrix from a raw data matrix.
#' @noRd
.npn_cor <- function(mat, npn) {
  n <- nrow(mat)
  if (npn == "skeptic") {
    # Spearman skeptic: S = 2 sin(pi/6 rho_s) (Liu et al. 2012).
    S <- 2 * sin(pi / 6 * stats::cor(mat, method = "spearman"))
    diag(S) <- 1
    return(S)
  }
  thresh <- 1 / (4 * n^0.25 * sqrt(pi * log(n)))
  R <- apply(mat, 2L, rank)
  Z <- if (npn == "shrinkage") {
    stats::qnorm(R / (n + 1))                      # avoids qnorm(1) = Inf
  } else {                                         # truncation
    stats::qnorm(pmin(pmax(R / n, thresh), 1 - thresh))
  }
  stats::cor(Z)                                    # cor is scale-invariant
}

#' Nonparanormal graphical model (huge)
#'
#' Estimates a Gaussian graphical model after a rank-based nonparanormal
#' transform that relaxes the multivariate-normal assumption, then selects the
#' L1 penalty by EBIC and refits to the certified optimum. Equivalent in purpose
#' to `huge::huge()` (nonparanormal) / `bootnet`'s `"huge"` default, but pure
#' base R and self-certified via [glasso_kkt()] on the transformed correlation.
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied (then the transform is skipped).
#' @param cor_matrix Optional pre-transformed correlation matrix; if given, `n`
#'   is required, `data` and `npn` are ignored.
#' @param n Sample size (required when `cor_matrix` is supplied).
#' @param npn Nonparanormal transform: `"shrinkage"` (default), `"truncation"`,
#'   or `"skeptic"` (Spearman).
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param nlambda Number of penalties on the path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed. Default 0.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the partial-correlation matrix,
#'   with `$precision`, `$lambda`, `$gamma`, `$cor_matrix` (the transformed
#'   correlation), `$npn`, `$ebic`, and `$kkt`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(300 * 5), 300, 5)
#' x <- exp(x %*% chol(0.4^abs(outer(1:5, 1:5, "-"))))   # break normality
#' huge_network(x)
#' @export
huge_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         npn = c("shrinkage", "truncation", "skeptic"),
                         gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,
                         threshold = 0, labels = NULL) {
  npn <- match.arg(npn)
  if (is.null(cor_matrix)) {
    mat <- .as_numeric_matrix(data)
    n   <- nrow(mat)
    if (is.null(labels)) labels <- colnames(mat)
    S   <- .npn_cor(mat, npn)
  } else {
    S <- as.matrix(cor_matrix)
    if (is.null(n)) stop("`n` is required when `cor_matrix` is supplied.",
                         call. = FALSE)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
    npn <- "precomputed"
  }
  stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma >= 0,
            nlambda >= 2L, lambda_min_ratio > 0, lambda_min_ratio < 1)

  # The skeptic correlation is not guaranteed positive definite; project it.
  if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < 1e-8) {
    S <- .nearest_pd_cor(S)
  }

  lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
  sel <- .select_ebic(S, lambda_path, n, gamma)

  pcor <- .precision_to_pcor(sel$wi)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)
  dimnames(S) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "huge",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = sel$wi, lambda = sel$lambda, gamma = gamma,
      cor_matrix = S, npn = npn, ebic = sel$ebic,
      kkt = glasso_kkt(sel$wi, S, sel$lambda)
    )
  )
}
