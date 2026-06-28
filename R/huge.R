# Nonparanormal graphical model (Liu, Lafferty & Wasserman 2009; Liu et al.
# 2012), clean-room base R. The nonparanormal relaxes the Gaussian assumption:
# each variable is assumed Gaussian only after an unknown monotone transform.
# We estimate that transform by normal-scoring the ranks (or via a rank
# correlation for the "skeptic"), then hand the resulting correlation matrix to
# the existing EBIC graphical-lasso path -- so the certificate is the ordinary
# glasso_kkt() on the transformed correlation. Equivalent in purpose to
# huge::huge() with the nonparanormal transform.

# Per-column nonparanormal normal-score transform (NA-preserving): ranks each
# column over its observed values. With no missing data this reproduces the
# classic `qnorm(rank/(n+1))` transform exactly.
#' @noRd
.npn_scores <- function(mat, npn) {
  Z <- vapply(seq_len(ncol(mat)), function(j) {
    x <- mat[, j]; ok <- !is.na(x); nj <- sum(ok); out <- x
    r <- rank(x[ok])
    out[ok] <- if (npn == "shrinkage") {
      stats::qnorm(r / (nj + 1))
    } else {                                       # truncation
      th <- 1 / (4 * nj^0.25 * sqrt(pi * log(nj)))
      stats::qnorm(pmin(pmax(r / nj, th), 1 - th))
    }
    out
  }, numeric(nrow(mat)))
  colnames(Z) <- colnames(mat)
  Z
}

# Nonparanormal correlation matrix + effective n from a (possibly incomplete)
# data matrix, honouring `na_method` (see .cor_input).
#' @noRd
.npn_cor <- function(mat, npn, na_method = c("pairwise", "listwise")) {
  na_method <- match.arg(na_method)
  if (anyNA(mat) && na_method == "listwise") {
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  }
  has_na <- anyNA(mat)
  use <- if (has_na) "pairwise.complete.obs" else "everything"
  if (npn == "skeptic") {
    S <- 2 * sin(pi / 6 * stats::cor(mat, use = use, method = "spearman"))
    diag(S) <- 1
  } else {
    S <- stats::cor(.npn_scores(mat, npn), use = use)
  }
  if (anyNA(S)) {
    stop("Pairwise nonparanormal correlation is undefined: a pair is never co-observed.",
         call. = FALSE)
  }
  if (has_na) S <- .nearest_pd_cor(S)
  n <- if (has_na) {
    co <- crossprod(!is.na(mat)); max(round(stats::median(co[upper.tri(co)])), 2L)
  } else nrow(mat)
  list(S = S, n = n)
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
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default, with the nonparanormal transform applied per column over observed
#'   values) or `"listwise"`. See [ebic_glasso()].
#' @param native Solver switch for the glasso path: `TRUE` (default) uses the
#'   pure-R solver; `FALSE` delegates to the `glasso` Fortran package (in
#'   `Suggests`). See [ebic_glasso()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the partial-correlation matrix,
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
                         threshold = 0, na_method = c("pairwise", "listwise"),
                         native = TRUE, labels = NULL) {
  npn <- match.arg(npn)
  na_method <- match.arg(na_method)
  engine <- .resolve_native(native, "glasso")
  if (is.null(cor_matrix)) {
    mat <- .as_numeric_matrix(data, drop_na = FALSE)
    if (is.null(labels)) labels <- colnames(mat)
    nc  <- .npn_cor(mat, npn, na_method); S <- nc$S; n <- nc$n
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (is.null(n)) stop("`n` is required when `cor_matrix` is supplied.",
                         call. = FALSE)
    stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n > 0)
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

  if (max(abs(S[upper.tri(S)])) <= 1e-12) {     # no association: empty optimum
    sel <- .empty_glasso(S, n)
  } else {
    lambda_path <- .compute_lambda_path(S, nlambda, lambda_min_ratio)
    sel <- .select_ebic(S, lambda_path, n, gamma, engine = engine)
  }

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
