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
  na_method <- .psn_match_arg(na_method)
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
    .psn_stop("Pairwise nonparanormal correlation is undefined: a pair is never co-observed.",
         call. = FALSE)
  }
  if (has_na) S <- .nearest_pd_cor(S)
  n <- if (has_na) {
    co <- crossprod(!is.na(mat)); max(round(stats::median(co[upper.tri(co)])), 2L)
  } else nrow(mat)
  list(S = S, n = n)
}

#' Estimate a nonparanormal Gaussian graphical model
#'
#' @description
#' `huge_network()` estimates a nonparanormal Gaussian graphical model: a
#' regularized partial-correlation network for variables that are normal only
#' after an unknown monotone transformation. The nonparanormal model (Liu et
#' al., 2009) relaxes the assumption of the Gaussian graphical model that the
#' variables are jointly normal, an assumption that skewed or heavy-tailed
#' questionnaire and symptom data often violate, so the network of partial
#' correlations is defined on the transformed scale and is unaffected by the
#' shape of each variable's distribution.
#'
#' The transformation is estimated from the ranks of the observations.
#' `huge_network()` transforms the data, computes the correlation matrix of
#' the transformed variables and estimates a regularized partial-correlation
#' network on it with the EBIC graphical lasso, as in [ebic_glasso()]. `npn`
#' chooses the transformation: normal scores of the ranks with shrinkage
#' (default) or truncation, or the "skeptic", a correlation derived from
#' Spearman's rank correlation.
#'
#' The returned network is the optimum of the penalized likelihood at the
#' selected penalty, which [certificate()] checks.
#'
#' For approximately normal data, [ebic_glasso()] fits the same model without
#' the transformation; its `cor_method = "auto"` treats ordinal items through
#' polychoric correlations instead. The model is the one fitted by
#' `huge::huge()` with the nonparanormal transformation.
#'
#' @details
#' The shrinkage transform replaces each observation by
#' \eqn{\Phi^{-1}(r / (n + 1))}, with \eqn{r} its rank among the \eqn{n}
#' observed values of the variable. The truncation transform uses
#' \eqn{\Phi^{-1}(r / n)} truncated to \eqn{[\delta, 1 - \delta]}, with
#' \eqn{\delta = 1 / (4 n^{1/4} \sqrt{\pi \log n})}. The skeptic uses
#' \eqn{2 \sin(\pi \rho_S / 6)}, with \eqn{\rho_S} Spearman's correlation. A
#' correlation matrix that is not positive definite (possible for the skeptic
#' and for pairwise data) is adjusted to the nearest positive-definite matrix
#' before estimation.
#'
#' The penalty path, the EBIC and the certificate are those of
#' [ebic_glasso()] with its default settings; the certificate is computed by
#' [glasso_kkt()] on the transformed correlation matrix. After a positive
#' `threshold` the certificate is missing; `print()` still shows the value for
#' the network before thresholding.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#'   Not needed when `cor_matrix` is given.
#' @param cor_matrix A correlation matrix of already transformed data, used
#'   instead of `data`; no transformation is applied and `n` is required.
#' @param n Sample size, required with `cor_matrix`: a single finite number
#'   >= 1. Ignored when `data` is given.
#' @param npn The transformation: `"shrinkage"` (default), `"truncation"` or
#'   `"skeptic"`. Ignored with `cor_matrix`.
#' @param gamma EBIC hyperparameter, a number >= 0; larger values favour
#'   sparser networks. Default `0.5`.
#' @param nlambda Number of penalties scanned, a whole number >= 2. Default
#'   `100`.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest, in
#'   (0, 1). Default `0.01`.
#' @param threshold Absolute partial correlation below which an edge is set to
#'   zero after estimation: a single number >= 0. Default `0`. A logical value
#'   is refused. A positive value leaves the certificate missing.
#' @param na_method Handling of missing values in `data`: `"pairwise"`
#'   (default) transforms each variable over its observed values and uses
#'   pairwise-complete correlations; `"listwise"` drops incomplete rows first.
#' @param native `TRUE` (default) uses the package's own graphical lasso
#'   solver; `FALSE` uses the `glasso` package, which must be installed. See
#'   [ebic_glasso()].
#' @param labels Optional node labels, one unique label per node of the
#'   network.
#' @return A fitted network of class `psychnet` ([psychnet-class]) whose edge
#'   weights are partial correlations on the transformed scale.
#'   `as.data.frame()` returns one row per edge with columns `from`, `to` and
#'   `weight`; [certificate()] reports the check. The object also stores the
#'   precision matrix (`precision`), the selected penalty (`lambda`),
#'   `gamma`, the transformed correlation matrix (`cor_matrix`), the
#'   transformation used (`npn`, `"precomputed"` with `cor_matrix`) and the
#'   EBIC of the selected network (`ebic`).
#' @section Conditions:
#'   `psychnet_bad_input` when an argument is out of range, `threshold` is
#'   logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
#'   matrix with a positive diagonal, `n` is missing with `cor_matrix`,
#'   `data` has fewer than two usable variables, a pair of variables is never
#'   observed together, `labels` do not match the nodes, or `native = FALSE`
#'   without the `glasso` package; `psychnet_constant_variable` (warning)
#'   when a constant or all-missing column is dropped.
#' @seealso [ebic_glasso()] for the untransformed model; [certificate()] and
#'   [glasso_kkt()] for the check; [psychnet()], which routes
#'   `method = "huge"` here; [psychnet-class] for the result.
#' @references Liu, H., Lafferty, J., & Wasserman, L. (2009). The
#'   nonparanormal: Semiparametric estimation of high dimensional undirected
#'   graphs. *Journal of Machine Learning Research*, 10, 2295-2328.
#'
#'   Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for
#'   Gaussian graphical models. *Advances in Neural Information Processing
#'   Systems*, 23, 604-612.
#'
#'   Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#'   covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
#'   432-441. \doi{10.1093/biostatistics/kxm045}
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
  npn <- .psn_match_arg(npn)
  na_method <- .psn_match_arg(na_method)
  engine <- .resolve_native(native, "glasso")
  .check_ebic_args(gamma, nlambda, lambda_min_ratio, threshold)
  if (is.null(cor_matrix)) {
    mat <- .as_numeric_matrix(data, drop_na = FALSE)
    .psn_check_labels(labels, ncol(mat))
    if (is.null(labels)) labels <- colnames(mat)
    nc  <- .npn_cor(mat, npn, na_method); S <- nc$S; n <- nc$n
  } else {
    inp <- .network_input(NULL, cor_matrix, n, "pearson", na_method, labels,
                          n_role = "required")
    S <- inp$S; n <- inp$n; labels <- inp$labels
    npn <- "precomputed"
  }

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
  fit_kkt <- glasso_kkt(sel$wi, S, sel$lambda)

  .new_psychnet(
    graph = pcor, labels = labels, method = "huge",
    directed = FALSE, n_obs = n,
    extra = list(
      precision = sel$wi, lambda = sel$lambda, gamma = gamma,
      cor_matrix = S, npn = npn, ebic = sel$ebic,
      kkt = if (threshold > 0) NA_real_ else fit_kkt,
      fit_kkt = fit_kkt,
      certificate_target = if (threshold > 0) "pre_threshold_fit" else
        "returned_network"
    )
  )
}
