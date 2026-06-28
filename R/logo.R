# LoGo: Local-Global sparse inverse covariance (Barfuss, Massara, Di Matteo &
# Aste 2016), clean-room base R. Given the chordal TMFG, the sparse precision is
# the closed-form maximum-likelihood (max-entropy) Gaussian Markov random field
# on that graph: sum of the inverse covariances of the 4-cliques minus the
# inverse covariances of the 3-clique separators, each zero-padded to p x p.
# Because the graph is chordal this reconstruction is exact, so W = Theta^{-1}
# matches S on the TMFG support -- certified by ggm_support_kkt().

# Assemble the LoGo precision from the TMFG clique / separator decomposition.
#' @noRd
.logo_precision <- function(S, cliques, separators) {
  p <- ncol(S)
  # The LoGo decomposition K = sum inv(S_clique) - sum inv(S_separator) is only
  # defined when every clique and separator covariance block is invertible. A
  # pseudo-inverse does NOT give the Gaussian decomposable MLE, so a singular
  # block (duplicate or collinear columns) is a hard error rather than a guess.
  inv <- function(B) tryCatch(solve(B), error = function(e)
    stop("LoGo: a clique or separator covariance block is singular (collinear ",
         "or duplicate variables); cannot form the chordal precision.",
         call. = FALSE))
  Theta <- matrix(0, p, p)
  for (cl in cliques) {
    Theta[cl, cl] <- Theta[cl, cl] + inv(S[cl, cl, drop = FALSE])
  }
  for (sp in separators) {
    Theta[sp, sp] <- Theta[sp, sp] - inv(S[sp, sp, drop = FALSE])
  }
  (Theta + t(Theta)) / 2
}

#' Local-Global sparse inverse covariance (LoGo)
#'
#' Estimates a sparse Gaussian graphical model whose conditional-independence
#' structure is the chordal TMFG: the precision is the closed-form Gaussian
#' Markov random field on that graph (Barfuss et al. 2016). Equivalent in
#' purpose to `NetworkToolbox::LoGo()` / `bootnet`'s `"LoGo"` default, pure base
#' R and self-certified via [ggm_support_kkt()] (the precision reproduces `S`
#' exactly on the TMFG support).
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix; if given, `n` is required.
#' @param n Sample size (recorded on the result; required with `cor_matrix`).
#' @param cor_method Correlation when `data` is supplied: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see [cor_auto()]).
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed. Default 0.
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the partial-correlation matrix,
#'   with `$precision`, `$support` (the TMFG graph), `$cor_matrix`, and `$kkt`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(300 * 6), 300, 6)
#' logo_network(x)
#' @export
logo_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         threshold = 0, na_method = c("pairwise", "listwise"),
                         labels = NULL) {
  cor_method <- match.arg(cor_method)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (is.null(n)) stop("`n` is required when `cor_matrix` is supplied.",
                         call. = FALSE)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
  }
  p <- ncol(S)
  if (p < 4L) stop("LoGo requires at least 4 variables.", call. = FALSE)

  built <- .tmfg_build(abs(S))
  theta <- .logo_precision(S, built$cliques, built$separators)
  support <- built$adj
  dimnames(theta) <- dimnames(support) <- dimnames(S) <- list(labels, labels)

  pcor <- .precision_to_pcor(theta)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "logo",
    directed = FALSE, n_obs = n,
    extra = list(precision = theta, support = support, cor_matrix = S,
                 kkt = ggm_support_kkt(theta, S, support))
  )
}
