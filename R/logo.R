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
  Theta <- matrix(0, p, p)
  for (cl in cliques) {
    Theta[cl, cl] <- Theta[cl, cl] + solve(S[cl, cl, drop = FALSE])
  }
  for (sp in separators) {
    Theta[sp, sp] <- Theta[sp, sp] - solve(S[sp, sp, drop = FALSE])
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
#' @param method Correlation method when `data` is supplied.
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed. Default 0.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the partial-correlation matrix,
#'   with `$precision`, `$support` (the TMFG graph), `$cor_matrix`, and `$kkt`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(300 * 6), 300, 6)
#' logo_network(x)
#' @export
logo_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         method = c("pearson", "spearman", "kendall"),
                         threshold = 0, labels = NULL) {
  method <- match.arg(method)
  if (is.null(cor_matrix)) {
    mat <- .as_numeric_matrix(data)
    S   <- stats::cor(mat, method = method)
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
  if (p < 4L) stop("LoGo requires at least 4 variables.", call. = FALSE)

  built <- .tmfg_build(abs(S))
  theta <- .logo_precision(S, built$cliques, built$separators)
  support <- built$adj
  dimnames(theta) <- dimnames(support) <- dimnames(S) <- list(labels, labels)

  pcor <- .precision_to_pcor(theta)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)

  .new_psychnet(
    graph = pcor, labels = labels, method = "LoGo",
    directed = FALSE, n_obs = n,
    extra = list(precision = theta, support = support, cor_matrix = S,
                 kkt = ggm_support_kkt(theta, S, support))
  )
}
