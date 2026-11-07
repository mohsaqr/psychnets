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
  # A block that solve() can still invert but whose reciprocal condition number
  # is below 1e-10 is refused too: its inverse has lost most of its digits, and
  # the resulting precision fails its own certificate (a near-duplicate column
  # gave a residual of 0.14 with no condition before this check).
  inv_block <- function(idx) {
    B <- S[idx, idx, drop = FALSE]
    rc <- rcond(B)
    out <- if (rc >= 1e-10) tryCatch(solve(B), error = function(e) NULL)
    if (is.null(out)) {
      nm <- if (is.null(colnames(S))) as.character(idx) else colnames(S)[idx]
      .psn_stop(sprintf(paste0(
        "LoGo: the covariance block of variables {%s} is singular or nearly ",
        "so (reciprocal condition number %.1e): collinear or duplicate ",
        "variables. Remove one of them; the chordal precision cannot be ",
        "formed."), paste(nm, collapse = ", "), rc),
        class = c("psychnet_singular", "psychnet_bad_input"))
    }
    out
  }
  Theta <- matrix(0, p, p)
  for (cl in cliques) {
    Theta[cl, cl] <- Theta[cl, cl] + inv_block(cl)
  }
  for (sp in separators) {
    Theta[sp, sp] <- Theta[sp, sp] - inv_block(sp)
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
#' @param n Sample size (recorded on the result; required with `cor_matrix`): a
#'   single finite number >= 1.
#' @param cor_method Correlation when `data` is supplied: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see [cor_auto()]).
#' @param threshold Partial correlations with absolute value below this are
#'   zeroed: a single number >= 0, default 0. With a positive value `$kkt` is `NA` for the returned
#'   graph and the pre-threshold residual is stored as `$fit_kkt`.
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels, one unique label per network node.
#' @return A `psychnet` object whose `$weights` is the partial-correlation matrix,
#'   with `$precision`, `$support` (the TMFG graph), `$cor_matrix`, and `$kkt`.
#'   A clique or separator block that is singular or nearly so (reciprocal
#'   condition number below 1e-10, e.g. a duplicate column), or a precision
#'   whose certificate exceeds 1e-6, is refused with an error of class
#'   `psychnet_singular` rather than returned as an uncertified network.
#' @references Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#'   Parsimonious modeling with information filtering networks. *Physical
#'   Review E*, 94(6), 062306. \doi{10.1103/PhysRevE.94.062306}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(300 * 6), 300, 6)
#' logo_network(x)
#' @export
logo_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         threshold = 0, na_method = c("pairwise", "listwise"),
                         labels = NULL) {
  cor_method <- .psn_match_arg(cor_method)
  na_method <- .psn_match_arg(na_method)
  .check_threshold(threshold)
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "required")
  S <- inp$S; n <- inp$n; labels <- inp$labels
  p <- ncol(S)
  if (p < 4L) .psn_stop("LoGo requires at least 4 variables.", call. = FALSE)

  built <- .tmfg_build(abs(S))
  theta <- .logo_precision(S, built$cliques, built$separators)
  support <- built$adj
  dimnames(theta) <- dimnames(support) <- dimnames(S) <- list(labels, labels)

  pcor <- .precision_to_pcor(theta)
  pcor[abs(pcor) < threshold] <- 0
  dimnames(pcor) <- list(labels, labels)
  fit_kkt <- ggm_support_kkt(theta, S, support)
  # Refuse rather than approximate: the closed form is exact on a chordal
  # graph, so a residual this large means the blocks were too ill-conditioned
  # for the inverse to be trusted even though each passed the rcond check.
  if (!is.finite(fit_kkt) || fit_kkt > 1e-6) {
    .psn_stop(sprintf(paste0(
      "LoGo: the chordal precision fails its own certificate (residual %.2e > ",
      "1e-6); the correlation matrix is too ill-conditioned (collinear or ",
      "near-duplicate variables) for the closed-form inverse."), fit_kkt),
      class = c("psychnet_singular", "psychnet_bad_input"))
  }

  .new_psychnet(
    graph = pcor, labels = labels, method = "logo",
    directed = FALSE, n_obs = n,
    extra = list(precision = theta, support = support, cor_matrix = S,
                 kkt = if (threshold > 0) NA_real_ else fit_kkt,
                 fit_kkt = fit_kkt,
                 certificate_target = if (threshold > 0)
                   "pre_threshold_fit" else "returned_network")
  )
}
