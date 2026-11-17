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

#' Estimate a sparse Gaussian graphical model on a filtered graph (LoGo)
#'
#' @description
#' `logo_network()` estimates a sparse Gaussian graphical model on a filtered
#' graph with the Local-Global (LoGo) method. LoGo works in two steps: it
#' first chooses the network structure by information filtering with the
#' triangulated maximally filtered graph (TMFG), and then estimates the
#' partial correlations as the maximum-likelihood Gaussian graphical model on
#' that graph (Barfuss et al., 2016).
#'
#' The TMFG retains the strongest correlations under the constraint that the
#' graph stays planar and is built from tetrahedra (four-variable cliques),
#' which fixes the number of edges at \eqn{3(p - 2)} for \eqn{p} variables.
#' Because the graph is decomposable, the maximum-likelihood estimate has a
#' closed form, assembled from the inverses of small local correlation
#' matrices, so no penalty or iterative fitting is needed. The model is
#' estimated from raw data or from a correlation matrix with its sample size.
#'
#' `logo_network()` returns the partial-correlation network of this model, the
#' maximum-likelihood network on the TMFG graph, which [certificate()] checks.
#'
#' The graph alone, with correlations as edge weights, is returned by
#' [tmfg_network()]. For a network whose sparsity is chosen by the data use
#' [ebic_glasso()] or [ggm_modselect()]. The model is the one fitted by
#' `NetworkToolbox::LoGo()`.
#'
#' @details
#' With \eqn{C} the four-variable cliques and \eqn{T} the three-variable
#' separators of the TMFG, the precision matrix is
#' \deqn{\Theta = \sum_{c \in C} [S_{cc}^{-1}]^0 - \sum_{t \in T} [S_{tt}^{-1}]^0,}
#' where \eqn{[\cdot]^0} pads a block with zeros to the full size. The
#' certificate is computed by [ggm_support_kkt()] with the TMFG as the graph.
#' After a positive `threshold` it is missing; `print()` still shows the value
#' for the network before thresholding.
#'
#' The closed form needs every clique and separator block to be invertible.
#' A block that is singular or nearly so (reciprocal condition number below
#' `1e-10`, as with duplicate or collinear variables) is refused rather than
#' approximated, as is a precision matrix whose certificate exceeds `1e-6`.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#'   Not needed when `cor_matrix` is given.
#' @param cor_matrix A correlation matrix to use instead of `data`; a
#'   covariance matrix is rescaled to a correlation matrix. `n` is then
#'   required.
#' @param n Sample size, required with `cor_matrix` and recorded on the
#'   result: a single finite number >= 1. Ignored when `data` is given.
#' @param cor_method Correlation computed from `data`: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric and polyserial
#'   correlations for ordinal items; see [cor_auto()]).
#' @param threshold Absolute partial correlation below which an edge is set to
#'   zero after estimation: a single number >= 0. Default `0`. A logical value
#'   is refused. A positive value leaves the certificate missing.
#' @param na_method Handling of missing values in `data`: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels, one unique label per node of the
#'   network.
#' @return A fitted network of class `psychnet` ([psychnet-class]) whose edge
#'   weights are partial correlations. `as.data.frame()` returns one row per
#'   edge with columns `from`, `to` and `weight`; [certificate()] reports the
#'   check. The object also stores the precision matrix (`precision`), the
#'   TMFG graph as a logical matrix (`support`) and the correlation matrix
#'   (`cor_matrix`).
#' @section Conditions:
#'   `psychnet_singular` (error, also of class `psychnet_bad_input`) when a
#'   clique or separator block is singular or nearly so, or the precision
#'   matrix fails its certificate; `psychnet_bad_input` when there are fewer
#'   than four variables, an argument is out of range, `threshold` is
#'   logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
#'   matrix with a positive diagonal, `n` is missing with `cor_matrix`, or
#'   `labels` do not match the nodes; `psychnet_constant_variable` (warning)
#'   when a constant or all-missing column is dropped.
#' @seealso [tmfg_network()] for the filtered graph; [certificate()] and
#'   [ggm_support_kkt()] for the check; [psychnet()], which routes
#'   `method = "logo"` here; [psychnet-class] for the result.
#' @references Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#'   Parsimonious modeling with information filtering networks. *Physical
#'   Review E*, 94(6), 062306. \doi{10.1103/PhysRevE.94.062306}
#'
#'   Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network filtering for big
#'   data: Triangulated maximally filtered graph. *Journal of Complex Networks*,
#'   5(2), 161-178. \doi{10.1093/comnet/cnw015}
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
