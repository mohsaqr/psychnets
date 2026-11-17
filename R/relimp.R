# Relative-importance network by LMG / Shapley decomposition (Lindeman, Merenda
# & Gold 1980; Gromping 2006), clean-room base R. Each node is regressed on all
# others and its model R-squared is decomposed into per-predictor shares by
# averaging each predictor's marginal contribution over every subset of the
# other predictors (the Shapley value). This produces a DIRECTED network: the
# edge predictor -> outcome carries the predictor's importance share. The
# combinatorial sum over 2^(p-1) subsets is the documented exception to the
# vectorisation rule. Certified by lmg_certificate(): a node's incoming shares
# sum to its full-model R-squared (Shapley efficiency), and that full-model
# R-squared equals the closed form 1 - 1/diag(solve(S)). Efficiency alone
# telescopes for ANY set function (it would pass even if every subset R^2 were
# wrong by the same factor); the closed form anchors the full model
# independently of the subset solver.

# Decomposition residual shared by relimp_network() and lmg_certificate():
# max of (a) |sum of incoming shares - stored R^2| and (b) |stored R^2 - closed
# form 1 - 1/diag(S^-1)|. (b) needs an invertible S; for a singular S (exactly
# collinear predictors) the closed form does not exist and only (a) is checked.
#' @noRd
.lmg_residual <- function(raw_importance, r2, S) {
  efficiency <- max(abs(colSums(raw_importance) - r2))
  rc <- rcond(S)
  if (!is.finite(rc) || rc < 1e-12) return(efficiency)
  r2_closed <- 1 - 1 / diag(solve(S))
  max(efficiency, max(abs(r2 - r2_closed)))
}

# Closed-form subset R-squared from the correlation matrix:
#   R^2_{j|A} = S[j,A] %*% solve(S[A,A]) %*% S[A,j].
# LMG shares for one outcome `j`, aligned to the predictor index `pred`.
#' @noRd
.lmg_node <- function(S, j) {
  p <- ncol(S)
  pred <- setdiff(seq_len(p), j)
  m <- length(pred)
  Spp <- S[pred, pred, drop = FALSE]
  sjp <- S[j, pred]

  r2 <- numeric(2^m)                                   # indexed by mask + 1
  for (mask in seq_len(2^m - 1L)) {
    A <- which(intToBits(mask)[seq_len(m)] != as.raw(0L))
    SA <- Spp[A, A, drop = FALSE]
    sol <- tryCatch(solve(SA, sjp[A]),
                    error = function(e) MASS_ginv(SA) %*% sjp[A])
    r2[mask + 1L] <- sum(sjp[A] * as.numeric(sol))
  }

  # Shapley weights via lgamma to stay finite up to m = 20.
  wt <- exp(lgamma(0:(m - 1L) + 1) + lgamma(m - 0:(m - 1L)) - lgamma(m + 1))
  lmg <- numeric(m)
  for (k in seq_len(m)) {
    bitk <- bitwShiftL(1L, k - 1L)
    for (mask in 0:(2^m - 1L)) {
      if (bitwAnd(mask, bitk) != 0L) next
      s <- sum(intToBits(mask)[seq_len(m)] != as.raw(0L))
      lmg[k] <- lmg[k] + wt[s + 1L] *
        (r2[bitwOr(mask, bitk) + 1L] - r2[mask + 1L])
    }
  }
  list(pred = pred, lmg = lmg, r2_full = r2[2^m])
}

#' Check the decomposition of a relative-importance network
#'
#' @description
#' `lmg_certificate()` checks the decomposition stored in a relative-importance
#' network from [relimp_network()]. Such a network divides the \eqn{R^2} of
#' each variable, regressed on all the others, among its predictors; it is
#' defined by a decomposition rather than as the optimum of an objective, and
#' two identities must hold for that decomposition to be correct.
#'
#' The shares a node receives must sum exactly to its \eqn{R^2}, which is the
#' efficiency property of the Shapley value. And that \eqn{R^2} must be the
#' one implied by the correlation matrix, which has a closed form. Together
#' they show that the stored shares add up to the variance each node's
#' predictors actually explain.
#'
#' `lmg_certificate()` returns the larger of the two deviations: a value near
#' zero confirms the decomposition, up to numerical precision.
#'
#' [certificate()] reports the same value for a fitted network. The check
#' does not verify the \eqn{R^2} of the individual subset
#' regressions between the empty and the full model, and the sum identity
#' alone would hold for any set of subset values; the closed form is what ties
#' the full model to the data.
#'
#' @details
#' For outcome \eqn{j} with shares \eqn{w_{kj}} and full-model \eqn{R^2_j},
#' the two deviations are \eqn{\max_j |\sum_k w_{kj} - R^2_j|} and
#' \eqn{\max_j |R^2_j - (1 - 1 / (S^{-1})_{jj})|}. When the correlation matrix
#' \eqn{S} is singular (reciprocal condition number below `1e-12`) the closed
#' form does not exist and only the first deviation is returned. For a
#' network with `normalized = TRUE` the unscaled shares are used.
#'
#' @param x A relative-importance network from [relimp_network()]
#'   ([psychnet-class]).
#' @return A single number, the larger of the two deviations; 0 for an exact
#'   decomposition.
#' @section Conditions:
#'   `psychnet_bad_input` when `x` is not a network from [relimp_network()].
#' @seealso [relimp_network()]; [certificate()]; [tmfg_certificate()], the
#'   other structural check.
#' @examples
#' S <- 0.4^abs(outer(1:5, 1:5, "-"))
#' lmg_certificate(relimp_network(cor_matrix = S))
#' @export
lmg_certificate <- function(x) {
  if (!inherits(x, "psychnet") || !identical(x$method, "relimp")) {
    .psn_stop("`x` must be a relimp network from relimp_network().", call. = FALSE)
  }
  weights <- if (!is.null(x$raw_importance)) x$raw_importance else x$weights
  .lmg_residual(weights, x$r2, x$cor_matrix)
}

#' Estimate a directed relative-importance network
#'
#' @description
#' `relimp_network()` estimates a directed relative-importance network: for
#' each variable, how much of its variance the other variables explain
#' together and how that explained variance is shared among them. Each
#' variable is regressed on all the others, and its \eqn{R^2} is divided among
#' the predictors by the LMG method (Lindeman et al., 1980; Groemping, 2006):
#' a predictor's share is its increase in \eqn{R^2} when it enters the
#' regression, averaged over all orders in which the predictors can enter,
#' which is the Shapley value of the predictor.
#'
#' The shares are non-negative and sum to the \eqn{R^2}. The network is
#' directed: the edge from a predictor to an outcome carries the predictor's
#' share of the outcome's explained variance, so the two edges between a pair
#' of variables can differ. `relimp_network()` computes these shares for every
#' variable from raw data or a correlation matrix. With `normalized = TRUE`
#' each outcome's shares are rescaled to sum to 1, giving proportions of
#' explained variance rather than proportions of total variance. The
#' computation considers every subset of predictors, so its cost doubles with
#' each added variable; `max_nodes` limits the network size.
#'
#' It returns the network of shares; the identity that defines the
#' decomposition is checked by [certificate()].
#'
#' For undirected conditional associations use [ebic_glasso()] or
#' [pcor_network()]. The shares are those of
#' `relaimpo::calc.relimp(type = "lmg")` applied to each variable in turn.
#'
#' @details
#' The \eqn{R^2} of each subset regression is computed from the correlation
#' matrix. The certificate, computed by [lmg_certificate()], checks that each
#' outcome's shares sum to its full-model \eqn{R^2} and that this \eqn{R^2}
#' agrees with its closed form from the inverse correlation matrix.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#'   Not needed when `cor_matrix` is given.
#' @param cor_matrix A correlation matrix to use instead of `data`; it must be
#'   symmetric and positive semi-definite, and a covariance matrix is rescaled
#'   to a correlation matrix.
#' @param n Accepted and ignored: the decomposition needs only the
#'   correlation matrix. It exists so that a call written for the other
#'   estimators, `(cor_matrix = , n = )`, also works here.
#' @param cor_method Correlation computed from `data`: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric and polyserial
#'   correlations for ordinal items; see [cor_auto()]).
#' @param max_nodes Largest number of variables accepted, a whole number >= 2.
#'   Default `21`; larger networks are refused.
#' @param normalized Whether to rescale each outcome's incoming shares to sum
#'   to 1, as `bootnet` and `relaimpo` report them. Default `FALSE`, which keeps
#'   shares that sum to the outcome's \eqn{R^2}.
#' @param na_method Handling of missing values in `data`: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels, one unique label per node of the
#'   network.
#' @return A fitted network of class `psychnet` ([psychnet-class]), directed,
#'   whose edge weights are importance shares. `as.data.frame()` returns one
#'   row per edge with columns `from` (the predictor), `to` (the outcome) and
#'   `weight`; [certificate()] reports the structural check. The object also
#'   stores each node's full-model \eqn{R^2} (`r2`), the correlation matrix
#'   (`cor_matrix`), whether the shares were normalized (`normalized`) and,
#'   with `normalized = TRUE`, the unscaled shares (`raw_importance`), from
#'   which the certificate is computed. The sample size is recorded when the
#'   network is built from `data` and is `NA` with `cor_matrix`.
#' @section Conditions:
#'   `psychnet_bad_input` when there are more variables than `max_nodes`,
#'   `max_nodes` or `normalized` is invalid, `cor_matrix` is not a finite,
#'   symmetric, positive semi-definite matrix with a positive diagonal, or
#'   `labels` do not match the nodes; `psychnet_constant_variable` (warning)
#'   when a constant or all-missing column is dropped.
#' @seealso [lmg_certificate()] and [certificate()] for the check;
#'   [ebic_glasso()] and [pcor_network()] for undirected networks;
#'   [psychnet()], which routes `method = "relimp"` here; [psychnet-class] for
#'   the result.
#' @references Lindeman, R. H., Merenda, P. F., & Gold, R. Z. (1980).
#'   *Introduction to bivariate and multivariate analysis*. Scott, Foresman.
#'
#'   Groemping, U. (2006). Relative importance for linear regression in R: The
#'   package relaimpo. *Journal of Statistical Software*, 17(1), 1-27.
#'   \doi{10.18637/jss.v017.i01}
#' @examples
#' S <- 0.4^abs(outer(1:5, 1:5, "-"))
#' relimp_network(cor_matrix = S)
#' @export
relimp_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                           cor_method = c("pearson", "spearman", "kendall", "auto"),
                           max_nodes = 21L, normalized = FALSE,
                           na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- .psn_match_arg(cor_method)
  na_method <- .psn_match_arg(na_method)
  .psn_check_number(max_nodes, "max_nodes", lower = 2, whole = TRUE)
  .psn_check_flag(normalized, "normalized")
  # The LMG decomposition needs a covariance (PSD) matrix; name that
  # requirement in the error rather than only the generic matrix check.
  if (!is.null(cor_matrix)) {
    tryCatch(.check_cor_matrix(cor_matrix), psychnet_bad_input = function(e)
      .psn_stop("relimp_network() requires a symmetric positive-semidefinite ",
                "matrix: ", conditionMessage(e)))
  }
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "ignored")
  S <- inp$S; labels <- inp$labels
  n_obs <- if (is.null(cor_matrix)) inp$n else NA_integer_
  p <- ncol(S)
  if (p > max_nodes) {
    .psn_stop(sprintf("relimp_network() refuses p = %d > max_nodes = %d (cost grows as 2^(p-1)).",
                 p, max_nodes), call. = FALSE)
  }

  G <- matrix(0, p, p, dimnames = list(labels, labels))
  r2 <- numeric(p)
  for (j in seq_len(p)) {
    res <- .lmg_node(S, j)
    G[res$pred, j] <- res$lmg                          # predictor -> outcome
    r2[j] <- res$r2_full
  }
  diag(G) <- 0
  raw_G <- G                                           # raw LMG shares, sum to r2
  if (isTRUE(normalized)) {
    cs <- colSums(G)
    zero <- cs < sqrt(.Machine$double.eps) * max(1, max(cs))
    cs[zero] <- 1                    # ~zero columns stay zero (no noise blow-up)
    G <- sweep(G, 2L, cs, "/")
  }
  dimnames(S) <- list(labels, labels)

  # Keep the pre-normalization shares only when they differ from $weights; the
  # certificate reads them, else falls back to $weights (raw == weights here).
  extra <- list(cor_matrix = S, r2 = stats::setNames(r2, labels),
                normalized = isTRUE(normalized),
                kkt = .lmg_residual(raw_G, r2, S))
  if (isTRUE(normalized)) extra$raw_importance <- raw_G

  .new_psychnet(
    graph = G, labels = labels, method = "relimp",
    directed = TRUE, n_obs = n_obs, extra = extra
  )
}
