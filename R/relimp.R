# Relative-importance network by LMG / Shapley decomposition (Lindeman, Merenda
# & Gold 1980; Gromping 2006), clean-room base R. Each node is regressed on all
# others and its model R-squared is decomposed into per-predictor shares by
# averaging each predictor's marginal contribution over every subset of the
# other predictors (the Shapley value). This produces a DIRECTED network: the
# edge predictor -> outcome carries the predictor's importance share. The
# combinatorial sum over 2^(p-1) subsets is the documented exception to the
# vectorisation rule. Certified by lmg_certificate(): a node's incoming shares
# sum to its full-model R-squared (Shapley efficiency).

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

#' Relative-importance (LMG / Shapley) certificate
#'
#' By Shapley efficiency, the importance shares a node receives from its
#' predictors sum exactly to that node's full-model R-squared. Returns the
#' maximum absolute deviation from that identity; near zero certifies the
#' decomposition.
#'
#' @param x A [psychnet] object produced by [relimp_network()].
#' @return Maximum absolute deviation of incoming-share sums from the per-node
#'   R-squared (scalar); 0 = exact decomposition.
#' @examples
#' S <- 0.4^abs(outer(1:5, 1:5, "-"))
#' lmg_certificate(relimp_network(cor_matrix = S))
#' @export
lmg_certificate <- function(x) {
  max(abs(colSums(x$graph) - x$r2))
}

#' Relative-importance network (LMG / Shapley)
#'
#' Builds a directed network in which the edge predictor -> outcome is the
#' predictor's LMG (Shapley) share of the outcome node's regression R-squared.
#' Equivalent in purpose to `relaimpo::calc.relimp(type = "lmg")` applied
#' nodewise / `bootnet`'s `"relimp"` default, pure base R and self-certified via
#' [lmg_certificate()].
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix.
#' @param method Correlation method when `data` is supplied.
#' @param max_nodes Refuse to run above this many nodes (the cost grows as
#'   `2^(p-1)` per node). Default 21.
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the directed importance matrix
#'   (`graph[k, j]` = importance of `k` for outcome `j`), with `$r2` (per-node
#'   full-model R-squared), `$cor_matrix`, and `$kkt` (the decomposition
#'   residual).
#' @examples
#' S <- 0.4^abs(outer(1:5, 1:5, "-"))
#' relimp_network(cor_matrix = S)
#' @export
relimp_network <- function(data = NULL, cor_matrix = NULL,
                           method = c("pearson", "spearman", "kendall"),
                           max_nodes = 21L,
                           na_method = c("pairwise", "listwise"), labels = NULL) {
  method <- match.arg(method)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = method, na_method = na_method)
    S <- ci$S; n_obs <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- as.matrix(cor_matrix)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
    n_obs <- NA_integer_
  }
  p <- ncol(S)
  if (p > max_nodes) {
    stop(sprintf("relimp_network() refuses p = %d > max_nodes = %d (cost grows as 2^(p-1)).",
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
  dimnames(S) <- list(labels, labels)

  .new_psychnet(
    graph = G, labels = labels, method = "relimp",
    directed = TRUE, n_obs = n_obs,
    extra = list(cor_matrix = S, r2 = stats::setNames(r2, labels),
                 kkt = max(abs(colSums(G) - r2)))
  )
}
