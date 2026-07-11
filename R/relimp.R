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
  if (!inherits(x, "psychnet") || !identical(x$method, "relimp")) {
    stop("`x` must be a relimp network from relimp_network().", call. = FALSE)
  }
  weights <- if (!is.null(x$raw_importance)) x$raw_importance else x$weights
  max(abs(colSums(weights) - x$r2))
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
#' @param cor_method Correlation when `data` is supplied: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see [cor_auto()]).
#' @param max_nodes Refuse to run above this many nodes (the cost grows as
#'   `2^(p-1)` per node). Default 21.
#' @param normalized If `TRUE`, rescale each outcome's incoming importance
#'   shares to sum to 1, matching `bootnet`/`relaimpo`'s normalized reporting.
#'   The default `FALSE` keeps raw LMG shares that sum to the outcome R-squared.
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the directed importance matrix
#'   (`weights[k, j]` = importance of `k` for outcome `j`), with `$r2` (per-node
#'   full-model R-squared), `$cor_matrix`, and `$kkt` (the decomposition
#'   residual). With `normalized = FALSE` (default) each outcome's incoming
#'   shares sum to its R-squared (`colSums($weights) == $r2`); with
#'   `normalized = TRUE` they are rescaled to sum to 1, `$raw_importance` holds
#'   the unscaled shares, and `$kkt` (like [lmg_certificate()]) is computed from
#'   those raw shares.
#' @examples
#' S <- 0.4^abs(outer(1:5, 1:5, "-"))
#' relimp_network(cor_matrix = S)
#' @export
relimp_network <- function(data = NULL, cor_matrix = NULL,
                           cor_method = c("pearson", "spearman", "kendall", "auto"),
                           max_nodes = 21L, normalized = FALSE,
                           na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- match.arg(cor_method)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n_obs <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- as.matrix(cor_matrix)
    if (nrow(S) != ncol(S) || any(abs(S - t(S)) > 1e-8) ||
        any(diag(S) <= 0) ||
        min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
      stop("`cor_matrix` must be a symmetric positive-semidefinite matrix with ",
           "positive diagonal; an indefinite matrix yields R-squared shares ",
           "outside [0, 1].", call. = FALSE)
    }
    S <- stats::cov2cor(S)        # the LMG R-squared formula needs unit diagonal
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
                kkt = max(abs(colSums(raw_G) - r2)))
  if (isTRUE(normalized)) extra$raw_importance <- raw_G

  .new_psychnet(
    graph = G, labels = labels, method = "relimp",
    directed = TRUE, n_obs = n_obs, extra = extra
  )
}
