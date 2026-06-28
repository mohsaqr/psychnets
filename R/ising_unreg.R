# Unregularized Ising model, clean-room base R. Each binary node is regressed
# on all others by *unpenalized* logistic regression (the lambda = 0 limit of
# the penalized IRLS kernel in R/lasso_glm.R), optionally pruning edges by a
# Wald p-value, then combined by the AND/OR rule and symmetrized. This is the
# unregularized counterpart of ising_fit(); it mirrors the estimator behind
# IsingSampler / IsingFit's unregularized mode. The certificate is the GLM
# stationarity residual at lambda = 0, i.e. the maximum-likelihood score ~ 0.

# Wald p-values for one node's unpenalized logistic fit on standardized X.
# Returns p-values for the slopes (intercept excluded), in column order of X.
#' @noRd
.logit_wald_p <- function(X, b0, beta, weights = NULL) {
  eta <- b0 + as.numeric(X %*% beta)
  pr  <- 1 / (1 + exp(-eta))
  w   <- pmax(pr * (1 - pr), 1e-10)
  if (!is.null(weights)) w <- weights * w          # weighted Fisher information
  Xd  <- cbind(1, X)                               # design with intercept
  info <- crossprod(Xd, w * Xd)
  vc <- tryCatch(solve(info), error = function(e) MASS_ginv(info))
  se <- sqrt(pmax(diag(vc)[-1L], 0))               # drop intercept
  z  <- ifelse(se > 0, beta / se, 0)
  2 * stats::pnorm(-abs(z))
}

# Base-R Moore-Penrose pseudoinverse fallback (no MASS dependency).
#' @noRd
MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > max(tol * s$d[1L], 0)
  s$v[, keep, drop = FALSE] %*% ((1 / s$d[keep]) * t(s$u[, keep, drop = FALSE]))
}

#' Unregularized Ising network for binary data
#'
#' Estimates an Ising model by *unpenalized* nodewise logistic regression, with
#' optional Wald p-value edge pruning, combined by the AND (default) or OR rule.
#' The unregularized counterpart of [ising_fit()]; self-certified by the
#' maximum-likelihood score residual (see [glm_lasso_kkt()] at `lambda = 0`).
#'
#' @param data Binary (0/1) data frame or matrix (rows = observations).
#' @param rule Edge-combination rule: `"AND"` (default) or `"OR"`.
#' @param alpha Significance level for Wald edge pruning; `NULL` (default) keeps
#'   every edge.
#' @param adjust Multiple-comparison adjustment for the edge p-values (any
#'   [stats::p.adjust] method). Default `"none"`.
#' @param min_sum Minimum row sum-score; rows below it are dropped before
#'   fitting. `NULL` (default) keeps every row.
#' @param weights Optional non-negative observation weights, one per retained
#'   row. `NULL` (default) is unweighted.
#' @param na_method Missing-data handling: `"pairwise"` (default, mode-impute) or
#'   `"listwise"`. See [ising_fit()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the symmetric weight matrix,
#'   with `$thresholds` (node intercepts), `$rule`, `$p_values`, `$nodewise`
#'   (for [net_predict()]), and `$kkt` (worst nodewise score residual).
#' @examples
#' set.seed(1)
#' z <- matrix(stats::rnorm(500 * 2), 500, 2)
#' x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(500 * 4), 500)
#' b <- (x > 0) * 1L
#' colnames(b) <- paste0("V", 1:4)
#' ising_sampler(b)
#' @export
ising_sampler <- function(data, rule = c("AND", "OR"), alpha = NULL,
                          adjust = "none", min_sum = NULL, weights = NULL,
                          na_method = c("pairwise", "listwise"),
                          labels = NULL) {
  rule <- match.arg(rule)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  na_method <- match.arg(na_method)
  stopifnot(is.null(alpha) ||
              (is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1))
  mat <- .as_binary_matrix(data, na_method)
  weights <- .check_weights(weights, nrow(mat))
  ms <- .apply_min_sum(mat, weights, min_sum)
  mat <- ms$mat; weights <- ms$weights
  p <- ncol(mat)
  if (!is.null(labels)) stopifnot(length(labels) == p)
  if (is.null(labels)) labels <- colnames(mat)

  std <- .standardize(mat, weights)
  Xs_full <- std$X

  fits <- lapply(seq_len(p), function(i) {
    Xi <- Xs_full[, -i, drop = FALSE]
    y  <- mat[, i]
    fit <- .glm_lasso_fit(Xi, y, "binomial", lambda = 0, weights = weights)
    kkt <- glm_lasso_kkt(Xi, y, fit$b0, fit$beta, 0, "binomial", weights = weights)
    pv  <- .logit_wald_p(Xi, fit$b0, fit$beta, weights = weights)
    list(b0 = fit$b0, beta_std = fit$beta, beta = fit$beta / std$scale[-i],
         kkt = kkt, p = pv)
  })

  # Raw-scale edge matrix B (the Ising interaction, consistent with ising_fit)
  # with its standardized counterpart B_std for predictability and the per-node
  # certificate. The pruning p-values are scale-invariant.
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  B_std <- matrix(0, p, p, dimnames = list(labels, labels))
  P <- matrix(1, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) {
    B[i, -i] <- fits[[i]]$beta; B_std[i, -i] <- fits[[i]]$beta_std
    P[i, -i] <- fits[[i]]$p
  }
  b0_std <- vapply(fits, function(f) f$b0, numeric(1))
  worst_kkt  <- max(vapply(fits, function(f) f$kkt, numeric(1)))

  if (!is.null(alpha)) {
    Padj <- P; off <- row(P) != col(P)
    Padj[off] <- stats::p.adjust(P[off], method = adjust)
    drop <- Padj > alpha
    B[drop] <- 0; B_std[drop] <- 0
  }

  # Raw-scale node thresholds from the (post-pruning) nodewise coefficients, so
  # they stay consistent with the retained edges and with net_predict(): for
  # B[i, i] = 0, (B %*% center)[i] = sum_{j != i} beta_raw_ij * center_j.
  thresholds <- b0_std - as.numeric(B %*% std$center)

  present <- if (rule == "AND") (B != 0) & (t(B) != 0) else (B != 0) | (t(B) != 0)
  W <- (B + t(B)) / 2
  W[!present] <- 0
  diag(W) <- 0

  .new_psychnet(
    W, labels, method = "ising_sampler", directed = FALSE, n_obs = nrow(mat),
    data = mat,
    extra = list(
      thresholds = stats::setNames(thresholds, labels), rule = rule,
      p_values = P, alpha = alpha, kkt = worst_kkt,
      nodewise = list(intercept = b0_std, beta_std = B_std,
                      families = rep("binomial", p),
                      center = std$center, scale = std$scale)
    )
  )
}
