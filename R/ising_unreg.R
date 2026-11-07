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
# A singular Fisher information (separation, or an exactly collinear predictor)
# has no inverse: the Wald standard errors are then read from the
# pseudo-inverse and are not trustworthy, which is surfaced as a classed
# `psychnet_singular_information` warning naming the node rather than hidden.
.logit_wald_p <- function(X, b0, beta, weights = NULL, node = NULL) {
  eta <- b0 + as.numeric(X %*% beta)
  pr  <- 1 / (1 + exp(-eta))
  w   <- pmax(pr * (1 - pr), 1e-10)
  if (!is.null(weights)) w <- weights * w          # weighted Fisher information
  Xd  <- cbind(1, X)                               # design with intercept
  info <- crossprod(Xd, w * Xd)
  vc <- tryCatch(solve(info), error = function(e) {
    .psn_warn(sprintf(paste0(
      "ising_sampler(): the Fisher information of node %s is singular (%s); ",
      "its Wald p-values use a pseudo-inverse and are not reliable."),
      node %||% "?", conditionMessage(e)),
      class = "psychnet_singular_information")
    MASS_ginv(info)
  })
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
#' @param data Binary data frame or matrix (rows = observations): numeric 0/1,
#'   logical, or two-level factor columns, as in [ising_fit()].
#' @param rule Edge-combination rule: `"AND"` (default) or `"OR"`.
#' @param alpha Significance level for Wald edge pruning; `NULL` (default) keeps
#'   every edge.
#' @param adjust Multiple-comparison adjustment for the edge p-values (any
#'   [stats::p.adjust] method). Default `"none"`.
#' @param min_sum Minimum row sum-score; rows below it are dropped before
#'   fitting. `NULL` (default) keeps every row.
#' @param weights Optional non-negative observation weights, one per input row
#'   after missing-data preparation. When `min_sum` is used, weights are filtered
#'   alongside their corresponding rows. Rescaled to a maximum of 1, so only
#'   relative sizes matter (see [ising_fit()]). `NULL` (default) is unweighted.
#' @param na_method Missing-data handling: `"pairwise"` (default, mode-impute) or
#'   `"listwise"`. See [ising_fit()].
#' @param labels Optional node labels: one per node, or one per input column.
#' @return A `psychnet` object whose `$weights` is the symmetric weight matrix,
#'   with `$thresholds` (node intercepts), `$rule`, `$p_values`, `$nodewise`
#'   (for [net_predict()]), and `$kkt` (worst nodewise score residual).
#'   If Wald pruning changes any coefficient, `$kkt` is `NA` for the altered
#'   fit; `$fit_kkt` retains the unpruned residual and `$certificate_target`
#'   identifies it as `"pre_pruning_fit"`.
#'
#'   Under (quasi-)complete separation -- a node perfectly predicted by the
#'   others -- the unpenalized maximum-likelihood estimate does not exist: the
#'   coefficients grow without bound and the score tends to zero, so a small
#'   residual would certify nothing. The fit then raises a
#'   `psychnet_separation` warning naming the nodes, `$kkt` is `NA` and
#'   `$certificate_target` is `"nonconverged_fit"`; use the penalized
#'   [ising_fit()] instead. Any other failure to converge raises
#'   `psychnet_no_converge` with the same certificate treatment. A singular
#'   Fisher information raises `psychnet_singular_information`. Input errors
#'   are of class `psychnet_bad_input`.
#' @examples
#' b <- dichotomize(SRL_GPT)
#' ising_sampler(b)                 # every edge, certified ML fit
#' pruned <- ising_sampler(b, alpha = 0.05)
#' as.data.frame(pruned)
#' certificate(pruned)              # NA: pruning altered the ML solution
#' @export
ising_sampler <- function(data, rule = c("AND", "OR"), alpha = NULL,
                          adjust = "none", min_sum = NULL, weights = NULL,
                          na_method = c("pairwise", "listwise"),
                          labels = NULL) {
  rule <- .psn_match_arg(rule)
  adjust <- .psn_match_arg(adjust, stats::p.adjust.methods)
  na_method <- .psn_match_arg(na_method)
  .psn_check_number(alpha, "alpha", lower = 0, upper = 1, open = TRUE,
                    null_ok = TRUE)
  .psn_check_number(min_sum, "min_sum", lower = 0, null_ok = TRUE)
  bin <- .as_binary_matrix(data, na_method, caller = "ising_sampler")
  mat <- bin$mat
  weights <- .check_weights(weights, nrow(mat))
  ms <- .apply_min_sum(mat, weights, min_sum)
  mat <- ms$mat; weights <- ms$weights
  p <- ncol(mat)
  labels <- .remap_labels(labels, bin$source_names, colnames(mat)) %||%
    colnames(mat)

  std <- .standardize(mat, weights)
  Xs_full <- std$X

  fits <- lapply(seq_len(p), function(i) {
    Xi <- Xs_full[, -i, drop = FALSE]
    y  <- mat[, i]
    fit <- .glm_lasso_fit(Xi, y, "binomial", lambda = 0, weights = weights)
    kkt <- glm_lasso_kkt(Xi, y, fit$b0, fit$beta, 0, "binomial", weights = weights)
    pv  <- .logit_wald_p(Xi, fit$b0, fit$beta, weights = weights,
                         node = labels[i])
    list(b0 = fit$b0, beta_std = fit$beta, beta = fit$beta / std$scale[-i],
         kkt = kkt, p = pv, converged = fit$converged,
         separation = fit$separation)
  })
  # A node whose IRLS never converged has no maximum-likelihood estimate to
  # certify: under separation the score -> 0 while |beta| -> Inf, so the
  # residual is small however wrong the coefficients are.
  separated <- vapply(fits, function(f) isTRUE(f$separation), logical(1))
  unconverged <- !vapply(fits, function(f) isTRUE(f$converged), logical(1))
  if (any(separated)) {
    .psn_warn(sprintf(paste0(
      "ising_sampler(): node(s) %s are (quasi-)completely separated -- ",
      "perfectly predicted by the other nodes -- so their unpenalized ",
      "maximum-likelihood estimate does not exist and the reported weights ",
      "are arbitrary. Use the penalized ising_fit(), or remove the redundant ",
      "node."), paste(labels[separated], collapse = ", ")),
      class = "psychnet_separation")
  } else if (any(unconverged)) {
    .psn_warn(sprintf(
      "ising_sampler(): the logistic fit(s) for %s did not converge; the returned weights are not certified.",
      paste(labels[unconverged], collapse = ", ")),
      class = "psychnet_no_converge")
  }
  nonconverged <- any(unconverged)

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
  pruned <- FALSE

  if (!is.null(alpha)) {
    Padj <- P; off <- row(P) != col(P)
    Padj[off] <- stats::p.adjust(P[off], method = adjust)
    drop <- Padj > alpha
    pruned <- any(B[drop] != 0)
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
      p_values = P, alpha = alpha,
      kkt = if (pruned || nonconverged) NA_real_ else worst_kkt,
      fit_kkt = worst_kkt,
      certificate_target = if (nonconverged) "nonconverged_fit" else
        if (pruned) "pre_pruning_fit" else "returned_network",
      nodewise = list(intercept = b0_std, beta_std = B_std,
                      families = rep("binomial", p),
                      center = std$center, scale = std$scale,
                      category_levels = bin$category_levels)
    )
  )
}
