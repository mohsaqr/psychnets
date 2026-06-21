# Ising model for binary data (van Borkulo et al. 2014), clean-room base R.
# Each node is regressed on all others by L1-penalized logistic regression
# (R/lasso_glm.R) with per-node EBIC selection; the asymmetric nodewise
# estimates are combined by the AND or OR rule and symmetrized.

# Coerce + validate a binary (0/1) matrix, then apply NA handling.
#' @noRd
.as_binary_matrix <- function(data, na_method = c("pairwise", "listwise")) {
  na_method <- match.arg(na_method)
  mat <- .as_numeric_matrix(data, drop_na = FALSE)
  u <- unique(as.vector(mat)); u <- u[!is.na(u)]
  if (!all(u %in% c(0, 1))) {
    stop("ising_fit() requires binary (0/1) data.", call. = FALSE)
  }
  .na_prep_nodewise(mat, na_method)
}

#' Ising network for binary data
#'
#' Estimates an Ising model by nodewise L1-penalized logistic regression with
#' EBIC selection, combined by the AND (default) or OR rule. Equivalent in
#' purpose to `IsingFit::IsingFit()`, but pure base R and self-certified: each
#' node's regression reports its stationarity (KKT) residual (see
#' [glm_lasso_kkt()]).
#'
#' @param data Binary (0/1) data frame or matrix (rows = observations).
#' @param gamma EBIC hyperparameter. Default 0.25.
#' @param rule Edge-combination rule: `"AND"` (default) or `"OR"`.
#' @param nlambda Number of penalties per nodewise path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param na_method Missing-data handling: `"pairwise"` (default) single-imputes
#'   each column over its observed values (mode for binary), keeping the full
#'   sample; `"listwise"` drops incomplete rows. Identical for complete data.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the symmetric weight matrix,
#'   with `$thresholds` (node intercepts) and `$kkt` (the worst nodewise
#'   stationarity residual).
#' @examples
#' set.seed(1)
#' z <- matrix(stats::rnorm(400 * 2), 400, 2)
#' x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(400 * 4), 400)
#' b <- (x > 0) * 1L
#' colnames(b) <- paste0("V", 1:4)
#' ising_fit(b)
#' @export
ising_fit <- function(data, gamma = 0.25, rule = c("AND", "OR"),
                      nlambda = 100L, lambda_min_ratio = 0.01,
                      na_method = c("pairwise", "listwise"), labels = NULL) {
  rule <- match.arg(rule)
  na_method <- match.arg(na_method)
  mat <- .as_binary_matrix(data, na_method)
  p <- ncol(mat)
  if (is.null(labels)) labels <- colnames(mat)

  fits <- lapply(seq_len(p), function(i) {
    .nodewise_ebic(mat[, -i, drop = FALSE], mat[, i], "binomial",
                   gamma, nlambda, lambda_min_ratio)
  })

  # Asymmetric matrix B: B[i, j] = effect of node j on node i on the raw 0/1
  # scale (the Ising interaction), with B_std its standardized counterpart.
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  B_std <- matrix(0, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) { B[i, -i] <- fits[[i]]$beta; B_std[i, -i] <- fits[[i]]$beta_std }
  std <- .standardize(mat)
  b0_std <- vapply(fits, function(f) f$b0, numeric(1))
  # Node thresholds on the raw scale: the standardized intercept de-centered by
  # the raw slopes, tau_i = b0_std_i - sum_j beta_raw_ij * center_j, so they are
  # coherent with the raw-scale edge weights in W.
  thresholds <- vapply(seq_len(p), function(i)
    b0_std[i] - sum(fits[[i]]$beta * std$center[-i]), numeric(1))
  worst_kkt  <- max(vapply(fits, function(f) f$kkt, numeric(1)))

  present <- if (rule == "AND") (B != 0) & (t(B) != 0) else (B != 0) | (t(B) != 0)
  W <- (B + t(B)) / 2
  W[!present] <- 0
  diag(W) <- 0

  .new_psychnet(W, labels, method = "ising", directed = FALSE,
                n_obs = nrow(mat),
                extra = list(thresholds = stats::setNames(thresholds, labels),
                             rule = rule, kkt = worst_kkt,
                             nodewise = list(intercept = b0_std,
                                             beta_std = B_std,
                                             families = rep("binomial", p),
                                             center = std$center,
                                             scale = std$scale)))
}
