# Mixed graphical model (Haslbeck & Waldorp 2020), clean-room base R. Each node
# is regressed on all others with the L1-penalized GLM matching its type --
# gaussian for continuous nodes, logistic for binary nodes -- and the asymmetric
# standardized estimates are combined by the AND rule and symmetrized. v0.1
# supports gaussian and binary nodes (the dominant mixed case); categorical
# nodes with more than two levels are not yet implemented and error explicitly.

# Detect node types: "c" (binary categorical, values in {0,1}) or "g" (gaussian).
#' @noRd
.detect_types <- function(mat) {
  vapply(seq_len(ncol(mat)), function(j) {
    u <- unique(mat[, j])
    if (length(u) == 2L && all(u %in% c(0, 1))) {
      "c"
    } else if (length(u) <= 10L && all(u == round(u)) && any(!u %in% c(0, 1))) {
      stop(sprintf(
        "Column '%s' looks categorical with >2 levels; mgm_fit() v0.1 supports only gaussian and binary nodes. One-hot encode it first.",
        colnames(mat)[j]), call. = FALSE)
    } else {
      "g"
    }
  }, character(1))
}

#' Mixed graphical model
#'
#' Estimates a mixed graphical model by nodewise L1-penalized regression -- a
#' gaussian (linear) lasso for continuous nodes and a logistic lasso for binary
#' nodes -- with per-node EBIC selection, combined by the AND rule. Equivalent
#' in purpose to `mgm::mgm()`, but pure base R and self-certified: each node's
#' regression reports its stationarity (KKT) residual (see [glm_lasso_kkt()]).
#'
#' @param data Numeric data frame or matrix (rows = observations); columns are
#'   continuous or binary (0/1).
#' @param gamma EBIC hyperparameter. Default 0.25.
#' @param types Optional character vector of node types (`"g"` gaussian, `"c"`
#'   binary); auto-detected if `NULL`.
#' @param nlambda Number of penalties per nodewise path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the symmetric standardized
#'   weight matrix, with `$types` and `$kkt` (the worst nodewise residual).
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(400)
#' g1 <- f + stats::rnorm(400); g2 <- f + stats::rnorm(400)
#' b1 <- (f + stats::rnorm(400) > 0) * 1L
#' d <- data.frame(g1 = g1, g2 = g2, b1 = b1, n = stats::rnorm(400))
#' mgm_fit(d)
#' @export
mgm_fit <- function(data, gamma = 0.25, types = NULL,
                    nlambda = 100L, lambda_min_ratio = 0.01, labels = NULL) {
  mat <- .as_numeric_matrix(data)
  p <- ncol(mat)
  if (is.null(labels)) labels <- colnames(mat)
  if (is.null(types))  types  <- .detect_types(mat)
  stopifnot(length(types) == p, all(types %in% c("g", "c")))

  fits <- lapply(seq_len(p), function(i) {
    fam <- if (types[i] == "c") "binomial" else "gaussian"
    .nodewise_ebic(mat[, -i, drop = FALSE], mat[, i], fam,
                   gamma, nlambda, lambda_min_ratio)
  })

  # Standardized asymmetric weights make cross-family edges comparable.
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) B[i, -i] <- fits[[i]]$beta_std
  intercepts <- vapply(fits, function(f) f$b0, numeric(1))
  worst_kkt <- max(vapply(fits, function(f) f$kkt, numeric(1)))
  std <- .standardize(mat)
  families <- ifelse(types == "c", "binomial", "gaussian")

  present <- (B != 0) & (t(B) != 0)            # AND rule
  W <- (B + t(B)) / 2
  W[!present] <- 0
  diag(W) <- 0

  .new_psychnet(W, labels, method = "mgm", directed = FALSE,
                n_obs = nrow(mat),
                extra = list(types = stats::setNames(types, labels),
                             kkt = worst_kkt,
                             nodewise = list(intercept = intercepts,
                                             beta_std = B,
                                             families = families,
                                             center = std$center,
                                             scale = std$scale)))
}
