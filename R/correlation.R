# Unregularized correlation and partial-correlation networks (pure base R).

# Coerce a data frame / matrix to a clean numeric matrix: keep numeric columns
# with non-zero variance, guarantee column names. `drop_na = TRUE` removes
# incomplete rows (listwise); FALSE keeps them (for pairwise handling upstream).
#' @noRd
.as_numeric_matrix <- function(data, drop_na = TRUE) {
  if (is.null(data)) stop("`data` is required.", call. = FALSE)
  if (is.matrix(data)) {
    mat <- data
    storage.mode(mat) <- "double"
  } else if (is.data.frame(data)) {
    num <- vapply(data, is.numeric, logical(1))
    if (!any(num)) stop("No numeric columns in `data`.", call. = FALSE)
    mat <- as.matrix(data[, num, drop = FALSE])
  } else {
    stop("`data` must be a data frame or numeric matrix.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  if (drop_na) mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  vars <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  if (any(vars == 0 | is.na(vars))) {
    mat <- mat[, vars > 0 & !is.na(vars), drop = FALSE]
  }
  if (ncol(mat) < 2L) stop("Need at least 2 usable variables.", call. = FALSE)
  mat
}

# Correlation matrix + effective sample size from possibly-incomplete data.
# `na_method = "listwise"` drops rows with any NA (the classic default, but it
# collapses catastrophically when missingness is spread across many columns).
# `na_method = "pairwise"` uses pairwise-complete correlations projected to the
# nearest positive-definite matrix, retaining far more information under MCAR/MAR
# (matching `qgraph::cor_auto` / `NetworkToolbox`'s `na.data`). With no missing
# data the two are identical. Returns the correlation, the effective n, and the
# usable labels.
#' @noRd
.cor_input <- function(data, method = "pearson",
                       na_method = c("pairwise", "listwise")) {
  na_method <- match.arg(na_method)
  mat <- .as_numeric_matrix(data, drop_na = FALSE)
  if (!anyNA(mat) || na_method == "listwise") {
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    if (nrow(mat) < 2L) stop("Need at least 2 complete observations.", call. = FALSE)
    return(list(S = stats::cor(mat, method = method), n = nrow(mat),
                labels = colnames(mat), na_method = "listwise"))
  }
  S <- stats::cor(mat, use = "pairwise.complete.obs", method = method)
  if (anyNA(S)) {
    stop("Pairwise correlation is undefined: a variable pair is never co-observed.",
         call. = FALSE)
  }
  S <- .nearest_pd_cor(S)                         # pairwise S may not be PD
  co <- crossprod(!is.na(mat))                    # pairwise co-observation counts
  n_eff <- max(round(stats::median(co[upper.tri(co)])), 2L)
  list(S = S, n = n_eff, labels = colnames(mat), na_method = "pairwise")
}

# Missing-data prep for the nodewise (Ising / mgm) estimators, which need a
# complete design. `na_method = "listwise"` drops incomplete rows; `"pairwise"`
# keeps every row and single-imputes each column over its observed values --
# the mode for binary columns, the mean otherwise -- retaining the full sample.
#' @noRd
.na_prep_nodewise <- function(mat, na_method = c("pairwise", "listwise")) {
  na_method <- match.arg(na_method)
  if (!anyNA(mat)) return(mat)
  if (na_method == "listwise") {
    return(mat[stats::complete.cases(mat), , drop = FALSE])
  }
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]; miss <- is.na(x)
    if (!any(miss)) next
    obs <- x[!miss]; u <- unique(obs)
    fill <- if (all(u %in% c(0, 1))) as.numeric(mean(obs) >= 0.5) else mean(obs)
    mat[miss, j] <- fill
  }
  mat
}

# Two-sided p-values for a (partial) correlation matrix. `k` is the number of
# variables partialled out (0 for marginal correlations, p - 2 for full-order
# partial correlations); df = n - 2 - k.
#' @noRd
.cor_pvalues <- function(r, n, k) {
  df <- max(n - 2L - k, 1L)
  tstat <- r * sqrt(df / pmax(1 - r^2, 1e-12))
  P <- 2 * stats::pt(-abs(tstat), df)
  diag(P) <- 1
  P
}

# Zero edges whose adjusted p-value exceeds alpha (symmetric upper-tri adjust).
#' @noRd
.apply_sig <- function(g, P, alpha, adjust) {
  ut <- upper.tri(P)
  padj <- P
  padj[ut] <- stats::p.adjust(P[ut], method = adjust)
  padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]
  g[padj > alpha] <- 0
  g
}

#' Correlation network
#'
#' Marginal (zero-order) association network: the Pearson correlation matrix
#' with the diagonal removed. Equivalent to `bootnet`'s `"cor"` default.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Correlation method: `"pearson"` (default), `"spearman"`, or
#'   `"kendall"`.
#' @param threshold Correlations with absolute value below this are set to zero.
#'   Default 0.
#' @param alpha Significance level; if set, correlations not significant at
#'   `alpha` are zeroed. `NULL` (default) keeps every edge.
#' @param adjust Multiple-comparison adjustment for the edge p-values (any
#'   [stats::p.adjust] method). Default `"none"`.
#' @param na_method Missing-data handling: `"pairwise"` (default) uses
#'   pairwise-complete correlations projected to the nearest positive-definite
#'   matrix; `"listwise"` drops rows with any `NA`. Identical when data is
#'   complete.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$graph` is the thresholded correlation
#'   matrix, with `$cor_matrix`, `$n_eff`, `$na_method` (and `$p_values` when
#'   `alpha` is used).
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' cor_network(x)
#' cor_network(x, alpha = 0.05, adjust = "BH")
#' @export
cor_network <- function(data, method = c("pearson", "spearman", "kendall"),
                        threshold = 0, alpha = NULL, adjust = "none",
                        na_method = c("pairwise", "listwise"), labels = NULL) {
  method <- match.arg(method)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  na_method <- match.arg(na_method)
  ci <- .cor_input(data, method = method, na_method = na_method)
  S <- ci$S; n <- ci$n
  if (is.null(labels)) labels <- ci$labels
  g <- S
  diag(g) <- 0
  g[abs(g) < threshold] <- 0
  extra <- list(cor_matrix = S, n_eff = n, na_method = ci$na_method)
  if (!is.null(alpha)) {
    P <- .cor_pvalues(g, n, k = 0L)
    g <- .apply_sig(g, P, alpha, adjust)
    extra$p_values <- P
  }
  .new_psychnet(g, labels, method = "cor", directed = FALSE,
                n_obs = n, extra = extra)
}

#' Partial correlation network
#'
#' Conditional (full-order) association network: each edge is the correlation
#' between two variables with all others partialled out, obtained from the
#' inverse correlation matrix. Equivalent to `bootnet`'s `"pcor"` default.
#'
#' @inheritParams cor_network
#' @return A `psychnet` object whose `$graph` is the thresholded
#'   partial-correlation matrix, with `$precision`, `$cor_matrix` (and
#'   `$p_values` when `alpha` is used).
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' pcor_network(x)
#' pcor_network(x, alpha = 0.05, adjust = "holm")
#' @export
pcor_network <- function(data, method = c("pearson", "spearman", "kendall"),
                         threshold = 0, alpha = NULL, adjust = "none",
                         na_method = c("pairwise", "listwise"), labels = NULL) {
  method <- match.arg(method)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  na_method <- match.arg(na_method)
  ci <- .cor_input(data, method = method, na_method = na_method)
  S <- ci$S; n <- ci$n
  if (is.null(labels)) labels <- ci$labels
  wi <- solve(S)
  g  <- .precision_to_pcor(wi)
  g[abs(g) < threshold] <- 0
  extra <- list(precision = wi, cor_matrix = S, n_eff = n, na_method = ci$na_method)
  if (!is.null(alpha)) {
    P <- .cor_pvalues(g, n, k = ncol(S) - 2L)   # full-order partials
    g <- .apply_sig(g, P, alpha, adjust)
    extra$p_values <- P
  }
  dimnames(g) <- list(labels, labels)
  .new_psychnet(g, labels, method = "pcor", directed = FALSE,
                n_obs = n, extra = extra)
}
