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

# Validate and normalize a user-supplied `cor_matrix`: square, at least two
# variables, finite, symmetric, positive diagonal, and positive semi-definite.
# A covariance matrix is normalized to unit diagonal (`cov2cor`) so the GGM
# estimators are scale-invariant on user input -- a proper correlation matrix
# passes through byte-unchanged (its diagonal is already 1). Malformed input is
# rejected early with a clear message instead of a cryptic downstream
# eigen/solve failure, a `max(numeric(0))` warning, or an all-NaN network.
#' @noRd
.check_cor_matrix <- function(cor_matrix) {
  S <- as.matrix(cor_matrix)
  if (nrow(S) != ncol(S)) {
    stop("`cor_matrix` must be a square matrix.", call. = FALSE)
  }
  if (ncol(S) < 2L) {
    stop("`cor_matrix` must have at least 2 variables.", call. = FALSE)
  }
  if (any(!is.finite(S))) {
    stop("`cor_matrix` must not contain missing or infinite values.",
         call. = FALSE)
  }
  if (any(abs(S - t(S)) > 1e-8)) {
    stop("`cor_matrix` must be symmetric.", call. = FALSE)
  }
  if (any(diag(S) <= 0)) {
    stop("`cor_matrix` must have a strictly positive diagonal.", call. = FALSE)
  }
  min_eig <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < -1e-8) {
    stop(sprintf(paste0("`cor_matrix` is not positive semi-definite (min ",
                        "eigenvalue %.2e); supply a PD matrix or pass raw `data`."),
                 min_eig), call. = FALSE)
  }
  stats::cov2cor(S)   # no-op for a unit-diagonal correlation matrix
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
                       na_method = c("pairwise", "listwise"),
                       ordinal_max_levels = 7L) {
  na_method <- match.arg(na_method)
  mat <- .as_numeric_matrix(data, drop_na = FALSE)
  if (method == "auto") {                           # polychoric / polyserial
    if (na_method == "listwise") {
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    }
    if (nrow(mat) < 2L) stop("Need at least 2 observations.", call. = FALSE)
    has_na <- anyNA(mat)
    S <- .cor_auto_matrix(mat, ordinal_max_levels)
    n <- if (has_na) {
      co <- crossprod(!is.na(mat)); max(round(stats::median(co[upper.tri(co)])), 2L)
    } else nrow(mat)
    return(list(S = S, n = n, labels = colnames(mat),
                na_method = if (has_na) "pairwise" else "listwise"))
  }
  if (!anyNA(mat) || na_method == "listwise") {
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    if (nrow(mat) < 2L) stop("Need at least 2 complete observations.", call. = FALSE)
    # Row deletion can turn a variable constant; that would give NA correlations.
    if (any(apply(mat, 2L, stats::sd) == 0)) {
      stop("After listwise deletion a variable has zero variance.", call. = FALSE)
    }
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

# Invert a correlation matrix, projecting to the nearest PD matrix first if it
# is singular (e.g. a complete-data correlation with n < p).
#' @noRd
.pd_solve <- function(S) {
  out <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(out)) out <- solve(.nearest_pd_cor(S))
  out
}

# Two-sided p-values for a (partial) correlation matrix. `k` is the number of
# variables partialled out (0 for marginal correlations, p - 2 for full-order
# partial correlations); df = n - 2 - k.
#' @noRd
.cor_pvalues <- function(r, n, k) {
  df <- n - 2L - k
  if (df <= 0L) {                       # no residual df: the test is undefined,
    P <- matrix(1, nrow(r), ncol(r))    # so nothing can be called significant
    return(P)
  }
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
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional precomputed correlation matrix; if given, `data`
#'   is ignored and `n` is required when `alpha` is used.
#' @param n Sample size (needed for significance testing when `cor_matrix` is
#'   supplied).
#' @param cor_method Correlation method: `"pearson"` (default), `"spearman"`,
#'   `"kendall"`, or `"auto"` (polychoric/polyserial for ordinal items, the
#'   `qgraph::cor_auto` default; see [cor_auto()]).
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
#' @return A `psychnet` object whose `$weights` is the thresholded correlation
#'   matrix, with `$cor_matrix`, `$n_eff`, `$na_method` (and `$p_values` when
#'   `alpha` is used).
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' cor_network(x)
#' cor_network(x, alpha = 0.05, adjust = "BH")
#' @export
cor_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                        cor_method = c("pearson", "spearman", "kendall", "auto"),
                        threshold = 0, alpha = NULL, adjust = "none",
                        na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- match.arg(cor_method)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n <- ci$n; na_used <- ci$na_method
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (!is.null(alpha) && is.null(n)) {
      stop("`n` is required when `cor_matrix` is supplied and `alpha` is set.",
           call. = FALSE)
    }
    na_used <- "none"
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
  }
  g <- S
  diag(g) <- 0
  r_full <- g                                  # p-values use the true correlations
  g[abs(g) < threshold] <- 0
  extra <- list(cor_matrix = S, n_eff = n, na_method = na_used)
  if (!is.null(alpha)) {
    P <- .cor_pvalues(r_full, n, k = 0L)
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
#' @return A `psychnet` object whose `$weights` is the thresholded
#'   partial-correlation matrix, with `$precision`, `$cor_matrix` (and
#'   `$p_values` when `alpha` is used).
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' pcor_network(x)
#' pcor_network(x, alpha = 0.05, adjust = "holm")
#' @export
pcor_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         threshold = 0, alpha = NULL, adjust = "none",
                         na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- match.arg(cor_method)
  adjust <- match.arg(adjust, stats::p.adjust.methods)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n <- ci$n; na_used <- ci$na_method
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (!is.null(alpha) && is.null(n)) {
      stop("`n` is required when `cor_matrix` is supplied and `alpha` is set.",
           call. = FALSE)
    }
    na_used <- "none"
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
  }
  wi <- .pd_solve(S)
  g  <- .precision_to_pcor(wi)
  r_full <- g                                   # p-values use the true partials
  g[abs(g) < threshold] <- 0
  extra <- list(precision = wi, cor_matrix = S, n_eff = n, na_method = na_used)
  if (!is.null(alpha)) {
    P <- .cor_pvalues(r_full, n, k = ncol(S) - 2L)   # full-order partials
    g <- .apply_sig(g, P, alpha, adjust)
    extra$p_values <- P
  }
  dimnames(g) <- list(labels, labels)
  .new_psychnet(g, labels, method = "pcor", directed = FALSE,
                n_obs = n, extra = extra)
}
