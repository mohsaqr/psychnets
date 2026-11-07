# Unregularized correlation and partial-correlation networks (pure base R).

# Package-private state. `quiet_columns` mutes the column-drop notices of
# .as_numeric_matrix() while a resampling verb refits a draw: the notice is
# raised once on the caller's full data, and repeating it per draw (a column
# that is constant in one bootstrap sample, a factor column a GGM ignores)
# would bury it under hundreds of copies. Set it only through
# .psn_with_quiet_columns(), which restores the previous value on exit.
.psn_state <- new.env(parent = emptyenv())
.psn_state$quiet_columns <- FALSE

# Evaluate `expr` with the column-drop notices of .as_numeric_matrix() muted.
#' @noRd
.psn_with_quiet_columns <- function(expr) {
  old <- .psn_state$quiet_columns
  .psn_state$quiet_columns <- TRUE
  on.exit(.psn_state$quiet_columns <- old, add = TRUE, after = FALSE)
  expr
}

# Classed message naming the data-frame columns an estimator cannot use.
#' @noRd
.psn_note_dropped <- function(cols) {
  if (!length(cols) || isTRUE(.psn_state$quiet_columns)) return(invisible(NULL))
  message(structure(
    class = c("psychnet_dropped_columns", "message", "condition"),
    list(message = sprintf(paste0(
      "Dropped %d non-numeric column(s): %s. Code an ordinal item as an ",
      "ordered factor or integer to keep it; one-hot encode a nominal one.\n"),
      length(cols), paste(cols, collapse = ", ")), call = NULL)))
}

# Classed warning naming the columns removed for having no variance.
#' @noRd
.psn_warn_constant <- function(cols, after_listwise = FALSE) {
  if (!length(cols) || isTRUE(.psn_state$quiet_columns)) return(invisible(NULL))
  .psn_warn(sprintf(paste0(
    "Dropped %d constant or all-missing column(s)%s: %s. They carry no ",
    "association and are not nodes of the network."),
    length(cols), if (after_listwise) " (after listwise deletion)" else "",
    paste(cols, collapse = ", ")), class = "psychnet_constant_variable")
}

# Coerce a data frame / matrix to a clean numeric matrix: keep numeric columns
# with non-zero variance, guarantee column names. `drop_na = TRUE` removes
# incomplete rows (listwise); FALSE keeps them (for pairwise handling upstream).
# An ordered factor is an ordinal item and becomes its integer level code (so
# cor_method = "auto" treats it as ordinal, as qgraph::cor_auto does). Other
# non-numeric data-frame columns (unordered factor, character, logical) are
# dropped with a `psychnet_dropped_columns` message; constant / all-missing
# columns are dropped with a `psychnet_constant_variable` warning. Numeric
# input is returned byte-identical to before -- only the notices are new.
#' @noRd
.as_numeric_matrix <- function(data, drop_na = TRUE) {
  if (is.null(data)) .psn_stop("`data` is required.", call. = FALSE)
  if (is.matrix(data)) {
    mat <- data
    storage.mode(mat) <- "double"
  } else if (is.data.frame(data)) {
    ordered_cols <- vapply(data, is.ordered, logical(1))
    data[ordered_cols] <- lapply(data[ordered_cols], as.integer)
    num <- vapply(data, is.numeric, logical(1))
    if (!any(num)) .psn_stop("No numeric columns in `data`.", call. = FALSE)
    .psn_note_dropped(names(data)[!num])
    mat <- as.matrix(data[, num, drop = FALSE])
  } else {
    .psn_stop("`data` must be a data frame or numeric matrix.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  n_rows <- nrow(mat)
  if (drop_na) mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  vars <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  usable <- vars > 0 & !is.na(vars)
  if (!all(usable)) {
    .psn_warn_constant(colnames(mat)[!usable],
                       after_listwise = nrow(mat) < n_rows)
    mat <- mat[, usable, drop = FALSE]
  }
  if (ncol(mat) < 2L) .psn_stop("Need at least 2 usable variables.", call. = FALSE)
  mat
}

# The rows a resampling verb (net_boot, net_stability, net_casedrop_reliability,
# net_split_reliability) draws from, and the estimator dots it forwards.
#
# `data` is kept as the caller gave it: a data.frame keeps its factor columns
# and its incomplete rows; a matrix is coerced once WITHOUT dropping rows. Every
# draw is therefore handed to the estimator exactly as the full data would be,
# and the estimator applies its own `na_method` and type handling per draw --
# so resampling `mgm_fit(native = FALSE)` keeps a factor node, resampling
# glasso drops it, and `na_method = "pairwise"` (the default) resamples the
# same rows a direct call would use instead of a listwise-reduced subset.
# Columns with fewer than two distinct observed values are removed once, as
# every estimator would remove them.
#
# When the estimator is mgm and the caller gave no `types`, node types are
# decided ONCE here, on the full data, and pinned for every fit: a draw could
# otherwise flip a borderline column (11 distinct integers in full, <= 10 in
# the draw) from gaussian to categorical, silently changing the model between
# draws. Deciding once also surfaces mgm's type-promotion warning exactly once
# per verb call rather than once per draw.
#' @noRd
.resample_input <- function(data, method, dots) {
  if (is.data.frame(data)) {
    keep <- .usable_columns(data)
    .psn_warn_constant(names(data)[!keep])
    data <- data[keep]
    if (ncol(data) < 2L) .psn_stop("Need at least 2 usable variables.", call. = FALSE)
  } else {
    data <- .as_numeric_matrix(data, drop_na = FALSE)
  }
  if (identical(.resolve_method(method), "mgm") && is.null(dots$types)) {
    dots$types <- .detect_types(data)$type
  }
  list(data = data, dots = dots)
}

# Merge the dedicated estimator-argument list used by resampling verbs. The
# separate channel is necessary because those verbs use names such as
# `threshold` for their own diagnostics. Explicit estimator_args win over `...`.
#' @noRd
.resample_dots <- function(dots, estimator_args = list()) {
  if (!is.list(estimator_args) || is.null(names(estimator_args)) && length(estimator_args)) {
    .psn_stop("`estimator_args` must be a named list.", call. = FALSE)
  }
  if (length(estimator_args)) dots[names(estimator_args)] <- estimator_args
  dots
}

# Back-compatible translation of net_boot(engine=) to the estimators' public
# native switch. Keeping this in one place prevents an internal `engine`
# argument from leaking through psychnet() as an unused argument.
#' @noRd
.resample_engine <- function(engine, method) {
  if (is.null(engine)) return(NULL)
  if (!is.character(engine) || length(engine) != 1L || is.na(engine))
    .psn_stop("`engine` must be NULL or one of 'base', 'glasso', or 'glmnet'.",
         call. = FALSE)
  resolved <- .resolve_method(method)
  allowed <- if (resolved %in% c("glasso", "huge")) c("base", "glasso")
             else if (resolved %in% c("ising", "mgm")) c("base", "glmnet")
             else character()
  if (!engine %in% allowed) {
    .psn_stop(sprintf("engine = '%s' is not available for method = '%s'.",
                 engine, resolved), call. = FALSE)
  }
  identical(engine, "base")
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
    .psn_stop("`cor_matrix` must be a square matrix.", call. = FALSE)
  }
  if (ncol(S) < 2L) {
    .psn_stop("`cor_matrix` must have at least 2 variables.", call. = FALSE)
  }
  if (any(!is.finite(S))) {
    .psn_stop("`cor_matrix` must not contain missing or infinite values.",
         call. = FALSE)
  }
  if (any(diag(S) <= 0)) {
    .psn_stop("`cor_matrix` must have a strictly positive diagonal.", call. = FALSE)
  }
  # Judge tolerances on correlation units: an indefinite covariance must not
  # pass merely because its measurement units make every entry tiny.
  S <- stats::cov2cor(S)
  if (any(abs(S - t(S)) > 1e-8)) {
    .psn_stop("`cor_matrix` must be symmetric.", call. = FALSE)
  }
  min_eig <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < -1e-8) {
    .psn_stop(sprintf(paste0("`cor_matrix` is not positive semi-definite (min ",
                        "eigenvalue %.2e); supply a PD matrix or pass raw `data`."),
                 min_eig), call. = FALSE)
  }
  S
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
  na_method <- .psn_match_arg(na_method)
  mat <- .as_numeric_matrix(data, drop_na = FALSE)
  if (method == "auto") {                           # polychoric / polyserial
    if (na_method == "listwise") {
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    }
    if (nrow(mat) < 2L) .psn_stop("Need at least 2 observations.", call. = FALSE)
    has_na <- anyNA(mat)
    S <- .cor_auto_matrix(mat, ordinal_max_levels)
    n <- if (has_na) {
      co <- crossprod(!is.na(mat)); max(round(stats::median(co[upper.tri(co)])), 2L)
    } else nrow(mat)
    co <- crossprod(!is.na(mat)); dimnames(co) <- list(colnames(mat), colnames(mat))
    return(list(S = S, n = n, n_pair = co, labels = colnames(mat),
                na_method = if (has_na) "pairwise" else "listwise"))
  }
  if (!anyNA(mat) || na_method == "listwise") {
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    if (nrow(mat) < 2L) .psn_stop("Need at least 2 complete observations.", call. = FALSE)
    # Row deletion can turn a variable constant; that would give NA correlations.
    if (any(apply(mat, 2L, stats::sd) == 0)) {
      .psn_stop("After listwise deletion a variable has zero variance.", call. = FALSE)
    }
    npm <- matrix(nrow(mat), ncol(mat), ncol(mat),
                  dimnames = list(colnames(mat), colnames(mat)))
    return(list(S = stats::cor(mat, method = method), n = nrow(mat),
                n_pair = npm,
                labels = colnames(mat), na_method = "listwise"))
  }
  S <- stats::cor(mat, use = "pairwise.complete.obs", method = method)
  if (anyNA(S)) {
    .psn_stop("Pairwise correlation is undefined: a variable pair is never co-observed.",
         call. = FALSE)
  }
  S <- .nearest_pd_cor(S)                         # pairwise S may not be PD
  co <- crossprod(!is.na(mat))                    # pairwise co-observation counts
  dimnames(co) <- list(colnames(mat), colnames(mat))
  n_eff <- max(round(stats::median(co[upper.tri(co)])), 2L)
  list(S = S, n = n_eff, n_pair = co, labels = colnames(mat),
       na_method = "pairwise")
}

# Missing-data prep for the nodewise (Ising / mgm) estimators, which need a
# complete design. `na_method = "listwise"` drops incomplete rows; `"pairwise"`
# keeps every row and single-imputes each column over its observed values --
# the mode for binary columns, the mean otherwise -- retaining the full sample.
#' @noRd
# `categorical` (optional logical, one per column) marks columns whose values are
# level codes: a missing entry there is filled with the modal code, never a mean
# that would fall between levels. A 0/1 column keeps its historical rule (mean
# rounded, ties to 1) whether or not it is flagged, so existing binary fits are
# byte-identical.
#' @noRd
.na_prep_nodewise <- function(mat, na_method = c("pairwise", "listwise"),
                              categorical = NULL) {
  na_method <- .psn_match_arg(na_method)
  if (!anyNA(mat)) return(mat)
  if (na_method == "listwise") {
    return(mat[stats::complete.cases(mat), , drop = FALSE])
  }
  if (is.null(categorical)) categorical <- rep(FALSE, ncol(mat))
  stopifnot("internal: `categorical` must be one logical per column" =
              is.logical(categorical) && length(categorical) == ncol(mat))
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]; miss <- is.na(x)
    if (!any(miss)) next
    obs <- x[!miss]; u <- unique(obs)
    fill <- if (all(u %in% c(0, 1))) {
      as.numeric(mean(obs) >= 0.5)
    } else if (categorical[j]) {
      tab <- table(obs)
      as.numeric(names(tab)[which.max(tab)])          # modal level code
    } else {
      mean(obs)
    }
    mat[miss, j] <- fill
  }
  mat
}

# Invert a correlation matrix. A singular matrix (e.g. a complete-data
# correlation with n <= p, or an exactly collinear column) has no inverse, so
# full-order partial correlations are not identified: the matrix is projected to
# the nearest PD matrix so a network is still returned, and a classed
# `psychnet_singular` warning says that its edges are not interpretable. A
# solvable but numerically singular matrix (reciprocal condition number below
# 1e-12) is inverted as before and warned about the same way.
#' @noRd
.pd_solve <- function(S) {
  out <- tryCatch(solve(S), error = function(e) {
    if (!grepl("singular", conditionMessage(e), fixed = TRUE)) stop(e)
    NULL
  })
  rc <- rcond(S)
  if (is.null(out) || rc < 1e-12) {
    .psn_warn(sprintf(paste0(
      "The correlation matrix is singular (reciprocal condition number %.1e): ",
      "partial correlations are not identified (fewer observations than ",
      "variables, or collinear / duplicate columns). %s; its edges are not ",
      "interpretable. Use a regularized estimator such as ebic_glasso()."),
      rc, if (is.null(out)) "The nearest positive-definite matrix was inverted"
          else "The returned network comes from a near-singular inverse"),
      class = "psychnet_singular")
  }
  if (is.null(out)) out <- solve(.nearest_pd_cor(S))
  out
}

# Two-sided p-values for a (partial) correlation matrix. `k` is the number of
# variables partialled out (0 for marginal correlations, p - 2 for full-order
# partial correlations); df = n - 2 - k. A pair whose df is not positive has no
# test at all, and treating it as p = 1 would silently delete every such edge,
# so it is refused (class `psychnet_degenerate_sample`).
#' @noRd
.cor_pvalues <- function(r, n, k) {
  if (length(n) == 1L) n <- matrix(n, nrow(r), ncol(r))
  if (!identical(dim(n), dim(r)))
    .psn_stop("Effective sample sizes must match the correlation matrix.", call. = FALSE)
  df <- n - 2L - k
  off <- row(df) != col(df)
  if (any(df[off] <= 0)) {
    .psn_stop(sprintf(paste0(
      "The sample size is too small for the significance test: it needs ",
      "n > %d (df = n - 2 - %d > 0) but the smallest %s is %s. Omit `alpha` ",
      "or supply more observations."),
      2L + k, k, if (length(unique(n[off])) > 1L)
        "pairwise sample size" else "n", format(min(n[off]))),
      class = c("psychnet_degenerate_sample", "psychnet_bad_input"))
  }
  P <- matrix(1, nrow(r), ncol(r), dimnames = dimnames(r))
  valid <- df > 0L & is.finite(r)
  tstat <- r[valid] * sqrt(df[valid] / pmax(1 - r[valid]^2, 1e-12))
  P[valid] <- 2 * stats::pt(-abs(tstat), df[valid])
  diag(P) <- 1
  P
}

# Absolute-weight cut-off shared by every thresholding verb: one finite number
# >= 0. A logical is qgraph's EBICglasso significance rule, which psychnets
# does not reimplement; it is refused with a pointer instead of being coerced
# (TRUE would otherwise become a cut-off of 1 and silently empty the network).
#' @noRd
.check_threshold <- function(threshold) {
  if (is.logical(threshold)) {
    .psn_stop(paste0(
      "`threshold` must be a single number >= 0 (an absolute edge-weight ",
      "cut-off). qgraph's logical `threshold = TRUE` rule is not offered; ",
      "give a numeric cut-off such as `threshold = 0.05`, or use `alpha` for ",
      "a significance rule in cor_network() / pcor_network()."))
  }
  .psn_check_number(threshold, "threshold", lower = 0)
}

# Significance level: NULL (no test) or one number strictly between 0 and 1.
#' @noRd
.check_alpha <- function(alpha) {
  .psn_check_number(alpha, "alpha", lower = 0, upper = 1, open = TRUE,
                    null_ok = TRUE)
}

# A user-supplied sample size that is only recorded (or enters a log-likelihood
# scale): NULL, or one finite number >= 1. Whole numbers are not required --
# an effective n from pairwise data may be fractional.
#' @noRd
.check_n_value <- function(n) {
  .psn_check_number(n, "n", lower = 1, null_ok = TRUE)
}

# Kendall's null distribution depends on ties in the original observations;
# a Pearson t test applied to tau tests a different statistic.
#' @noRd
.kendall_pvalues <- function(data, na_method) {
  mat <- .as_numeric_matrix(data, drop_na = identical(na_method, "listwise"))
  p <- ncol(mat)
  out <- matrix(1, p, p, dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    ok <- stats::complete.cases(mat[, c(i, j), drop = FALSE])
    a <- mat[ok, i]; b <- mat[ok, j]
    if (length(a) >= 2L && length(unique(a)) >= 2L &&
        length(unique(b)) >= 2L) {
      out[i, j] <- out[j, i] <- stats::cor.test(a, b, method = "kendall",
                                               exact = FALSE)$p.value
    }
  }
  out
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

# Correlation input shared by the correlation-based verbs: from raw `data`
# (via .cor_input()) or from a validated `cor_matrix`. `n_role` says what the
# verb does with `n` when a matrix is given: "required" (it enters the
# likelihood), "optional" (only a significance test needs it), or "ignored"
# (a structural filter). Labels are checked against the number of nodes that
# survive column handling, so a dropped column is reported by name instead of
# as a raw dimnames error. Returns list(S, n, n_pair, na_used, labels).
#' @noRd
.network_input <- function(data, cor_matrix, n, cor_method, na_method, labels,
                           n_role = c("required", "optional", "ignored")) {
  n_role <- .psn_match_arg(n_role)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    .psn_check_labels(labels, ncol(ci$S))
    return(list(S = ci$S, n = ci$n, n_pair = ci$n_pair, na_used = ci$na_method,
                labels = if (is.null(labels)) ci$labels else labels))
  }
  S <- .check_cor_matrix(cor_matrix)
  if (identical(n_role, "required") && is.null(n)) {
    .psn_stop("`n` is required when `cor_matrix` is supplied.", call. = FALSE)
  }
  if (!identical(n_role, "ignored")) .check_n_value(n)
  .psn_check_labels(labels, ncol(S))
  if (is.null(labels)) {
    labels <- colnames(S)
    if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
  }
  list(S = S, n = n, n_pair = NULL, na_used = "none", labels = labels)
}

#' Correlation network
#'
#' Marginal (zero-order) association network: the correlation matrix of the
#' chosen `cor_method` (Pearson by default) with the diagonal removed.
#' Equivalent to `bootnet`'s `"cor"` default.
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied. Ordered-factor columns enter as their integer
#'   level codes; other non-numeric columns are dropped with a
#'   `psychnet_dropped_columns` message, and constant / all-missing columns
#'   with a `psychnet_constant_variable` warning.
#' @param cor_matrix Optional precomputed correlation matrix; if given, `data`
#'   is ignored and `n` is required when `alpha` is used.
#' @param n Sample size (needed for significance testing when `cor_matrix` is
#'   supplied): a single number >= 1, and a whole number of at least 3 when
#'   `alpha` is set.
#' @param cor_method Correlation method: `"pearson"` (default), `"spearman"`,
#'   `"kendall"`, or `"auto"` (polychoric/polyserial for ordinal items, the
#'   `qgraph::cor_auto` default; see [cor_auto()]).
#' @param threshold Correlations with absolute value below this are set to zero.
#'   A single number >= 0; default 0. (qgraph's logical `threshold = TRUE` is
#'   refused -- use `alpha` for a significance rule.)
#' @param alpha Significance level in (0, 1); if set, correlations not
#'   significant at `alpha` are zeroed. `NULL` (default) keeps every edge.
#'   Kendall marginal correlations use a tie-corrected asymptotic test and
#'   require raw data. Significance testing of full-order partial Kendall
#'   correlations is not supported. A sample too small for the test
#'   (df = n - 2 - k <= 0, with k the number of variables partialled out) is an
#'   error of class `psychnet_degenerate_sample`, never an empty network.
#' @param adjust Multiple-comparison adjustment for the edge p-values (any
#'   [stats::p.adjust] method). Default `"none"`.
#' @param na_method Missing-data handling: `"pairwise"` (default) uses
#'   pairwise-complete correlations projected to the nearest positive-definite
#'   matrix; `"listwise"` drops rows with any `NA`. Identical when data is
#'   complete.
#' @param labels Optional node labels, one unique label per network node.
#' @return A `psychnet` object whose `$weights` is the thresholded correlation
#'   matrix, with `$cor_matrix`, `$n_eff`, `$na_method` (and `$p_values` when
#'   `alpha` is used). `$n_pair` (the pairwise co-observation counts) is present
#'   only when the network is estimated from `data`. Node order in `$weights`
#'   and `$nodes` is always the column order of the input `data` /
#'   `cor_matrix`, never sorted. Invalid arguments raise errors of class
#'   `psychnet_bad_input`.
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' cor_network(x)
#' cor_network(x, alpha = 0.05, adjust = "BH")
#' @export
cor_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                        cor_method = c("pearson", "spearman", "kendall", "auto"),
                        threshold = 0, alpha = NULL, adjust = "none",
                        na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- .psn_match_arg(cor_method)
  adjust <- .psn_match_arg(adjust, stats::p.adjust.methods)
  na_method <- .psn_match_arg(na_method)
  .check_threshold(threshold)
  .check_alpha(alpha)
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "optional")
  S <- inp$S; n <- inp$n; n_pair <- inp$n_pair; labels <- inp$labels
  if (!is.null(alpha) && !is.null(cor_matrix)) {
    if (is.null(n)) {
      .psn_stop("`n` is required when `cor_matrix` is supplied and `alpha` is set.",
           call. = FALSE)
    }
    .psn_check_n(n, min = 3)
  }
  g <- S
  diag(g) <- 0
  r_full <- g                                  # p-values use the true correlations
  g[abs(g) < threshold] <- 0
  extra <- c(list(cor_matrix = S, n_eff = n, na_method = inp$na_used),
             if (!is.null(n_pair)) list(n_pair = n_pair))
  if (!is.null(alpha)) {
    if (identical(cor_method, "kendall") && !is.null(cor_matrix))
      .psn_stop("Kendall significance testing requires raw `data` to account for ties.",
           call. = FALSE)
    P <- if (identical(cor_method, "kendall")) .kendall_pvalues(data, na_method)
    else .cor_pvalues(r_full, if (!is.null(n_pair)) n_pair else n, k = 0L)
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
#' Full-order partial correlations need an invertible correlation matrix. With
#' fewer observations than variables, or collinear / duplicate columns, they
#' are not identified: the nearest positive-definite matrix is inverted so a
#' network is still returned, and a `psychnet_singular` warning says its edges
#' are not interpretable (use a regularized estimator such as [ebic_glasso()]).
#'
#' @inheritParams cor_network
#' @return A `psychnet` object whose `$weights` is the thresholded
#'   partial-correlation matrix, with `$precision`, `$cor_matrix` (and
#'   `$p_values` when `alpha` is used). Node order in `$weights` and `$nodes` is
#'   always the column order of the input `data` / `cor_matrix`, never sorted.
#' @examples
#' x <- matrix(stats::rnorm(200 * 4), 200, 4)
#' pcor_network(x)
#' pcor_network(x, alpha = 0.05, adjust = "holm")
#' @export
pcor_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         threshold = 0, alpha = NULL, adjust = "none",
                         na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- .psn_match_arg(cor_method)
  .check_threshold(threshold)
  .check_alpha(alpha)
  if (!is.null(alpha) && identical(cor_method, "kendall"))
    .psn_stop("Significance testing of full-order partial Kendall correlations is not supported; omit `alpha` or use Pearson correlations.",
         call. = FALSE)
  adjust <- .psn_match_arg(adjust, stats::p.adjust.methods)
  na_method <- .psn_match_arg(na_method)
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "optional")
  S <- inp$S; n <- inp$n; n_pair <- inp$n_pair; labels <- inp$labels
  na_used <- inp$na_used
  if (!is.null(alpha) && !is.null(cor_matrix)) {
    if (is.null(n)) {
      .psn_stop("`n` is required when `cor_matrix` is supplied and `alpha` is set.",
           call. = FALSE)
    }
    .psn_check_n(n, min = 3)
  }
  wi <- .pd_solve(S)
  g  <- .precision_to_pcor(wi)
  r_full <- g                                   # p-values use the true partials
  g[abs(g) < threshold] <- 0
  extra <- c(list(precision = wi, cor_matrix = S, n_eff = n,
                  na_method = na_used),
             if (!is.null(n_pair)) list(n_pair = n_pair))
  if (!is.null(alpha)) {
    if (identical(na_used, "pairwise") &&
        length(unique(n_pair[upper.tri(n_pair)])) > 1L) {
      .psn_stop("`alpha` inference for full-order partial correlations is undefined when pairwise-complete sample sizes differ; use `na_method = 'listwise'` or omit `alpha`.",
           call. = FALSE)
    }
    P <- .cor_pvalues(r_full, n, k = ncol(S) - 2L)   # full-order partials
    g <- .apply_sig(g, P, alpha, adjust)
    extra$p_values <- P
  }
  dimnames(g) <- list(labels, labels)
  .new_psychnet(g, labels, method = "pcor", directed = FALSE,
                n_obs = n, extra = extra)
}
