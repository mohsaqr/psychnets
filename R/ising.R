# Ising model for binary data (van Borkulo et al. 2014), clean-room base R.
# Each node is regressed on all others by L1-penalized logistic regression
# (R/lasso_glm.R) with per-node EBIC selection; the asymmetric nodewise
# estimates are combined by the AND or OR rule and symmetrized.

# Coerce + validate a binary (0/1) matrix, then apply NA handling.
#
# A binary node may arrive as a 0/1 number, a logical, or a two-level factor;
# all three are coded 0/1 here (a factor's FIRST level is 0), so no column is
# ever dropped for its storage type. Anything else -- a factor with more than
# two levels, a character column, a numeric column not coded 0/1 -- is refused
# by name. Returns list(mat, source_names): `source_names` are the column names
# BEFORE constant columns were dropped, so the caller can remap a `labels`
# vector given for every input column.
#' @noRd
.as_binary_matrix <- function(data, na_method = c("pairwise", "listwise"),
                              caller = "ising_fit") {
  na_method <- .psn_match_arg(na_method)
  if (is.null(data)) .psn_stop("`data` is required.")
  if (is.data.frame(data)) {
    if (is.null(names(data)) || ncol(data) == 0L) {
      .psn_stop(sprintf("%s() needs a data frame with at least two columns.", caller))
    }
    nms <- names(data)
    # The training coding of each non-numeric column, kept so net_predict()
    # encodes new data the same way (by label, not by the levels present then).
    codings <- lapply(data, function(x) {
      if (is.logical(x)) c("FALSE", "TRUE")
      else if (is.factor(x)) levels(droplevels(x))
      else NULL
    })
    to01 <- function(x, nm) {
      if (is.logical(x)) return(as.numeric(x))
      if (is.factor(x)) {
        x <- droplevels(x)
        if (nlevels(x) > 2L) {
          .psn_stop(sprintf(paste0(
            "%s() requires binary data; factor column '%s' has %d levels (%s)."),
            caller, nm, nlevels(x), paste(levels(x), collapse = ", ")))
        }
        return(as.numeric(as.integer(x) - 1L))   # first level = 0
      }
      if (is.numeric(x)) return(as.numeric(x))
      .psn_stop(sprintf(paste0(
        "%s() requires binary data; column '%s' is %s. Code it 0/1, as a ",
        "logical, or as a two-level factor."), caller, nm, class(x)[1L]))
    }
    data <- as.data.frame(Map(to01, data, nms), col.names = nms,
                          check.names = FALSE, stringsAsFactors = FALSE)
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    codings <- NULL
    if (!(is.numeric(data) || is.logical(data))) {
      .psn_stop(sprintf("%s() requires a numeric or logical 0/1 matrix.", caller))
    }
  } else {
    .psn_stop("`data` must be a data frame or a numeric/logical matrix.")
  }
  storage.mode(data) <- "double"
  if (is.null(colnames(data))) colnames(data) <- paste0("V", seq_len(ncol(data)))
  source_names <- colnames(data)
  not01 <- vapply(seq_len(ncol(data)), function(j) {
    u <- data[, j]; u <- u[!is.na(u)]
    !all(u == 0 | u == 1)
  }, logical(1))
  if (any(not01)) {
    vals <- vapply(which(not01), function(j) {
      u <- sort(unique(data[!is.na(data[, j]), j]))
      paste(utils::head(format(u, trim = TRUE), 4L), collapse = ", ")
    }, character(1))
    .psn_stop(sprintf(paste0(
      "%s() requires binary (0/1) data; column(s) %s are not coded 0/1 ",
      "(values: %s). Recode them (e.g. x - 1 for 1/2 coding)."),
      caller, paste(source_names[not01], collapse = ", "),
      paste(vals, collapse = "; ")))
  }
  mat <- .as_numeric_matrix(data, drop_na = FALSE)
  kept <- match(colnames(mat), source_names)
  category_levels <- if (is.null(codings)) NULL else unname(codings[kept])
  list(mat = .na_prep_nodewise(mat, na_method), source_names = source_names,
       category_levels = category_levels)
}

# Remap a `labels` vector given for every INPUT column onto the columns that
# survived (constant columns are dropped), then validate it against p.
#' @noRd
.remap_labels <- function(labels, source_names, kept_names) {
  if (is.null(labels)) return(NULL)
  if (length(labels) == length(source_names) &&
      length(kept_names) < length(source_names)) {
    labels <- labels[match(kept_names, source_names)]
  }
  .psn_check_labels(labels, length(kept_names))
  as.character(labels)
}

# Validate the arguments shared by the nodewise lasso estimators.
#' @noRd
.check_nodewise_args <- function(gamma, nlambda, lambda_min_ratio) {
  .psn_check_number(gamma, "gamma", lower = 0)
  .psn_check_number(nlambda, "nlambda", lower = 2, whole = TRUE)
  .psn_check_number(lambda_min_ratio, "lambda_min_ratio", lower = 0, upper = 1,
                    open = TRUE)
  invisible(TRUE)
}

# Drop rows whose total sum-score is below min_sum, keeping weights aligned.
#' @noRd
.apply_min_sum <- function(mat, weights, min_sum) {
  if (is.null(min_sum)) return(list(mat = mat, weights = weights))
  keep <- rowSums(mat) >= min_sum
  if (!any(keep)) .psn_stop("`min_sum` removed every row.", call. = FALSE)
  # Weights are rescaled to max 1 (mgm's convention) only after filtering, so a
  # removed row's weight cannot change the retained rows' scale.
  w <- if (is.null(weights)) NULL else weights[keep]
  if (!is.null(w) && sum(w) <= 0)
    .psn_stop("`min_sum` kept only rows with zero weight.")
  list(mat = mat[keep, , drop = FALSE],
       weights = if (is.null(w)) NULL else w / max(w))
}

#' Ising network for binary data
#'
#' Estimates an Ising model by nodewise L1-penalized logistic regression with
#' EBIC selection, combined by the AND (default) or OR rule. Equivalent in
#' purpose to `IsingFit::IsingFit()`, but pure base R and self-certified: each
#' node's regression reports its stationarity (KKT) residual (see
#' [glm_lasso_kkt()]).
#'
#' @param data Binary data frame or matrix (rows = observations). Columns may
#'   be numeric 0/1, logical, or two-level factors (the first level is coded 0).
#'   Any other column -- a numeric column not coded 0/1 (e.g. 1/2), a factor
#'   with more than two levels, a character column -- is an error naming it,
#'   never a silently dropped node. A constant column carries no information
#'   and is dropped.
#' @param gamma EBIC hyperparameter. Default 0.25.
#' @param rule Edge-combination rule: `"AND"` (default) or `"OR"`.
#' @param nlambda Number of penalties per nodewise path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#'   With `native = FALSE`, omitting the argument retains glmnet's reference
#'   default for compatibility; an explicitly supplied value is honored.
#' @param min_sum Minimum row sum-score (number of endorsed items); rows below
#'   it are dropped before fitting. `NULL` (default) keeps every row.
#' @param weights Optional non-negative observation weights, one per input row
#'   after missing-data preparation. When `min_sum` is used, weights are filtered
#'   alongside their corresponding rows. Weights are rescaled to a maximum of 1
#'   (the `mgm::mgm()` convention), so only their relative sizes matter:
#'   `weights = w` and `weights = 7 * w` give the same network. `NULL`
#'   (default) is unweighted.
#' @param na_method Missing-data handling: `"pairwise"` (default) single-imputes
#'   each column over its observed values (mode for binary), keeping the full
#'   sample; `"listwise"` drops incomplete rows. Identical for complete data.
#' @param native Solver switch. `TRUE` (default) uses psychnet's own pure-R,
#'   dependency-free, self-certified L1 logistic path (KKT ~1e-9). `FALSE`
#'   delegates each per-node fit to the `glmnet` package with the IsingFit EBIC
#'   path, so the returned `$weights`/`$thresholds` byte-match
#'   `IsingFit::IsingFit()` (to ~1e-16) at the cost of glmnet's looser
#'   self-certificate. `native = FALSE` needs the optional `glmnet` package
#'   (Suggests). `min_sum` is supported by both engines; observation `weights`
#'   require `native = TRUE`.
#' @param labels Optional node labels: one per node, or one per input column
#'   (labels of dropped constant columns are then discarded).
#' @return A `psychnet` object whose `$weights` is the symmetric weight matrix,
#'   with `$thresholds` (node intercepts), `$rule`, `$nodewise` (the per-node
#'   fits used by [net_predict()]) and `$kkt` (the worst nodewise stationarity
#'   residual, read with [certificate()]).
#'
#'   Errors of class `psychnet_bad_input` for non-binary columns, invalid
#'   `gamma`/`nlambda`/`lambda_min_ratio`/`weights`, or mislengthed `labels`.
#'   A nodewise fit that does not reach its optimum raises a
#'   `psychnet_no_converge` warning naming the node.
#' @references van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F.,
#'   Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
#'   constructing networks from binary data. *Scientific Reports*, 4, 5918.
#'   \doi{10.1038/srep05918}
#' @examples
#' # Ising models need 0/1 data: split the five SRL construct scores at the median.
#' b <- dichotomize(SRL_GPT)
#' fit <- ising_fit(b)
#' fit
#' as.data.frame(fit)
#' certificate(fit)
#' @export
ising_fit <- function(data, gamma = 0.25, rule = c("AND", "OR"),
                      nlambda = 100L, lambda_min_ratio = 0.01, min_sum = NULL,
                      weights = NULL,
                      na_method = c("pairwise", "listwise"),
                      native = TRUE, labels = NULL) {
  ratio_missing <- missing(lambda_min_ratio)
  rule <- .psn_match_arg(rule)
  na_method <- .psn_match_arg(na_method)
  engine <- .resolve_native(native, "glmnet")
  .check_nodewise_args(gamma, nlambda, lambda_min_ratio)
  .psn_check_number(min_sum, "min_sum", lower = 0, null_ok = TRUE)
  bin <- .as_binary_matrix(data, na_method, caller = "ising_fit")
  mat <- bin$mat
  weights <- .check_weights(weights, nrow(mat))
  ms <- .apply_min_sum(mat, weights, min_sum)
  mat <- ms$mat; weights <- ms$weights
  p <- ncol(mat)
  labels <- .remap_labels(labels, bin$source_names, colnames(mat)) %||%
    colnames(mat)

  if (engine == "glmnet") {
    if (!is.null(weights)) {
      .psn_stop("native = FALSE does not support `weights`; use native = TRUE.",
           call. = FALSE)
    }
    out <- .ising_fit_glmnet(mat, gamma, rule, nlambda,
                             if (ratio_missing) NULL else lambda_min_ratio,
                             labels)
    out$nodewise$category_levels <- bin$category_levels
    return(out)
  }

  fits <- lapply(seq_len(p), function(i) {
    .nodewise_ebic(mat[, -i, drop = FALSE], mat[, i], "binomial",
                   gamma, nlambda, lambda_min_ratio, weights = weights)
  })
  .warn_nodewise_convergence(fits, labels, "ising_fit")

  # Asymmetric matrix B: B[i, j] = effect of node j on node i on the raw 0/1
  # scale (the Ising interaction), with B_std its standardized counterpart.
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  B_std <- matrix(0, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) { B[i, -i] <- fits[[i]]$beta; B_std[i, -i] <- fits[[i]]$beta_std }
  # Use the same weighted centers/scales the nodewise fits were standardized on
  # (lasso_glm.R `.standardize(X, weights)`), so the stored center/scale and the
  # de-centered thresholds stay coherent under non-uniform `weights`.
  std <- .standardize(mat, weights)
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
                n_obs = nrow(mat), data = mat,
                extra = list(thresholds = stats::setNames(thresholds, labels),
                             rule = rule, kkt = worst_kkt, native = native,
                             nlambda = nlambda,
                             lambda_min_ratio = lambda_min_ratio,
                             nodewise = list(intercept = b0_std,
                                             beta_std = B_std,
                                             families = rep("binomial", p),
                                             center = std$center,
                                             scale = std$scale,
                                             category_levels =
                                               bin$category_levels)))
}
