# Mixed graphical model (Haslbeck & Waldorp 2020), clean-room base R. Each node
# is regressed on all others with the L1-penalized GLM matching its type --
# gaussian for continuous nodes, logistic for binary nodes -- and the asymmetric
# mgm-scale estimates are combined by the AND rule and symmetrized. The base
# kernel supports gaussian and binary nodes (the dominant mixed case); a
# categorical node with more than two levels needs the multinomial glmnet engine
# (`native = FALSE`) and errors explicitly, by name, on the base kernel.

# Detect node types, following mgm::mgm's own rule: a factor/character column is
# categorical, and so is a numeric column with <= 10 distinct integer values;
# everything else is gaussian. Returns a list of parallel `type` ("c"/"g") and
# `level` vectors. `level` is the number of distinct non-NA values of EVERY
# column (not just the categorical ones), so a caller that overrides a
# detected-gaussian column to "c" via `types=` still gets its true level count.
#
# The <= 10-distinct-integers rule is what mgm actually uses, and Likert / count
# items hit it constantly, so a numeric column caught by it is warned about by
# name: it silently switches that node from gaussian to multinomial, which is a
# different model rather than a different encoding. The warning carries the
# condition class "psychnet_type_promotion" so a resampling verb, which has
# already surfaced it once on the observed fit, can muffle the repeats.
#' @noRd
.detect_types <- function(data, warn = TRUE) {
  cols <- if (is.data.frame(data)) as.list(data) else
    lapply(seq_len(ncol(data)), function(j) data[, j])
  nms <- if (is.data.frame(data)) names(data) else colnames(data)

  n_lev <- vapply(cols, function(x) length(unique(x[!is.na(x)])), integer(1))
  is_cat <- vapply(cols, function(x) {
    if (is.factor(x) || is.character(x)) return(TRUE)
    u <- unique(x[!is.na(x)])
    length(u) <= 10L && all(u == round(u))
  }, logical(1))

  # A character column with more than 10 distinct values is almost always an
  # identifier or free text, not a node; refuse it rather than fit a k-hundred
  # class multinomial. An explicit factor() is the caller saying "I mean it".
  id_like <- vapply(cols, is.character, logical(1)) & n_lev > 10L
  if (any(id_like)) {
    .psn_stop(sprintf(
      "Column(s) %s are character with more than 10 distinct values and look like identifiers, not nodes. Drop them, or convert to factor() to model them as categorical.",
      paste(nms[id_like], collapse = ", ")), call. = FALSE)
  }

  # Warn only for numeric columns promoted by the integer rule, and only when
  # they have more than two values: a two-valued column is a binary node
  # however it is coded (0/1, 1/2, ...) -- it is recoded to its 0/1 level index
  # downstream, so treating it as categorical changes nothing.
  promoted <- vapply(seq_along(cols), function(j) {
    x <- cols[[j]]
    if (!is_cat[j] || is.factor(x) || is.character(x)) return(FALSE)
    length(unique(x[!is.na(x)])) > 2L
  }, logical(1))
  if (warn && any(promoted)) {
    warning(structure(
      class = c("psychnet_type_promotion", "warning", "condition"),
      list(message = sprintf(
        "Column(s) %s are numeric with <= 10 distinct integer values and are treated as categorical (mgm's rule), not gaussian. Pass types = to override.",
        paste(nms[promoted], collapse = ", ")), call = NULL)))
  }

  list(type = ifelse(is_cat, "c", "g"), level = n_lev)
}

# Columns with fewer than two distinct non-NA values carry no information and
# are dropped by `.as_numeric_matrix()` downstream. Drop them from the ORIGINAL
# data first, so type detection, `types`, `levels_v` and the numeric matrix all
# describe the same set of columns.
#' @noRd
.usable_columns <- function(data) {
  cols <- if (is.data.frame(data)) as.list(data) else
    lapply(seq_len(ncol(data)), function(j) data[, j])
  vapply(cols, function(x) length(unique(x[!is.na(x)])) >= 2L, logical(1))
}

# Canonical data.frame view for the nodewise fits: every categorical node is a
# factor (so model.matrix() expands it into level dummies), every gaussian node
# is numeric. This is the single place a node's storage type is decided, so the
# glmnet path and the numeric view below can never disagree about it.
#' @noRd
.mgm_as_frame <- function(data, types, levels_v) {
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  cat_ix <- which(types == "c")
  for (j in cat_ix) df[[j]] <- droplevels(as.factor(df[[j]]))
  df
}

# Resolve `moderators` (a column index or a column name) to a single column
# index of `data`, refusing anything else by class.
#' @noRd
.mgm_moderator_index <- function(moderators, data) {
  nms <- colnames(data)
  if (is.character(moderators)) {
    if (length(moderators) != 1L || is.na(moderators) ||
        !(moderators %in% nms)) {
      .psn_stop(sprintf(
        "`moderators` must name one column of `data`; '%s' is not one.",
        paste(moderators, collapse = "', '")))
    }
    return(match(moderators, nms))
  }
  .psn_check_number(moderators, "moderators", lower = 1, upper = ncol(data),
                    whole = TRUE)
  as.integer(round(moderators))
}

# Numeric view of that frame: a factor becomes its 0-based level index, which
# reproduces the 0/1 coding for a binary node and gives a k-level node the
# 0..k-1 integer coding. Supplies nrow(), the up-front gaussian scaling, and the
# stored $data, none of which can hold a factor.
#' @noRd
.mgm_frame_numeric <- function(df) {
  out <- lapply(df, function(x) if (is.factor(x)) as.integer(x) - 1L else x)
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Estimate a mixed graphical model for continuous and binary variables
#'
#' @description
#' `mgm_fit()` estimates a mixed graphical model: one network for variables
#' of different types, such as continuous symptom ratings and binary
#' indicators. An edge connects two variables that remain associated after
#' conditioning on all other variables, and an absent edge indicates
#' conditional independence (Haslbeck & Waldorp, 2020), so the network shows
#' which variables are directly related when they cannot all be treated as
#' Gaussian or all as binary.
#'
#' The function regresses each variable on all others with the L1-penalized
#' (lasso) regression suited to its type, linear for continuous and logistic
#' for binary variables, and selects each variable's penalty by the extended
#' Bayesian information criterion (EBIC). Small coefficients are then removed
#' by a threshold (`threshold`, default `"LW"`), and the two estimates of
#' each pair are combined by the AND rule (default) or the OR rule into an
#' undirected edge whose weight is their mean. Variable types are detected
#' from the data or given in `types`. Each regression is checked against its
#' optimality (Karush–Kuhn–Tucker) conditions; [certificate()] reports the
#' largest violation. A threshold that removes a coefficient changes the
#' network after fitting, so its certificate is then reported as missing;
#' `threshold = "none"` keeps it. With `moderators`, the function instead
#' fits a moderated model in which one variable can change the strength of
#' every edge; see [condition()].
#'
#' The result is an undirected network of class `psychnet`
#' ([psychnet-class]), or a `psychnet_moderated` fit with `moderators`.
#'
#' For binary data only, use [ising_fit()]; for continuous data only,
#' [ebic_glasso()]. Categorical variables with more than two levels need
#' `native = FALSE`. The model is the one fitted by `mgm::mgm()`.
#'
#' @details
#' Continuous variables are standardized before fitting, binary predictors
#' enter as 0/1 indicators, and the coefficient of a binary response is
#' expressed on the two-class multinomial scale of `mgm::mgm()` (half the
#' logistic coefficient). Edge magnitudes are therefore comparable with those
#' of `mgm::mgm()`. A binary-binary edge keeps the sign of its logistic
#' coefficient, whereas `mgm::mgm()` reports such an edge without a sign;
#' compare those edges on their absolute value. Weak edges near the selection
#' or threshold boundary can differ between the two engines.
#'
#' With \eqn{\beta_i} the selected coefficients of variable \eqn{i}, \eqn{q}
#' the number of its predictor columns (\eqn{p - 1} for \eqn{p} continuous or
#' binary variables) and \eqn{n} the sample size, the `"LW"` threshold
#' removes coefficients smaller in absolute value than
#' \eqn{\|\beta_i\|_2 \sqrt{\log(q) / n}}, and the `"HW"` threshold those
#' smaller than \eqn{\sqrt{\log(q) / n}}.
#'
#' Node order is the column order of `data`, and a categorical variable with
#' several levels remains one node.
#'
#' @param data A data frame or matrix with one row per observation. Columns
#'   are continuous, binary, or (with `native = FALSE`) categorical with more
#'   than two levels. With `native = TRUE` every column must be numeric.
#'   Columns with fewer than two distinct observed values are dropped. A
#'   character column with more than 10 distinct values is refused as an
#'   identifier; convert it with `factor()` if it is a categorical variable.
#' @param gamma EBIC hyperparameter, a single number >= 0; larger values give
#'   sparser networks. Default `0.25`.
#' @param types Node types, one per column of `data`: `"g"` (Gaussian) or
#'   `"c"` (categorical). Default `NULL` detects them as `mgm::mgm()` does:
#'   factor and character columns are categorical, and so are numeric columns
#'   with at most 10 distinct integer values. Such a numeric column with more
#'   than two values is named in a warning, because it changes the model for
#'   that node; Likert items often meet this rule. A non-factor column
#'   declared `"c"` with more than 10 distinct values is refused.
#' @param nlambda Number of penalty values on each nodewise path, a whole
#'   number >= 2. Default `100`.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest, in
#'   (0, 1). Default `0.01`. With `native = FALSE`, leaving it unset uses
#'   glmnet's own default, as `mgm::mgm()` does.
#' @param threshold Threshold on the selected coefficients: `"LW"` (default),
#'   `"HW"` or `"none"` (see Details), as in `mgm::mgm()`.
#' @param rule `"AND"` (default) keeps an edge when both nodewise estimates
#'   are nonzero; `"OR"` when either is.
#' @param moderators Optional moderator: one column of `data`, by name or by
#'   position. The moderated model lets this variable change the strength of
#'   every pairwise edge and returns a `psychnet_moderated` fit, read with
#'   [condition()] and `as.data.frame()`. Needs at least three variables and
#'   does not accept `weights`. Default `NULL`.
#' @param weights Optional non-negative observation weights, one per row after
#'   missing-data handling. Only their relative sizes matter, as in
#'   `mgm::mgm()`. Requires `native = TRUE`. Default `NULL` (unweighted).
#' @param na_method `"pairwise"` (default) fills each missing value with the
#'   variable's mean (continuous) or most frequent level (categorical) and
#'   keeps every row; `"listwise"` drops rows with a missing value. The two
#'   agree for complete data.
#' @param native `TRUE` (default) uses the package's own solver, for
#'   continuous and binary variables. `FALSE` fits each regression with the
#'   `glmnet` package (which must be installed), following `mgm::mgm()`; it
#'   also allows categorical variables with more than two levels, but not
#'   `weights`, and its certificate reflects glmnet's convergence tolerance.
#'   A network with a multi-level categorical node cannot be passed to
#'   [net_predict()]. One-hot encoding such a variable is not equivalent: it
#'   turns one node into several.
#' @param labels Optional node labels: one per node, or one per input column
#'   (labels of dropped constant columns are then discarded).
#' @return An undirected network of class `psychnet` ([psychnet-class]); with
#'   `moderators`, a `psychnet_moderated` fit (see [condition()]).
#'   `as.data.frame()` gives one row per edge (`from`, `to`, `weight`),
#'   [certificate()] the largest optimality violation of the nodewise
#'   regressions (missing when the threshold removed a coefficient), and
#'   [net_predict()] the predictability of each node. The object also stores
#'   the node types (`types`), the number of levels per node (`levels`, 1 for
#'   continuous) and the certificate of the regressions before thresholding
#'   (`fit_kkt`).
#' @section Conditions:
#'   `psychnet_bad_input` when `data` is not a data frame or matrix, `types`
#'   is not one `"g"`/`"c"` per column, a character column looks like an
#'   identifier, a non-factor column declared `"c"` has more than 10 values,
#'   a categorical variable has more than two levels or a column is not
#'   numeric with `native = TRUE`, `moderators` does not name one column or
#'   names a dropped constant column, a moderated model has fewer than three
#'   variables or is given `weights`, `weights` is combined with
#'   `native = FALSE`, `native = FALSE` is requested without the `glmnet`
#'   package, or `gamma`, `nlambda`, `lambda_min_ratio`, `weights` or
#'   `labels` is invalid; `psychnet_type_promotion` (warning) when a numeric
#'   column with more than two values is detected as categorical;
#'   `psychnet_no_converge` (warning) when a nodewise regression does not
#'   reach its optimum; `psychnet_degenerate_sample` (warning) when a binary
#'   variable takes a single value on every row with positive weight. The
#'   last two are raised with `native = TRUE`.
#' @seealso [condition()] for the moderated model; [ising_fit()] for binary
#'   data; [ebic_glasso()] for continuous data; [certificate()] and
#'   [glm_lasso_kkt()] for the optimality check; [net_predict()];
#'   [psychnet()] with `method = "mgm"`.
#' @references Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
#'   time-varying mixed graphical models in high-dimensional data. *Journal of
#'   Statistical Software*, 93(8), 1-46. \doi{10.18637/jss.v093.i08}
#'
#'   Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
#'   for Gaussian graphical models. *Advances in Neural Information Processing
#'   Systems*, 23, 604-612.
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(400)
#' g1 <- f + stats::rnorm(400); g2 <- f + stats::rnorm(400)
#' b1 <- (f + stats::rnorm(400) > 0) * 1L
#' d <- data.frame(g1 = g1, g2 = g2, b1 = b1, n = stats::rnorm(400))
#' fit <- mgm_fit(d)
#' fit
#' as.data.frame(fit)
#' @export
mgm_fit <- function(data, gamma = 0.25, types = NULL,
                    nlambda = 100L, lambda_min_ratio = 0.01,
                    threshold = c("LW", "HW", "none"), rule = c("AND", "OR"),
                    moderators = NULL, weights = NULL,
                    na_method = c("pairwise", "listwise"),
                    native = TRUE, labels = NULL) {
  ratio_missing <- missing(lambda_min_ratio)
  threshold <- .psn_match_arg(threshold)
  rule <- .psn_match_arg(rule)
  na_method <- .psn_match_arg(na_method)
  engine <- .resolve_native(native, "glmnet")
  .check_nodewise_args(gamma, nlambda, lambda_min_ratio)
  if (!is.data.frame(data) && !is.matrix(data)) {
    .psn_stop("`data` must be a data frame or matrix.")
  }
  if (!is.null(moderators)) moderators <- .mgm_moderator_index(moderators, data)
  if (!is.null(types) &&
      (!is.character(types) || length(types) != ncol(data) ||
       !all(types %in% c("g", "c")))) {
    .psn_stop(sprintf(
      "`types` must be a character vector of \"g\"/\"c\", one per column of `data` (%d).",
      ncol(data)))
  }
  # Listwise deletion happens FIRST, on the original data, so that column
  # usability, type detection and level counts describe the rows that are
  # actually modelled: a level seen only in an incomplete row is not a level,
  # and a node left with one level after deletion is a constant column, dropped
  # like any other. (Pairwise imputation cannot remove a level -- the modal
  # fill is an observed value -- so it stays where it was.)
  if (na_method == "listwise") {
    data <- data[stats::complete.cases(data), , drop = FALSE]
  }
  # Constant / all-NA columns are dropped up front (they would be dropped by
  # `.as_numeric_matrix()` anyway), so everything derived from `data` below --
  # `types`, `levels_v`, the frame and the numeric matrix -- stays aligned.
  usable <- .usable_columns(data)
  if (!all(usable)) {
    orig_nms <- if (is.data.frame(data)) names(data) else colnames(data)
    # Everything indexed against the ORIGINAL columns must follow the drop:
    # `types`, a `labels` vector given for all original columns, and above all
    # a positional `moderators` index -- left untranslated it would silently
    # point at a different variable, and a moderator that is itself constant
    # cannot moderate anything.
    if (!is.null(types)) types <- types[usable]
    if (!is.null(labels) && length(labels) == length(usable)) labels <- labels[usable]
    if (!is.null(moderators)) {
      if (!usable[moderators]) {
        .psn_stop(sprintf(
          "moderators = %d refers to column '%s', which has fewer than two distinct values and was dropped.",
          as.integer(moderators), orig_nms[moderators]), call. = FALSE)
      }
      moderators <- match(moderators, which(usable))
    }
    data <- if (is.data.frame(data)) data[usable] else data[, usable, drop = FALSE]
  }
  # Type/level detection runs on the ORIGINAL data, before any numeric coercion,
  # so a factor column is still recognisable as one. `.as_numeric_matrix()` would
  # drop it and silently remove a node.
  # The promotion warning is only meaningful when detection is actually making
  # the call; an explicit `types` means nothing is being decided silently.
  det <- .detect_types(data, warn = is.null(types))
  if (is.null(types)) types <- det$type
  levels_v <- ifelse(types == "c", det$level, 1L)

  # A NON-factor column declared categorical via `types` but carrying more than
  # 10 distinct values is a continuous variable mislabelled, not a k-class node
  # (10 is the bound mgm's own detection rule uses); refuse it rather than point
  # the caller at a hundred-class multinomial. An explicit `factor` is exempt:
  # that is the caller saying "I mean it", and it is the documented way to model
  # a legitimate > 10-level categorical.
  nms <- if (is.data.frame(data)) names(data) else colnames(data)
  is_fac <- if (is.data.frame(data)) vapply(data, is.factor, logical(1)) else
    rep(FALSE, ncol(data))
  too_many <- types == "c" & levels_v > 10L & !is_fac
  if (any(too_many)) {
    .psn_stop(sprintf(
      "types declares column(s) %s as categorical ('c'), but they have more than 10 distinct values and are not coded as levels; declare them 'g' (gaussian) or recode them.",
      paste(nms[too_many], collapse = ", ")), call. = FALSE)
  }

  # Only the glmnet engine has the multinomial machinery for a >2-level node;
  # the base kernel is a binomial solver. Reject early and by name, rather than
  # letting a k-level column reach a logistic fit that diverges silently.
  multi <- types == "c" & levels_v > 2L
  if (any(multi) && engine != "glmnet") {
    .psn_stop(sprintf(
      "Column(s) %s are categorical with more than 2 levels; the base kernel is a binomial solver and supports gaussian and binary (0/1) nodes only. Pass native = FALSE to use the multinomial glmnet engine.",
      paste(nms[multi], collapse = ", ")), call. = FALSE)
  }
  # A factor column can only travel the glmnet path; on the base path it would
  # be dropped by the numeric coercion below.
  if (is.data.frame(data)) {
    nonnum <- !vapply(data, is.numeric, logical(1))
    if (any(nonnum) && engine != "glmnet") {
      .psn_stop(sprintf(
        "mgm_fit() requires numeric columns when native = TRUE; column(s) %s are non-numeric. Recode a binary column to 0/1, or pass native = FALSE to model it as a categorical node.",
        paste(names(data)[nonnum], collapse = ", ")), call. = FALSE)
    }
  }

  # `mat` is the numeric view: every categorical node as its 0..k-1 level code,
  # missing values handled per `na_method` (mode-imputed for categorical nodes,
  # mean-imputed for gaussian ones, or rows dropped listwise). It supplies
  # nrow(), the up-front scaling, and the stored $data. The glmnet path needs a
  # data.frame with real factors so model.matrix() can expand a k-level node
  # into k-1 dummies; that frame is rebuilt FROM `mat`, after imputation, so the
  # two views can never disagree about rows or values.
  training_frame <- .mgm_as_frame(data, types, levels_v)
  source_names <- names(training_frame)
  category_levels <- lapply(seq_along(types), function(j)
    if (types[j] == "c") levels(training_frame[[j]]) else NULL)
  mat <- .na_prep_nodewise(
    .as_numeric_matrix(.mgm_frame_numeric(training_frame),
                       drop_na = FALSE),
    na_method, categorical = types == "c")
  df <- .mgm_as_frame(mat, types, levels_v)
  p <- ncol(mat)
  .psn_check_labels(labels, p)
  labels <- if (is.null(labels)) colnames(mat) else as.character(labels)
  stopifnot("internal: one type per node" = length(types) == p)

  # Moderated MGM: a chosen variable moderates every edge. Dispatched before the
  # binary 0/1 enforcement because the moderated path treats every categorical as
  # a factor; the glmnet engine additionally allows multi-level ones (matching
  # mgm). Honours `native` like every other path: the default base kernel is
  # pure R and KKT-certified. Weights unsupported. Returns a psychnet_moderated,
  # read via condition().
  if (!is.null(moderators)) {
    if (!is.null(weights)) {
      .psn_stop("moderated MGM does not support `weights`.", call. = FALSE)
    }
    return(.mmg_estimate(mat, types, moderators, gamma = gamma, rule = rule,
                         threshold = threshold, labels = labels,
                         engine = engine, nlambda = nlambda,
                         lambda_min_ratio = if (engine == "glmnet" && ratio_missing)
                           NULL else lambda_min_ratio,
                         moderator_levels = category_levels[[moderators]]))
  }

  # A user-declared binary ('c') column must actually be 0/1 on the base path,
  # else the logistic nodewise fit diverges silently (auto-detected 'c' columns
  # pass trivially). The glmnet path recodes via factor levels, so any coding
  # works there.
  cbin <- which(types == "c" & levels_v <= 2L)
  if (length(cbin) && engine != "glmnet") {
    bad <- cbin[!vapply(cbin, function(j) all(mat[, j] %in% c(0, 1)), logical(1))]
    if (length(bad)) {
      .psn_stop(sprintf(
        "types declares column(s) %s as binary ('c'), but they are not coded 0/1.",
        paste(labels[bad], collapse = ", ")), call. = FALSE)
    }
  }
  weights <- .check_weights(weights, nrow(mat))

  if (engine == "glmnet") {
    if (!is.null(weights)) {
      .psn_stop("native = FALSE does not support `weights`; use native = TRUE.",
           call. = FALSE)
    }
    return(.mgm_fit_glmnet(df, types, gamma, threshold, rule, nlambda,
                           if (ratio_missing) NULL else lambda_min_ratio, labels,
                           levels_v = levels_v, mat = mat,
                           source_names = source_names,
                           category_levels = category_levels))
  }

  # Scale continuous columns to unit variance up front (mgm::mgm scale = TRUE):
  # a gaussian response on its raw scale inflates its edge weights by its SD, so
  # the gaussian-node magnitudes only match mgm once each continuous variable is
  # standardized. resp_center/resp_scale record this so net_predict() can put
  # the gaussian response on the same scale as the fitted linear predictor.
  raw_mat <- mat
  resp_center <- numeric(p); resp_scale <- rep(1, p)
  gix <- which(types == "g")
  if (length(gix)) {
    G <- mat[, gix, drop = FALSE]
    if (is.null(weights)) {
      resp_center[gix] <- colMeans(G)
      resp_scale[gix]  <- pmax(apply(G, 2L, stats::sd), 1e-12)
    } else {
      # frequency-weighted center/scale; the (sum(w) - 1) divisor makes uniform
      # weights reduce exactly to the unweighted sample-sd scaling above.
      sw <- sum(weights)
      resp_center[gix] <- colSums(weights * G) / sw
      Gc <- sweep(G, 2L, resp_center[gix], "-")
      resp_scale[gix]  <- pmax(sqrt(colSums(weights * Gc^2) / max(sw - 1, 1e-12)), 1e-12)
    }
    mat[, gix] <- sweep(sweep(mat[, gix, drop = FALSE], 2L, resp_center[gix], "-"),
                        2L, resp_scale[gix], "/")
  }

  fits <- lapply(seq_len(p), function(i) {
    fam <- if (types[i] == "c") "binomial" else "gaussian"
    .nodewise_ebic(mat[, -i, drop = FALSE], mat[, i], fam,
                   gamma, nlambda, lambda_min_ratio, weights = weights)
  })
  .warn_nodewise_convergence(fits, labels, "mgm_fit")

  # Asymmetric graph weights follow mgm::mgm's reported coefficient convention:
  # coefficients are on the model-matrix scale (binary predictors are unscaled
  # 0/1 dummies), and a two-class categorical response reports one class
  # coefficient, i.e. half of the binary-logit coefficient. The standardized
  # coefficients are still stored for KKT/predictability reconstruction.
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  B_std <- matrix(0, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) {
    beta_graph <- fits[[i]]$beta
    if (types[i] == "c") beta_graph <- beta_graph / 2
    B[i, -i] <- beta_graph
    B_std[i, -i] <- fits[[i]]$beta_std
  }
  npar <- p - 1L
  tau <- rep(0, p)
  Bt <- B
  if (threshold != "none") {
    for (i in seq_len(p)) {
      beta_i <- B[i, -i]
      tau[i] <- if (threshold == "LW") {
        sqrt(sum(beta_i^2)) * sqrt(log(npar) / nrow(mat))
      } else {
        sqrt(log(npar) / nrow(mat))
      }
      Bt[i, abs(Bt[i, ]) < tau[i]] <- 0
    }
  }
  B_std_t <- B_std
  B_std_t[Bt == 0] <- 0
  intercepts <- vapply(fits, function(f) f$b0, numeric(1))
  worst_kkt <- max(vapply(fits, function(f) f$kkt, numeric(1)))
  # The threshold edits the penalized solution after the fact; when it removed
  # anything the returned graph is not the certified optimum (glasso/huge rule).
  thresholded <- any(B != 0 & Bt == 0)
  # Weighted centers/scales, matching the nodewise fits' standardization, so the
  # composite back-transform below is coherent under non-uniform `weights`.
  std <- .standardize(mat, weights)
  # Composite raw -> standardized-predictor transform for net_predict(), which
  # receives the user's raw data: raw -> (up-front scaling) -> nodewise standard.
  comp_center <- resp_center + resp_scale * std$center
  comp_scale  <- resp_scale * std$scale
  families <- ifelse(types == "c", "binomial", "gaussian")

  present <- if (rule == "AND") (Bt != 0) & (t(Bt) != 0)
             else (Bt != 0) | (t(Bt) != 0)
  W <- (Bt + t(Bt)) / 2
  W[!present] <- 0
  diag(W) <- 0

  .new_psychnet(W, labels, method = "mgm", directed = FALSE,
                n_obs = nrow(mat), data = raw_mat,
                extra = list(types = stats::setNames(types, labels),
                             levels = stats::setNames(as.integer(levels_v), labels),
                             kkt = if (thresholded) NA_real_ else worst_kkt,
                             fit_kkt = worst_kkt,
                             certificate_target = if (thresholded)
                               "pre_threshold_fit" else "returned_network",
                             native = native,
                             threshold = threshold,
                             nlambda = nlambda,
                             lambda_min_ratio = lambda_min_ratio,
                             nodewise = list(intercept = intercepts,
                                             source_names = source_names,
                                             category_levels = category_levels,
                                             beta_std = B_std,
                                             beta_std_thresholded = B_std_t,
                                             beta_graph = B,
                                             beta_graph_thresholded = Bt,
                                             tau = stats::setNames(tau, labels),
                                             families = families,
                                             center = comp_center,
                                             scale = comp_scale,
                                             resp_center = resp_center,
                                             resp_scale = resp_scale)))
}
