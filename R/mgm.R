# Mixed graphical model (Haslbeck & Waldorp 2020), clean-room base R. Each node
# is regressed on all others with the L1-penalized GLM matching its type --
# gaussian for continuous nodes, logistic for binary nodes -- and the asymmetric
# mgm-scale estimates are combined by the AND rule and symmetrized. v0.1
# supports gaussian and binary nodes (the dominant mixed case); categorical
# nodes with more than two levels are not yet implemented and error explicitly.

# Detect node types, following mgm::mgm's own rule: a factor/character column is
# categorical, and so is a numeric column with <= 10 distinct integer values;
# everything else is gaussian. Returns a list of parallel `type` ("c"/"g") and
# `level` (number of levels for "c", else 1) vectors.
#
# The <= 10-distinct-integers rule is what mgm actually uses, and Likert / count
# items hit it constantly, so a numeric column caught by it is warned about by
# name: it silently switches that node from gaussian to multinomial, which is a
# different model rather than a different encoding.
#' @noRd
.detect_types <- function(data, warn = TRUE) {
  cols <- if (is.data.frame(data)) as.list(data) else
    lapply(seq_len(ncol(data)), function(j) data[, j])
  nms <- if (is.data.frame(data)) names(data) else colnames(data)

  is_cat <- vapply(cols, function(x) {
    if (is.factor(x) || is.character(x)) return(TRUE)
    u <- unique(x[!is.na(x)])
    length(u) <= 10L && all(u == round(u))
  }, logical(1))

  # Warn only for numeric columns promoted by the integer rule, and only when
  # they are not the plain 0/1 binaries that were always categorical.
  promoted <- vapply(seq_along(cols), function(j) {
    x <- cols[[j]]
    if (!is_cat[j] || is.factor(x) || is.character(x)) return(FALSE)
    u <- unique(x[!is.na(x)])
    !(length(u) == 2L && all(u %in% c(0, 1)))
  }, logical(1))
  if (warn && any(promoted)) {
    warning(sprintf(
      "Column(s) %s are numeric with <= 10 distinct integer values and are treated as categorical (mgm's rule), not gaussian. Pass types = to override.",
      paste(nms[promoted], collapse = ", ")), call. = FALSE)
  }

  list(
    type = ifelse(is_cat, "c", "g"),
    level = vapply(seq_along(cols), function(j) {
      if (!is_cat[j]) return(1L)
      length(unique(cols[[j]][!is.na(cols[[j]])]))
    }, integer(1))
  )
}

# Canonical data.frame view for the nodewise fits: every categorical node is a
# factor (so model.matrix() expands it into level dummies), every gaussian node
# is numeric. This is the single place a node's storage type is decided, so the
# glmnet path and the numeric view below can never disagree about it.
#' @noRd
.mgm_as_frame <- function(data, types, levels_v) {
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  cat_ix <- which(types == "c")
  for (j in cat_ix) df[[j]] <- as.factor(df[[j]])
  df
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

#' Mixed graphical model
#'
#' Estimates a mixed graphical model by nodewise L1-penalized regression -- a
#' gaussian (linear) lasso for continuous nodes and a logistic lasso for binary
#' nodes -- with per-node EBIC selection, combined by the AND rule. Equivalent
#' in purpose to `mgm::mgm()`, but pure base R and self-certified: each node's
#' regression reports its stationarity (KKT) residual (see [glm_lasso_kkt()]).
#'
#' @param data Data frame or matrix (rows = observations). Columns are
#'   continuous, binary, or -- with `native = FALSE` -- multi-level categorical
#'   (a `factor`/`character` column, or a numeric one caught by the detection
#'   rule below).
#' @param gamma EBIC hyperparameter. Default 0.25.
#' @param types Optional character vector of node types (`"g"` gaussian, `"c"`
#'   categorical); auto-detected if `NULL`. Detection follows `mgm::mgm()`'s own
#'   rule: a `factor`/`character` column is categorical, and so is a numeric
#'   column with 10 or fewer distinct integer values. A numeric column promoted
#'   by that integer rule is warned about by name, because it switches the node
#'   from gaussian to categorical -- a different model, not a different encoding.
#'   Likert and count items hit this rule routinely.
#' @param nlambda Number of penalties per nodewise path. Default 100.
#' @param lambda_min_ratio Smallest penalty as a fraction of the largest.
#' @param threshold Post-selection coefficient threshold: `"LW"` (default),
#'   `"HW"`, or `"none"`, matching `mgm::mgm()`.
#' @param rule Edge-combination rule: `"AND"` (default) or `"OR"`.
#' @param moderators Optional single column index of a moderator variable. When
#'   supplied, fits a *moderated* MGM (that variable moderates every pairwise
#'   edge) and returns a `psychnet_moderated` object to be read with
#'   [condition()]. Honours `native` like the unmoderated path: the default base
#'   kernel is pure R and KKT-certified, and covers gaussian and binary nodes;
#'   `native = FALSE` uses `glmnet` and additionally allows multi-level
#'   categorical nodes. `weights` are not supported in this mode.
#' @param weights Optional non-negative observation weights, one per row of the
#'   (NA-prepared) data. `NULL` (default) is unweighted.
#' @param na_method Missing-data handling: `"pairwise"` (default) single-imputes
#'   each column over its observed values (mean for continuous, mode for binary),
#'   keeping the full sample; `"listwise"` drops incomplete rows.
#' @param native Solver switch. `TRUE` (default) uses psychnet's own pure-R,
#'   dependency-free, self-certified L1 path (KKT ~1e-9). `FALSE` delegates each
#'   per-node fit to the `glmnet` package with mgm's exact EBIC/LW path (gaussian
#'   lasso for continuous nodes, multinomial lasso for categorical ones), so
#'   the returned edge magnitudes byte-match `abs(mgm::mgm()$pairwise$wadj)`
#'   (to ~1e-6) at the cost of glmnet's looser self-certificate. `native = FALSE`
#'   needs the optional `glmnet` package (Suggests); `weights` are supported with
#'   `native = TRUE` only.
#'
#'   Only `native = FALSE` supports **multi-level** categorical nodes: the base
#'   kernel is a binomial solver, so a node with more than two levels errors
#'   there by name. One-hot encoding is not an equivalent workaround -- it turns
#'   one k-level node into k separate nodes, which is a different model. A
#'   network containing a multi-level node cannot be passed to [net_predict()],
#'   because a k-level predictor occupies k-1 design columns and so has no
#'   one-coefficient-per-variable representation.
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the symmetric standardized
#'   weight matrix, with `$types`, `$levels` (number of levels per node, 1 for
#'   gaussian) and `$kkt` (the worst nodewise residual). Node order in
#'   `$weights` and `$nodes` is always the column order of the input `data`,
#'   never sorted, and a k-level categorical stays ONE node. A
#'   binary-binary edge carries the sign of its nodewise-logistic coefficient;
#'   `mgm::mgm()` reports the same edge as a magnitude only (its sign is
#'   undefined for a categorical-categorical interaction), so compare such edges
#'   on `abs()`. Continuous columns are standardized internally, binary
#'   predictors enter the graph on their 0/1 dummy scale, and binary-response
#'   logit coefficients are converted to `mgm`'s two-class multinomial scale
#'   before edge aggregation. With these conventions the edge magnitudes match
#'   `mgm::mgm` closely for gaussian-gaussian, gaussian-binary, and binary-binary
#'   edges alike; weak edges near the EBIC/threshold boundary can still differ in
#'   support because the penalty is selected on an independent base-R path.
#' @examples
#' set.seed(1)
#' f <- stats::rnorm(400)
#' g1 <- f + stats::rnorm(400); g2 <- f + stats::rnorm(400)
#' b1 <- (f + stats::rnorm(400) > 0) * 1L
#' d <- data.frame(g1 = g1, g2 = g2, b1 = b1, n = stats::rnorm(400))
#' mgm_fit(d)
#' @export
mgm_fit <- function(data, gamma = 0.25, types = NULL,
                    nlambda = 100L, lambda_min_ratio = 0.01,
                    threshold = c("LW", "HW", "none"), rule = c("AND", "OR"),
                    moderators = NULL, weights = NULL,
                    na_method = c("pairwise", "listwise"),
                    native = TRUE, labels = NULL) {
  threshold <- match.arg(threshold)
  rule <- match.arg(rule)
  na_method <- match.arg(na_method)
  engine <- .resolve_native(native, "glmnet")
  stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma >= 0,
            nlambda >= 2L, lambda_min_ratio > 0, lambda_min_ratio < 1)
  # Type/level detection runs on the ORIGINAL data, before any numeric coercion,
  # so a factor column is still recognisable as one. `.as_numeric_matrix()` would
  # drop it and silently remove a node.
  # The promotion warning is only meaningful when detection is actually making
  # the call; an explicit `types` means nothing is being decided silently.
  det <- .detect_types(data, warn = is.null(types))
  if (is.null(types)) {
    types <- det$type
    levels_v <- det$level
  } else {
    levels_v <- ifelse(types == "c", det$level, 1L)
  }

  # Only the glmnet engine has the multinomial machinery for a >2-level node;
  # the base kernel is a binomial solver. Reject early and by name, rather than
  # letting a k-level column reach a logistic fit that diverges silently.
  multi <- types == "c" & levels_v > 2L
  nms <- if (is.data.frame(data)) names(data) else colnames(data)
  if (any(multi) && engine != "glmnet") {
    stop(sprintf(
      "Column(s) %s are categorical with more than 2 levels; the base kernel is a binomial solver and supports gaussian and binary (0/1) nodes only. Pass native = FALSE to use the multinomial glmnet engine.",
      paste(nms[multi], collapse = ", ")), call. = FALSE)
  }
  # A factor column can only travel the glmnet path; on the base path it would
  # be dropped by the numeric coercion below.
  if (is.data.frame(data)) {
    nonnum <- !vapply(data, is.numeric, logical(1))
    if (any(nonnum) && engine != "glmnet") {
      stop(sprintf(
        "mgm_fit() requires numeric columns when native = TRUE; column(s) %s are non-numeric. Recode a binary column to 0/1, or pass native = FALSE to model it as a categorical node.",
        paste(names(data)[nonnum], collapse = ", ")), call. = FALSE)
    }
  }

  # The glmnet path keeps a data.frame (factors intact) so model.matrix() can
  # expand a k-level node into k-1 dummy columns; the base path keeps the numeric
  # matrix it has always used. `mat` stays the numeric view in both cases -- it
  # is what supplies nrow(), the up-front scaling, and the stored $data.
  df <- .mgm_as_frame(data, types, levels_v)
  mat <- .na_prep_nodewise(.as_numeric_matrix(.mgm_frame_numeric(df),
                                              drop_na = FALSE), na_method)
  p <- ncol(mat)
  if (!is.null(labels)) stopifnot(length(labels) == p)
  if (is.null(labels)) labels <- colnames(mat)
  stopifnot(length(types) == p, all(types %in% c("g", "c")))

  # Moderated MGM: a chosen variable moderates every edge. Dispatched before the
  # binary 0/1 enforcement because the moderated path treats every categorical as
  # a factor; the glmnet engine additionally allows multi-level ones (matching
  # mgm). Honours `native` like every other path: the default base kernel is
  # pure R and KKT-certified. Weights unsupported. Returns a psychnet_moderated,
  # read via condition().
  if (!is.null(moderators)) {
    if (!is.null(weights)) {
      stop("moderated MGM does not support `weights`.", call. = FALSE)
    }
    thr_mod <- if (threshold == "LW") "LW" else "none"
    return(.mmg_estimate(mat, types, moderators, gamma = gamma, rule = rule,
                         threshold = thr_mod, labels = labels,
                         engine = engine, nlambda = nlambda,
                         lambda_min_ratio = lambda_min_ratio))
  }

  # A user-declared binary ('c') column must actually be 0/1 on the base path,
  # else the logistic nodewise fit diverges silently (auto-detected 'c' columns
  # pass trivially). The glmnet path recodes via factor levels, so any coding
  # works there.
  cbin <- which(types == "c" & levels_v <= 2L)
  if (length(cbin) && engine != "glmnet") {
    bad <- cbin[!vapply(cbin, function(j) all(mat[, j] %in% c(0, 1)), logical(1))]
    if (length(bad)) {
      stop(sprintf(
        "types declares column(s) %s as binary ('c'), but they are not coded 0/1.",
        paste(labels[bad], collapse = ", ")), call. = FALSE)
    }
  }
  weights <- .check_weights(weights, nrow(mat))

  if (engine == "glmnet") {
    if (!is.null(weights)) {
      stop("native = FALSE does not support `weights`; use native = TRUE.",
           call. = FALSE)
    }
    return(.mgm_fit_glmnet(df, types, gamma, threshold, rule, nlambda, labels,
                           levels_v = levels_v, mat = mat))
  }

  # Scale continuous columns to unit variance up front (mgm::mgm scale = TRUE):
  # a gaussian response on its raw scale inflates its edge weights by its SD, so
  # the gaussian-node magnitudes only match mgm once each continuous variable is
  # standardized. resp_center/resp_scale record this so net_predict() can put
  # the gaussian response on the same scale as the fitted linear predictor.
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
                n_obs = nrow(mat), data = mat,
                extra = list(types = stats::setNames(types, labels),
                             kkt = worst_kkt, native = native,
                             threshold = threshold,
                             nodewise = list(intercept = intercepts,
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
