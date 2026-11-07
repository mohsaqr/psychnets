# Moderated Mixed Graphical Model (Haslbeck "group differences via moderation").
# Ported verbatim in its numerical logic from Nestimate's verified implementation
# (which matches mgm::mgm(..., moderators = k) + mgm::condition() to machine
# precision); only the argument names are adapted to psychnet conventions
# (gamma, rule, types) and the outputs are wrapped as psychnet objects. A chosen
# variable moderates every pairwise edge; group-structure differences appear as
# nonzero moderator x endpoint interactions. Two engines: the default base
# engine runs every nodewise fit on psychnets' own penalized-IRLS kernel
# (R/lasso_glm.R, KKT-certified, gaussian + binary nodes); `native = FALSE`
# uses glmnet, the reference path, which also covers multi-level categoricals.

# Per-node main-effect magnitude for V_target (mean |coef| over its dummies).
#' @noRd
.mmg_main_mag <- function(node_fit, target) {
  cn <- node_fit$col_names
  idx <- which(!grepl(":", cn, fixed = TRUE) &
                 grepl(paste0("V", target, "."), cn, fixed = TRUE))
  if (!length(idx)) return(0)
  vals <- if (node_fit$multinomial)
    unlist(lapply(node_fit$beta_list, function(b) b[idx])) else node_fit$beta[idx]
  mean(abs(vals))
}

# Per-node interaction magnitude for the (target, other) pair.
#' @noRd
.mmg_int_mag <- function(node_fit, target, other) {
  cn <- node_fit$col_names
  has_int <- grepl(":", cn, fixed = TRUE)
  idx <- which(grepl(paste0("V", target, "."), cn, fixed = TRUE) & has_int &
               grepl(paste0("V", other, "."), cn, fixed = TRUE))
  if (!length(idx)) return(0)
  vals <- if (node_fit$multinomial)
    unlist(lapply(node_fit$beta_list, function(b) b[idx])) else node_fit$beta[idx]
  mean(abs(vals))
}

# Conditioning indicator for one interaction coefficient: for a categorical
# moderator, 1 when the dummy level equals mod_value, else 0; for a continuous
# moderator, mod_value itself.
#' `mod_value` is on the INTERNAL scale: the 0-based level code of a
#' categorical moderator, the standardized (z) value of a continuous one --
#' condition() maps the user's value onto it with `.mmg_moderator_value()`.
#' The model-matrix dummy of a categorical moderator is named "V<k>.<code>".
#' @noRd
.mmg_cond_indicator <- function(coef_name, moderator, mod_value, type_mod) {
  parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
  mod_part <- parts[grepl(paste0("V", moderator, "."), parts, fixed = TRUE)]
  if (!length(mod_part)) return(0)
  stripped <- sub("^V", "", mod_part[1])
  if (type_mod == "c") {
    dot_parts <- strsplit(stripped, ".", fixed = TRUE)[[1]]
    if (length(dot_parts) < 2 || dot_parts[2] == "") return(mod_value)
    if (!grepl("^[0-9]+$", dot_parts[2])) return(0)
    lev <- as.integer(dot_parts[2])
    return(if (identical(lev, as.integer(mod_value))) 1 else 0)
  }
  mod_value
}

# Map a user-supplied moderator value onto the internal scale the fit uses:
# the 0-based level code for a categorical moderator (matched against the
# moderator's ORIGINAL levels, so a 1/2-coded or labelled moderator is read in
# the user's own coding), the z value (value - mean) / sd for a continuous one
# (the moderator is standardized internally, mgm's scale = TRUE). Anything that
# is not one of those is refused by class.
#' @noRd
.mmg_moderator_value <- function(object, value) {
  mod <- object$moderator
  mod_label <- object$labels[mod]
  if (length(value) != 1L || is.na(value)) {
    .psn_stop(sprintf(
      "`value` must be a single non-missing moderator value for '%s'.",
      mod_label))
  }
  if (object$types[mod] == "c") {
    levs <- object$moderator_levels
    key <- if (is.factor(value) || is.logical(value)) as.character(value) else value
    if (!is.character(key) && !is.numeric(key)) {
      .psn_stop("`value` must be a number or a level label.")
    }
    # factor() labels numeric levels with as.character(), so the same coercion
    # of a numeric value finds its level exactly.
    idx <- match(as.character(key), levs)
    if (is.na(idx)) {
      .psn_stop(sprintf(paste0(
        "`value` = %s is not a level of the categorical moderator '%s'; ",
        "its levels are: %s."), format(value), mod_label,
        paste(levs, collapse = ", ")))
    }
    return(idx - 1L)
  }
  if (!is.numeric(value) || !is.finite(value)) {
    .psn_stop(sprintf(paste0(
      "`value` for the continuous moderator '%s' must be a single finite ",
      "number on its original scale."), mod_label))
  }
  (value - object$moderator_center) / object$moderator_scale
}

# Apply the requested mgm threshold to a coefficient vector.
#' @noRd
.mmg_threshold <- function(b, threshold, npar, n) {
  if (threshold == "none") return(b)
  tau <- if (threshold == "LW") {
    sqrt(2L) * sqrt(sum(b^2)) * sqrt(log(npar) / n)
  } else {
    sqrt(log(npar) / n)
  }
  b[abs(b) < tau] <- 0
  b
}

# One nodewise moderated fit on the BASE kernel (penalized IRLS + coordinate
# descent, R/lasso_glm.R) -- no glmnet, no compiled code.
#
# The design matrix X and the aggregation downstream are identical to the glmnet
# path; only the solver differs. Two equivalences make this a drop-in rather than
# an approximation:
#
#  * EBIC. The glmnet path minimises -2*LL + df*log(n) + 2*gamma*df*log(npar).
#    For a gaussian node -2*LL_lambda collapses to RSS - 2*LL_sat with LL_sat
#    constant along the path, and for a categorical node -2*LL IS the deviance,
#    so the argmin equals .nodewise_ebic()'s dev + df*log(n) + 2*gamma*df*log(p).
#    The two criteria select on the same objective.
#
#  * The 2-class multinomial convention. glmnet fits a binary categorical node
#    with family = "multinomial", whose two class coefficient vectors are
#    symmetric (-a, +a); the ordinary binary logit coefficient is b = 2a. The
#    aggregation helpers average |coef| over the class list, so returning
#    beta_list = list(-b/2, +b/2) reproduces mgm's directed parameter exactly.
#    The LW threshold is unaffected by the halving: thresholding a = b/2 at
#    tau(a) is identical to thresholding b at tau(b), because tau scales with the
#    vector's own norm.
#' @noRd
.mmg_fit_node_base <- function(X, y, is_cat, gamma, threshold, nlambda,
                               lambda_min_ratio) {
  npar <- ncol(X)
  n <- nrow(X)
  if (is_cat) {
    fit <- .nodewise_ebic(X, y - 1, "binomial", gamma, nlambda,
                          lambda_min_ratio)
    b <- .mmg_threshold(fit$beta, threshold, npar, n)
    # Split the logit coefficient into glmnet's symmetric 2-class convention.
    return(list(beta = NULL, beta_list = list(-b / 2, b / 2),
                col_names = colnames(X), multinomial = TRUE,
                lambda = fit$lambda, kkt = fit$kkt, converged = fit$converged,
                thresholded = any(fit$beta != 0 & b == 0)))
  }
  fit <- .nodewise_ebic(X, y, "gaussian", gamma, nlambda, lambda_min_ratio)
  b <- .mmg_threshold(fit$beta, threshold, npar, n)
  list(beta = b, col_names = colnames(X), multinomial = FALSE,
       lambda = fit$lambda, kkt = fit$kkt, converged = fit$converged,
       thresholded = any(fit$beta != 0 & b == 0))
}

# Fit the moderated MGM (nodewise lasso + EBIC + LW threshold + AND/OR
# 2-side main / 3-side interaction aliveness). Returns a psychnet_moderated.
#
# engine = "base"   -- psychnets' own penalized-IRLS kernel (default). Gaussian
#                      and binary nodes; every nodewise fit is KKT-certified.
# engine = "glmnet" -- the reference path (native = FALSE). Additionally supports
#                      multi-level categorical nodes via multinomial glmnet.
#' @noRd
.mmg_estimate <- function(mat, types, moderator, gamma, rule, threshold,
                          labels, level = NULL, engine = "base",
                          nlambda = 100L, lambda_min_ratio = 0.01,
                          moderator_levels = NULL) {
  if (engine == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
    .psn_stop("native = FALSE needs the 'glmnet' package; install it, or keep ",
         "native = TRUE (the default, which uses the base-R kernel).",
         call. = FALSE)
  }
  stopifnot("internal: `moderator` is one column index" =
              length(moderator) == 1L && moderator >= 1 && moderator <= ncol(mat))
  if (ncol(mat) < 3L)
    .psn_stop("Moderated MGM requires at least three variables (two endpoints and a moderator).",
         call. = FALSE)
  # The base kernel is a binomial solver: it covers gaussian + binary nodes, the
  # documented mgm_fit() scope. A multi-level categorical needs the multinomial
  # solver, which only the glmnet engine has.
  if (engine == "base") {
    n_lev <- vapply(which(types == "c"),
                    function(j) length(unique(mat[, j])), integer(1))
    if (length(n_lev) && any(n_lev > 2L)) {
      bad <- which(types == "c")[n_lev > 2L]
      .psn_stop(sprintf(paste0(
        "Column(s) %s are categorical with more than 2 levels; the base kernel ",
        "supports gaussian and binary (0/1) nodes only. One-hot encode them, or ",
        "pass native = FALSE to use the multinomial glmnet engine."),
        paste(colnames(mat)[bad], collapse = ", ")), call. = FALSE)
    }
  }
  data <- as.data.frame(mat)
  p <- ncol(data)
  scale_on <- any(types == "g")
  # The moderator's own coding, kept so condition() can read a value on the
  # user's scale: original level labels (in code order) for a categorical
  # moderator, the mean/sd it is standardized by for a continuous one.
  if (types[moderator] == "c") {
    moderator_levels <- moderator_levels %||%
      as.character(sort(unique(mat[, moderator])))
    moderator_center <- NA_real_; moderator_scale <- NA_real_
  } else {
    moderator_levels <- NULL
    moderator_center <- mean(mat[, moderator])
    moderator_scale <- stats::sd(mat[, moderator])
  }
  # Standardize continuous columns to unit variance up front (mgm scale = TRUE):
  # without this the gaussian main effects and edges carry each variable's SD and
  # diverge from mgm. Must happen here -- the moderated path bypasses mgm_fit's
  # own response scaling.
  if (scale_on)
    for (i in which(types == "g")) data[[i]] <- as.numeric(scale(data[[i]]))
  for (i in which(types == "c")) data[[i]] <- as.factor(data[[i]])
  colnames(data) <- paste0("V", seq_len(p), ".")
  mod_name <- colnames(data)[moderator]
  n <- nrow(data)

  fits <- lapply(seq_len(p), function(v) {
    is_mod <- (v == moderator)
    main_rhs <- paste(colnames(data)[-v], collapse = " + ")
    if (is_mod) {
      other <- colnames(data)[-v]
      int_terms <- vapply(utils::combn(other, 2, simplify = FALSE),
                          function(pr) paste(pr, collapse = "*"), character(1))
    } else {
      nonmod_preds <- setdiff(colnames(data)[-v], mod_name)
      int_terms <- paste0(nonmod_preds, "*", mod_name)
    }
    form <- stats::as.formula(paste(colnames(data)[v], "~", main_rhs, "+",
                                    paste(int_terms, collapse = " + ")))
    X <- stats::model.matrix(form, data = data)[, -1, drop = FALSE]

    # Standardize interaction columns whose involved variables are ALL Gaussian,
    # after the model matrix is built (mgm main lines 216-226).
    if (scale_on) {
      l_split <- strsplit(colnames(X), ":", fixed = TRUE)
      all_gauss <- vapply(l_split, function(parts) {
        nums <- as.numeric(sub("\\.[0-9]*$", "", sub("^V", "", parts)))
        all(!is.na(nums) & types[nums] == "g")
      }, logical(1))
      if (any(all_gauss))
        X[, all_gauss] <- apply(X[, all_gauss, drop = FALSE], 2, scale)
    }

    npar <- ncol(X)
    y <- as.numeric(data[[v]])

    if (engine == "base") {
      return(.mmg_fit_node_base(X, y, is_cat = types[v] == "c", gamma = gamma,
                                threshold = threshold, nlambda = nlambda,
                                lambda_min_ratio = lambda_min_ratio))
    }

    if (types[v] == "c") {
      glm_args <- list(x = X, y = y, family = "multinomial", alpha = 1,
                       nlambda = nlambda, intercept = TRUE)
      if (!is.null(lambda_min_ratio))
        glm_args$lambda.min.ratio <- lambda_min_ratio
      fit <- do.call(glmnet::glmnet, glm_args)
      beta_list <- lapply(fit$beta, as.matrix)
      nz <- Reduce("+", lapply(beta_list, function(B) (B != 0) * 1)) > 0
      n_nb <- colSums(nz)
      tab <- tabulate(y, nbins = max(y)); pj <- tab / n
      LL_null <- n * sum(pj[pj > 0] * log(pj[pj > 0]))
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      LL_lambda <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
      EBIC <- -2 * LL_lambda + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
      idx <- which.min(EBIC)
      beta_raw <- lapply(beta_list, function(B) B[, idx])
      beta_sel <- lapply(beta_raw, .mmg_threshold, threshold = threshold,
                         npar = npar, n = n)
      return(list(beta = NULL, beta_list = beta_sel, col_names = colnames(X),
                  multinomial = TRUE, lambda = fit$lambda[idx],
                  thresholded = any(unlist(beta_raw) != 0 &
                                      unlist(beta_sel) == 0)))
    }

    glm_args <- list(x = X, y = y, family = "gaussian", alpha = 1,
                     nlambda = nlambda, intercept = TRUE)
    if (!is.null(lambda_min_ratio))
      glm_args$lambda.min.ratio <- lambda_min_ratio
    fit <- do.call(glmnet::glmnet, glm_args)
    beta_path <- as.matrix(fit$beta)
    n_nb <- colSums(beta_path != 0)
    LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y))^2)) + 1)
    LL_sat  <- 0.5 * fit$nulldev + LL_null
    LL_lambda <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
    EBIC <- -2 * LL_lambda + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
    idx <- which.min(EBIC)
    b_raw <- beta_path[, idx]
    b <- .mmg_threshold(b_raw, threshold, npar, n)
    list(beta = b, col_names = colnames(X), multinomial = FALSE,
         lambda = fit$lambda[idx], thresholded = any(b_raw != 0 & b == 0))
  })

  pw_alive <- matrix(FALSE, p, p)
  int_alive <- matrix(FALSE, p, p)
  if (rule == "AND") {
    all_pairs <- utils::combn(p, 2)
    for (k in seq_len(ncol(all_pairs))) {
      i <- all_pairs[1, k]; j <- all_pairs[2, k]
      pw_alive[i, j] <- pw_alive[j, i] <-
        (.mmg_main_mag(fits[[i]], j) > 0 && .mmg_main_mag(fits[[j]], i) > 0)
    }
    non_mod <- setdiff(seq_len(p), moderator)
    if (length(non_mod) >= 2L) {
      nm_pairs <- utils::combn(non_mod, 2)
      for (k in seq_len(ncol(nm_pairs))) {
        i <- nm_pairs[1, k]; j <- nm_pairs[2, k]
        int_alive[i, j] <- int_alive[j, i] <-
          (.mmg_int_mag(fits[[i]], j, moderator) > 0 &&
           .mmg_int_mag(fits[[j]], i, moderator) > 0 &&
           .mmg_int_mag(fits[[moderator]], i, j) > 0)
      }
    }
  } else {
    pw_alive[] <- TRUE
    non_mod <- setdiff(seq_len(p), moderator)
    int_alive[non_mod, non_mod] <- TRUE
  }

  # Worst nodewise KKT residual: the base kernel certifies every nodewise fit
  # against its own convex objective. The glmnet path carries no such residual.
  # The LW/HW threshold edits the penalized solution afterwards, so when it
  # removed a coefficient the network is not the certified optimum: `kkt` is
  # NA and the residual of the penalized fits is kept as `fit_kkt`.
  fit_kkt <- if (engine == "base")
    max(vapply(fits, function(f) f$kkt %||% NA_real_, numeric(1))) else NA_real_
  thresholded <- any(vapply(fits, function(f) isTRUE(f$thresholded), logical(1)))
  if (engine == "base") .warn_nodewise_convergence(fits, labels, "mgm_fit")

  structure(list(fits = fits, moderator = moderator, p = p, types = types,
                 level = level, labels = labels, n = n, data = mat,
                 pw_alive = pw_alive, int_alive = int_alive,
                 kkt = if (thresholded) NA_real_ else fit_kkt,
                 fit_kkt = fit_kkt,
                 certificate_target = if (thresholded) "pre_threshold_fit" else
                   "returned_network",
                 moderator_levels = moderator_levels,
                 moderator_center = moderator_center,
                 moderator_scale = moderator_scale,
                 engine = engine,
                 params = list(gamma = gamma, rule = rule, threshold = threshold,
                               nlambda = nlambda,
                               lambda_min_ratio = lambda_min_ratio)),
            class = "psychnet_moderated")
}

# Conditioned p x p weight matrix at a moderator value (mgm::condition()).
#' @noRd
.mmg_condition_wadj <- function(fit, mod_value, rule) {
  fits <- fit$fits; p <- fit$p; moderator <- fit$moderator
  pw_alive <- fit$pw_alive; int_alive <- fit$int_alive
  type_mod <- fit$types[moderator]
  side_mag <- matrix(0, p, p)

  for (v in seq_len(p)) {
    if (v == moderator) next
    nf <- fits[[v]]; cn <- nf$col_names
    for (j in seq_len(p)) {
      if (j == v || j == moderator) next
      m_idx <- which(!grepl(":", cn, fixed = TRUE) &
                       grepl(paste0("V", j, "."), cn, fixed = TRUE))
      has_int <- grepl(":", cn, fixed = TRUE)
      i_idx <- which(grepl(paste0("V", j, "."), cn, fixed = TRUE) & has_int &
                     grepl(paste0("V", moderator, "."), cn, fixed = TRUE))

      .cond_dummies <- function(beta_vec) {
        vapply(m_idx, function(d) {
          val <- beta_vec[d]
          if (!pw_alive[v, j]) val <- 0
          if (length(i_idx) && int_alive[v, j]) {
            dn <- gsub(".", "\\.", cn[d], fixed = TRUE)
            assoc <- i_idx[grepl(paste0("^", dn, ":"), cn[i_idx]) |
                             grepl(paste0(":", dn, "$"), cn[i_idx])]
            for (ii in assoc) {
              ind <- .mmg_cond_indicator(cn[ii], moderator, mod_value, type_mod)
              val <- val + beta_vec[ii] * ind
            }
          }
          val
        }, numeric(1))
      }

      if (!nf$multinomial) {
        cvals <- .cond_dummies(nf$beta)
        side_mag[v, j] <- if (length(cvals)) mean(abs(cvals)) else 0
      } else {
        per_class <- vapply(nf$beta_list, function(b) {
          cvals <- .cond_dummies(b)
          if (length(cvals)) mean(abs(cvals)) else 0
        }, numeric(1))
        side_mag[v, j] <- mean(per_class)
      }
    }
  }

  wadj <- matrix(0, p, p)
  pairs <- utils::combn(p, 2)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    if (i == moderator || j == moderator) next
    m_par <- c(side_mag[i, j], side_mag[j, i])
    edge <- if (rule == "AND") { if (any(m_par == 0)) 0 else mean(m_par) }
            else mean(m_par)
    wadj[i, j] <- wadj[j, i] <- edge
  }
  wadj
}

#' Condition a moderated network at a moderator value
#'
#' Extracts the effective pairwise network implied by a moderated MGM
#' ([mgm_fit()] with `moderators`) at a given value of the moderator, mirroring
#' `mgm::condition()`: it applies the AND-rule pre-filter, absorbs the moderator
#' value into the main-effect coefficients, and re-aggregates the pairwise edges.
#'
#' `value` is always read on the moderator's **own** scale. For a categorical
#' moderator it is one of the levels observed in the data, in the coding the
#' data used: `1`/`2` for a moderator coded 1/2, `"treat"` for a factor level
#' `"treat"` (the internal 0-based level code is never exposed). For a
#' continuous moderator it is a value on the raw scale; the moderator is
#' standardized internally (mgm's `scale = TRUE`), and the value is mapped to
#' `(value - mean) / sd` before conditioning, so `value = mean(moderator)` gives
#' the network at the average moderator. (`mgm::condition()` instead takes the
#' value on its internal standardized scale: `condition(fit, v)` here equals
#' `mgm::condition(ref, list(k = (v - mean) / sd))`.)
#'
#' @param object A `psychnet_moderated` object from `mgm_fit(..., moderators=)`.
#' @param value A single moderator value to condition on: a level of a
#'   categorical moderator (number, label, factor or logical, as coded in the
#'   data), or a finite number on the raw scale of a continuous one.
#' @param rule Symmetrization rule; defaults to the rule used at fit time.
#' @param x A `psychnet_moderated` object.
#' @param row.names,optional Standard [as.data.frame()] arguments; ignored.
#' @param ... Unused.
#' @return `condition()` returns a `psychnet` network object (the moderator node
#'   carries no edges), with `$mod_value` the value conditioned on.
#'
#'   `as.data.frame()` on the moderated fit returns, with `value = NULL`, one
#'   row per pair of nodes joined by a main effect or a moderation effect, with
#'   columns `from`, `to`, `main_effect` (aggregated magnitude of the pairwise
#'   main effect, mgm convention: mean absolute nodewise estimate, zero when the
#'   pair does not survive `rule`), `interaction` (aggregated magnitude of the
#'   moderator x pair three-way interaction, `NA` for a pair that includes the
#'   moderator itself), `moderated` (logical: `interaction > 0`) and
#'   `involves_moderator`. With a `value`, it returns the tidy edge list
#'   (`from`, `to`, `weight`) of `condition(x, value)`.
#'
#'   `print()` on the moderated fit shows the number of nodes, the moderator
#'   and the fit settings, and returns `x` invisibly.
#'
#'   Errors of class `psychnet_bad_input` when `value` is missing, not length
#'   one, not a level of a categorical moderator, or not a finite number for a
#'   continuous one.
#' @examples
#' set.seed(1)
#' x1 <- stats::rnorm(400); x2 <- stats::rnorm(400)
#' mod <- rep(1:2, each = 200)                 # moderator coded 1/2
#' y <- x1 * (mod == 2) + stats::rnorm(400)   # x1-y edge only when mod == 2
#' d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
#' fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = "mod")
#' fit
#' as.data.frame(fit)            # main and moderation effects per edge
#' condition(fit, value = 1)     # no x1-y edge
#' condition(fit, value = 2)     # the x1-y edge appears
#' as.data.frame(fit, value = 2) # the same network as a tidy edge list
#' @export
condition <- function(object, value, rule = NULL) {
  if (!inherits(object, "psychnet_moderated")) {
    .psn_stop("`object` must be a moderated fit from mgm_fit(..., moderators = ).")
  }
  if (is.null(rule)) rule <- object$params$rule
  rule <- .psn_match_arg(rule, c("AND", "OR"))
  internal <- .mmg_moderator_value(object, value)
  wadj <- .mmg_condition_wadj(object, internal, rule)
  .new_psychnet(wadj, object$labels, method = "mgm_moderated",
                directed = FALSE, n_obs = object$n, data = object$data,
                extra = list(types = stats::setNames(object$types, object$labels),
                             moderator = object$moderator, mod_value = value,
                             rule = rule))
}

#' @describeIn condition Tidy table of the moderated fit: per-edge main and
#'   moderation effects, or (with `value`) the conditioned edge list.
#' @export
as.data.frame.psychnet_moderated <- function(x, row.names = NULL,
                                             optional = FALSE, ...,
                                             value = NULL) {
  if (!is.null(value)) return(as.data.frame(condition(x, value)))
  p <- x$p; mod <- x$moderator; fits <- x$fits
  pairs <- utils::combn(p, 2)
  rows <- lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]; j <- pairs[2, k]
    side <- c(.mmg_main_mag(fits[[i]], j), .mmg_main_mag(fits[[j]], i))
    main <- if (x$pw_alive[i, j]) mean(side) else 0
    with_mod <- i == mod || j == mod
    inter <- if (with_mod) NA_real_ else if (x$int_alive[i, j]) {
      mean(c(.mmg_int_mag(fits[[i]], j, mod), .mmg_int_mag(fits[[j]], i, mod),
             .mmg_int_mag(fits[[mod]], i, j)))
    } else 0
    data.frame(from = x$labels[i], to = x$labels[j], main_effect = main,
               interaction = inter, moderated = if (with_mod) NA else inter > 0,
               involves_moderator = with_mod, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  keep <- out$main_effect > 0 | (!is.na(out$interaction) & out$interaction > 0)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @describeIn condition The moderated fit: nodes, moderator and settings.
#' @export
print.psychnet_moderated <- function(x, ...) {
  cat(sprintf("<psychnet_moderated> %d nodes, moderator = %s\n",
              x$p, x$labels[x$moderator]))
  cat(sprintf("  gamma = %.3g, rule = %s, threshold = %s\n",
              x$params$gamma, x$params$rule, x$params$threshold))
  if (identical(x$types[x$moderator], "c") && length(x$moderator_levels)) {
    cat(sprintf("  moderator levels: %s\n",
                paste(x$moderator_levels, collapse = ", ")))
  }
  cat("  as.data.frame(fit) for main and moderation effects;",
      "condition(fit, value) for the network at a moderator value.\n")
  invisible(x)
}
