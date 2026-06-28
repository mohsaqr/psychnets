# Moderated Mixed Graphical Model (Haslbeck "group differences via moderation").
# Ported verbatim in its numerical logic from Nestimate's verified implementation
# (which matches mgm::mgm(..., moderators = k) + mgm::condition() to machine
# precision); only the argument names are adapted to psychnet conventions
# (gamma, rule, types) and the outputs are wrapped as psychnet objects. A chosen
# variable moderates every pairwise edge; group-structure differences appear as
# nonzero moderator x endpoint interactions. The fit is glmnet-based (the
# reference path); there is no base engine for the moderated model.

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
#' @noRd
.mmg_cond_indicator <- function(coef_name, moderator, mod_value, type_mod) {
  parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
  mod_part <- parts[grepl(paste0("V", moderator, "."), parts, fixed = TRUE)]
  if (!length(mod_part)) return(0)
  stripped <- sub("^V", "", mod_part[1])
  if (type_mod == "c") {
    dot_parts <- strsplit(stripped, ".", fixed = TRUE)[[1]]
    if (length(dot_parts) < 2 || dot_parts[2] == "") return(mod_value)
    lev <- suppressWarnings(as.numeric(dot_parts[2]))
    if (is.na(lev)) return(0)
    return(if (lev == mod_value) 1 else 0)
  }
  mod_value
}

# Fit the moderated MGM (glmnet nodewise lasso + EBIC + LW threshold + AND/OR
# 2-side main / 3-side interaction aliveness). Returns a psychnet_moderated.
#' @noRd
.mmg_estimate <- function(mat, types, moderator, gamma, rule, threshold,
                          labels, level = NULL) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE),
            length(moderator) == 1L, moderator >= 1, moderator <= ncol(mat))
  data <- as.data.frame(mat)
  p <- ncol(data)
  scale_on <- any(types == "g")
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

    if (types[v] == "c") {
      fit <- glmnet::glmnet(X, y, family = "multinomial", alpha = 1,
                            intercept = TRUE)
      beta_list <- lapply(fit$beta, as.matrix)
      nz <- Reduce("+", lapply(beta_list, function(B) (B != 0) * 1)) > 0
      n_nb <- colSums(nz)
      tab <- tabulate(y, nbins = max(y)); pj <- tab / n
      LL_null <- n * sum(pj[pj > 0] * log(pj[pj > 0]))
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      LL_lambda <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
      EBIC <- -2 * LL_lambda + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
      idx <- which.min(EBIC)
      beta_sel <- lapply(beta_list, function(B) B[, idx])
      if (threshold == "LW") {
        beta_sel <- lapply(beta_sel, function(bb) {
          tau <- sqrt(2L) * sqrt(sum(bb^2)) * sqrt(log(npar) / n)
          bb[abs(bb) < tau] <- 0; bb
        })
      }
      return(list(beta = NULL, beta_list = beta_sel, col_names = colnames(X),
                  multinomial = TRUE, lambda = fit$lambda[idx]))
    }

    fit <- glmnet::glmnet(X, y, family = "gaussian", alpha = 1, intercept = TRUE)
    beta_path <- as.matrix(fit$beta)
    n_nb <- colSums(beta_path != 0)
    LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y))^2)) + 1)
    LL_sat  <- 0.5 * fit$nulldev + LL_null
    LL_lambda <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
    EBIC <- -2 * LL_lambda + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
    idx <- which.min(EBIC)
    b <- beta_path[, idx]
    if (threshold == "LW") {
      tau <- sqrt(2L) * sqrt(sum(b^2)) * sqrt(log(npar) / n)
      b[abs(b) < tau] <- 0
    }
    list(beta = b, col_names = colnames(X), multinomial = FALSE,
         lambda = fit$lambda[idx])
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

  structure(list(fits = fits, moderator = moderator, p = p, types = types,
                 level = level, labels = labels, n = n, data = mat,
                 pw_alive = pw_alive, int_alive = int_alive,
                 params = list(gamma = gamma, rule = rule, threshold = threshold)),
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
#' @param object A `psychnet_moderated` object from `mgm_fit(..., moderators=)`.
#' @param value Moderator value to condition on (e.g. `0` or `1` for a binary
#'   moderator, or any numeric for a continuous one).
#' @param rule Symmetrization rule; defaults to the rule used at fit time.
#' @return A `psychnet` network object (the moderator node carries no edges).
#' @examples
#' set.seed(1)
#' x1 <- stats::rnorm(400); x2 <- stats::rnorm(400)
#' mod <- rep(0:1, each = 200)
#' y <- x1 * (mod == 1) + stats::rnorm(400)   # x1-y edge only when mod == 1
#' d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
#' fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4)
#' condition(fit, value = 1)
#' @export
condition <- function(object, value, rule = NULL) {
  stopifnot(inherits(object, "psychnet_moderated"))
  if (is.null(rule)) rule <- object$params$rule
  rule <- match.arg(rule, c("AND", "OR"))
  wadj <- .mmg_condition_wadj(object, value, rule)
  .new_psychnet(wadj, object$labels, method = "mgm_moderated",
                directed = FALSE, n_obs = object$n, data = object$data,
                extra = list(types = stats::setNames(object$types, object$labels),
                             moderator = object$moderator, mod_value = value,
                             rule = rule))
}

#' Print a moderated MGM fit
#'
#' @param x A `psychnet_moderated` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_moderated <- function(x, ...) {
  cat(sprintf("<psychnet_moderated> %d nodes, moderator = %s\n",
              x$p, x$labels[x$moderator]))
  cat(sprintf("  gamma = %.3g, rule = %s, threshold = %s\n",
              x$params$gamma, x$params$rule, x$params$threshold))
  cat("  condition(fit, value) to extract the network at a moderator value.\n")
  invisible(x)
}
