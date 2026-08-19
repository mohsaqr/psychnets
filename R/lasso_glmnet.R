# Optional glmnet solver engine for the nodewise estimators (ising_fit, mgm_fit).
#
# The default engine = "base" is the pure-R, dependency-free, self-certified
# nodewise lasso in R/lasso_glm.R (internally optimal, KKT ~1e-9). It agrees
# with IsingFit / mgm on structure and sign, but selects a slightly different
# EBIC lambda along an independent path, so its weights are not byte-identical
# to those packages. engine = "glmnet" swaps the per-node solve to
# glmnet::glmnet() with the reference's exact lambda path, EBIC formula, df
# count, thresholding, and AND/OR symmetrization, so the returned network
# byte-matches IsingFit::IsingFit() / mgm::mgm() (to ~1e-16). This mirrors the
# engine = "glasso" option on the GGM side (R/glasso.R): an opt-in Suggests
# dependency that reproduces the reference solver exactly, at the cost of the
# self-certification guarantee (the certificate then reflects glmnet's own
# tolerance, not the base path's certified optimum).

# The public `native` switch (TRUE = pure-R FHT-2010 solver, default; FALSE =
# delegate each per-node fit to glmnet for byte-identical IsingFit / mgm parity)
# is resolved to the internal engine name by .resolve_native() in glasso.R.

# Numerically stable log(1 + exp(x)).
#' @noRd
.log1pexp <- function(x) {
  out <- x
  small <- x < 0
  out[small]  <- log1p(exp(x[small]))
  out[!small] <- x[!small] + log1p(exp(-x[!small]))
  out
}

# Self-certify a glmnet lasso fit without knowing its internal lambda scale.
# A lasso optimum has all ACTIVE gradients equal in magnitude (= the effective
# penalty) and all INACTIVE gradients no larger. We recover that penalty as the
# mean active |gradient| on the package's own standardized predictor scale and
# report the worst stationarity violation -- near-zero certifies the fit solves
# its lasso objective, regardless of glmnet's reported lambda units.
#' @noRd
.self_lambda_kkt <- function(X, y, b0, beta, family) {
  std <- .standardize(X)
  Xs  <- std$X
  beta_std <- beta * std$scale
  n <- nrow(X)
  eta <- b0 + as.numeric(X %*% beta)
  mu  <- if (family == "binomial") 1 / (1 + exp(-eta)) else eta
  grad <- as.numeric(crossprod(Xs, y - mu)) / n
  active <- abs(beta_std) > 1e-8
  # With an empty active set the penalty cannot be read off the active
  # gradients; the smallest lambda consistent with the empty solution is
  # max|grad| (glmnet's own lambda_max), at which the empty model is exactly
  # optimal. Reporting lam_eff = 0 there would grade the whole gradient as a
  # violation of a fit that is in fact optimal.
  lam_eff <- if (any(active)) mean(abs(grad[active])) else max(abs(grad))
  v_0 <- abs(sum(y - mu) / n)
  v_a <- if (any(active))  max(abs(grad[active] - lam_eff * sign(beta_std[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(grad[!active]) - lam_eff, 0)) else 0
  max(v_0, v_a, v_i)
}

# Multinomial counterpart of .self_lambda_kkt for a k-class glmnet fit with the
# (default, ungrouped) per-class lasso. The stationarity condition of the
# multinomial objective uses the SOFTMAX class probabilities
#   p_cl = exp(eta_cl) / sum_m exp(eta_m),
# not a per-class sigmoid: grading each class as an independent logistic
# regression evaluates a different objective and reports a residual of ~0.2 for
# a fit glmnet has in fact converged (true residual ~1e-5). `a0` is the length-k
# intercept vector, `beta_list` the k per-class coefficient vectors.
#' @noRd
.self_lambda_kkt_multinomial <- function(X, y, a0, beta_list) {
  std <- .standardize(X)
  Xs  <- std$X
  n <- nrow(X); k <- length(beta_list)
  B <- do.call(cbind, lapply(beta_list, as.numeric))          # npar x k
  B_std <- B * std$scale
  eta <- sweep(X %*% B, 2L, a0, "+")
  eta <- eta - apply(eta, 1L, max)                            # overflow guard
  P <- exp(eta); P <- P / rowSums(P)
  Y <- outer(as.integer(y), seq_len(k), "==") * 1
  grad <- crossprod(Xs, Y - P) / n
  active <- abs(B_std) > 1e-8
  lam_eff <- if (any(active)) mean(abs(grad[active])) else max(abs(grad))  # see .self_lambda_kkt
  v_0 <- max(abs(colSums(Y - P) / n))
  v_a <- if (any(active))  max(abs(grad[active] - lam_eff * sign(B_std[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(grad[!active]) - lam_eff, 0)) else 0
  max(v_0, v_a, v_i)
}

# One nodewise glmnet fit (single-response: gaussian or binomial) with EBIC
# selection on glmnet's default lambda path. Returns the same fields the
# nodewise callers consume from .nodewise_ebic(), plus the self-lambda KKT.
#
# EBIC convention matches the references exactly:
#   binomial (IsingFit): -2*loglik + J*log(n) + 2*gamma*J*log(p_pred)
#   gaussian (mgm)     : -2*LL    + J*log(n) + 2*gamma*J*log(p_pred)
# where J = #nonzero coefficients and p_pred = number of predictor columns.
#' @noRd
.nodewise_glmnet <- function(X, y, family, gamma, p_pred, nlambda,
                             lambda_min_ratio) {
  n <- nrow(X)
  fam <- if (family == "binomial") "binomial" else "gaussian"
  glm_args <- list(x = X, y = y, family = fam, alpha = 1, nlambda = nlambda,
                   intercept = TRUE, standardize = TRUE)
  if (!is.null(lambda_min_ratio))
    glm_args$lambda.min.ratio <- lambda_min_ratio
  fit <- do.call(glmnet::glmnet, glm_args)
  beta_path <- as.matrix(fit$beta)
  J <- colSums(beta_path != 0)

  if (fam == "binomial") {
    ll <- vapply(seq_along(fit$lambda), function(k) {
      eta <- fit$a0[k] + as.numeric(X %*% beta_path[, k])
      sum(y * eta - .log1pexp(eta))
    }, numeric(1))
    ebic <- -2 * ll + J * log(n) + 2 * gamma * J * log(p_pred)
  } else {
    LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y))^2)) + 1)
    LL_sat  <- 0.5 * fit$nulldev + LL_null
    LL      <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
    ebic <- -2 * LL + J * log(n) + 2 * gamma * J * log(p_pred)
  }

  idx  <- which.min(ebic)
  b0   <- fit$a0[idx]
  beta <- as.numeric(beta_path[, idx])
  kkt  <- .self_lambda_kkt(X, y, b0, beta, fam)
  list(b0 = b0, beta = beta, lambda = fit$lambda[idx], ebic = ebic[idx],
       kkt = kkt)
}

# --- Ising via glmnet (byte-matches IsingFit::IsingFit) -----------------------
# Per-node L1-penalized logistic regression on the raw 0/1 predictors with the
# IsingFit EBIC, AND/OR symmetrization, and raw-scale node thresholds (the a0
# intercepts). Predictors are 0/1, so the stored nodewise center/scale are
# 0/1 and beta_std holds the raw coefficients -- giving net_predict() the
# correct logistic linear predictor with no transform.
#' @noRd
.ising_fit_glmnet <- function(mat, gamma, rule, nlambda, lambda_min_ratio,
                              labels) {
  n <- nrow(mat); p <- ncol(mat)
  fits <- lapply(seq_len(p), function(i) {
    .nodewise_glmnet(mat[, -i, drop = FALSE], mat[, i], "binomial",
                     gamma, p_pred = p - 1L, nlambda = nlambda,
                     lambda_min_ratio = lambda_min_ratio)
  })
  B <- matrix(0, p, p, dimnames = list(labels, labels))
  b0_raw <- vapply(fits, function(f) f$b0, numeric(1))
  for (i in seq_len(p)) B[i, -i] <- fits[[i]]$beta
  worst_kkt <- max(vapply(fits, function(f) f$kkt, numeric(1)))

  present <- if (rule == "AND") (B != 0) & (t(B) != 0) else (B != 0) | (t(B) != 0)
  W <- (B + t(B)) / 2
  W[!present] <- 0
  diag(W) <- 0

  .new_psychnet(W, labels, method = "ising", directed = FALSE,
                n_obs = n, data = mat,
                extra = list(thresholds = stats::setNames(b0_raw, labels),
                             rule = rule, kkt = worst_kkt, native = FALSE,
                             nlambda = nlambda,
                             lambda_min_ratio = lambda_min_ratio,
                             nodewise = list(intercept = b0_raw,
                                             beta_std = B,
                                             families = rep("binomial", p),
                                             center = rep(0, p),
                                             scale = rep(1, p))))
}

# --- MGM via glmnet (byte-matches mgm::mgm magnitudes) ------------------------
# Ports the mgm main-effects path: continuous columns scaled to unit variance,
# binary columns entered as 0/1 model-matrix dummies, gaussian responses fit by
# a gaussian lasso and binary responses by a 2-class multinomial lasso, EBIC
# selection with npar = #predictor columns, LW/HW thresholding, and AND/OR
# magnitude symmetrization. The reported edge magnitude byte-matches
# mgm$pairwise$wadj; the stored sign is recovered from a gaussian endpoint
# where one exists (a binary-binary edge sign is undefined, stored positive).
# `df` is the factor-typed frame REBUILT from the imputed numeric matrix `mat`
# by mgm_fit(), so both views share rows and values whatever `na_method` did.
#' @noRd
.mgm_fit_glmnet <- function(df, types, gamma, threshold, rule, nlambda,
                            lambda_min_ratio, labels, levels_v = NULL,
                            mat = NULL) {
  n <- nrow(df); p <- ncol(df)
  gix <- which(types == "g")
  if (is.null(levels_v)) levels_v <- ifelse(types == "c", 2L, 1L)
  if (is.null(mat)) mat <- as.matrix(.mgm_frame_numeric(df))
  has_multi <- any(levels_v > 2L)

  # Column-wise center/scale mapping raw data to the glmnet predictor scale:
  # gaussian columns are standardized (mean, sample sd, matching base scale());
  # categorical columns enter as model-matrix level dummies.
  #
  # Scaling is applied to the whole frame ONCE, up front, before any model
  # matrix is built -- so a gaussian node is standardized when it is the
  # response too, not only when it is a predictor. mgm::mgm scales at the top of
  # its main function; scaling per-node instead diverges.
  resp_center <- numeric(p); resp_scale <- rep(1, p)
  if (length(gix)) {
    resp_center[gix] <- colMeans(mat[, gix, drop = FALSE])
    resp_scale[gix]  <- pmax(apply(mat[, gix, drop = FALSE], 2L, stats::sd), 1e-12)
    for (j in gix) df[[j]] <- (mat[, j] - resp_center[j]) / resp_scale[j]
  }

  # Per-node fit. Returns the selected (unthresholded) coefficients on the
  # scaled-predictor scale, a signed per-predictor side value (for sign
  # recovery), the magnitude side value (for the graph), an effective-binomial
  # representation for prediction, and the self-lambda KKT.
  node_fit <- function(v) {
    # A k-level factor contributes k-1 dummy columns; `assign` maps each design
    # column back to its source variable, because the pairwise magnitudes are
    # aggregated over a variable's COLUMNS, not over single coefficients.
    mm  <- stats::model.matrix(~ ., data = df[, -v, drop = FALSE])
    asg <- attr(mm, "assign")[-1]
    X   <- mm[, -1, drop = FALSE]
    pred_var <- seq_len(p)[-v][asg]
    # npar is ncol(X) -- DUMMY COLUMNS, not variables. It enters both the EBIC
    # penalty and the LW threshold, so a 3-level factor inflates it by 2, not 1.
    npar <- ncol(X)

    if (types[v] == "c") {
      # Categorical responses use multinomial INCLUDING the binary case: a
      # 2-class multinomial, not family = "binomial".
      y <- factor(df[[v]])
      glm_args <- list(x = X, y = y, family = "multinomial", alpha = 1,
                       nlambda = nlambda, intercept = TRUE)
      if (!is.null(lambda_min_ratio))
        glm_args$lambda.min.ratio <- lambda_min_ratio
      fit <- do.call(glmnet::glmnet, glm_args)
      beta_list <- lapply(fit$beta, as.matrix)                  # one per class
      # n_neighbors counts non-zero coefficient COLUMNS aggregated across
      # response classes, not distinct predictor variables.
      nz <- Reduce("+", lapply(beta_list, function(B) (B != 0) * 1)) > 0
      n_nb <- colSums(nz)
      tab  <- tabulate(as.integer(y), nbins = length(beta_list))
      pj   <- tab / n
      LL_null <- n * sum(pj[pj > 0] * log(pj[pj > 0]))
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      LL      <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
      ebic <- -2 * LL + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
      idx  <- which.min(ebic)
      beta_sel <- lapply(beta_list, function(B) B[, idx])       # per-class, unthresholded
      a0_sel   <- fit$a0[, idx]
      k <- length(beta_list)
      if (k == 2L) {
        # Effective 2-class binomial logit for the second (positive) class:
        # P(y = lev2) = sigma(eta2 - eta1). glmnet's class order is the factor
        # level order, so this is positional and works for any 2-level coding.
        beta_eff <- beta_sel[[2L]] - beta_sel[[1L]]
        b0_eff   <- a0_sel[2L] - a0_sel[1L]
        y01 <- as.numeric(as.integer(y) - 1L)
        kkt <- .self_lambda_kkt(X, y01, b0_eff, beta_eff, "binomial")
        fam <- "binomial"
      } else {
        # A k > 2 response has no single effective logit; there is nothing
        # honest to hand net_predict, so it is marked and refused there rather
        # than collapsed into a binomial that would predict wrongly. Its
        # certificate is the multinomial (softmax) stationarity residual.
        beta_eff <- numeric(npar); b0_eff <- 0
        kkt <- .self_lambda_kkt_multinomial(X, y, a0_sel, beta_sel)
        fam <- "multinomial"
      }
      list(multinomial = TRUE, pred_var = pred_var, npar = npar,
           beta_sel = beta_sel, beta_eff = beta_eff, b0_eff = b0_eff,
           family = fam, kkt = kkt)
    } else {
      y <- as.numeric(df[[v]])
      glm_args <- list(x = X, y = y, family = "gaussian", alpha = 1,
                       nlambda = nlambda, intercept = TRUE)
      if (!is.null(lambda_min_ratio))
        glm_args$lambda.min.ratio <- lambda_min_ratio
      fit <- do.call(glmnet::glmnet, glm_args)
      beta_path <- as.matrix(fit$beta)
      n_nb <- colSums(beta_path != 0)
      LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y))^2)) + 1)
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      LL      <- -0.5 * ((1 - fit$dev.ratio) * fit$nulldev) + LL_sat
      ebic <- -2 * LL + n_nb * log(n) + 2 * gamma * n_nb * log(npar)
      idx  <- which.min(ebic)
      b <- beta_path[, idx]
      kkt <- .self_lambda_kkt(X, y, fit$a0[idx], b, "gaussian")
      list(multinomial = FALSE, pred_var = pred_var, npar = npar,
           beta_sel = b, beta_eff = b, b0_eff = fit$a0[idx],
           family = "gaussian", kkt = kkt)
    }
  }
  fits <- lapply(seq_len(p), node_fit)

  # LW / HW post-selection thresholding of the graph-side coefficients (mgm).
  thr_side <- function(nf) {
    if (threshold == "none") return(nf$beta_sel)
    if (nf$multinomial) {
      lapply(nf$beta_sel, function(b) {
        tau <- if (threshold == "LW") sqrt(sum(b^2)) * sqrt(log(nf$npar) / n)
               else sqrt(log(nf$npar) / n)
        b[abs(b) < tau] <- 0; b
      })
    } else {
      b <- nf$beta_sel
      tau <- if (threshold == "LW") sqrt(sum(b^2)) * sqrt(log(nf$npar) / n)
             else sqrt(log(nf$npar) / n)
      b[abs(b) < tau] <- 0; b
    }
  }
  beta_thr <- lapply(fits, thr_side)

  # Signed magnitude of node `i`'s thresholded effect of predictor `j`:
  #   magnitude = mean(|coef columns for j|) (mgm side_mag),
  #   sign      = sign of the single coef for a gaussian node, NA otherwise.
  side <- function(i, j) {
    nf <- fits[[i]]; cols <- which(nf$pred_var == j)
    if (!length(cols)) return(c(mag = 0, sgn = NA_real_))
    if (nf$multinomial) {
      vals <- unlist(lapply(beta_thr[[i]], function(b) b[cols]))
      c(mag = mean(abs(vals)), sgn = NA_real_)
    } else {
      v <- beta_thr[[i]][cols]
      # A sign is only meaningful when the predictor is a single column; a
      # multi-dummy categorical has no signed direction (mgm leaves it undefined).
      sgn <- if (length(v) == 1L && v[1] != 0) sign(v[1]) else NA_real_
      c(mag = mean(abs(v)), sgn = sgn)
    }
  }

  W <- matrix(0, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      sij <- side(i, j); sji <- side(j, i)
      mags <- c(sij["mag"], sji["mag"])
      mag <- if (rule == "AND") { if (any(mags == 0)) 0 else mean(mags) } else mean(mags)
      if (mag == 0) next
      sgn <- c(sij["sgn"], sji["sgn"])
      sgn <- sgn[!is.na(sgn)]
      edge_sign <- if (length(sgn)) sign(sgn[1]) else 1
      W[i, j] <- W[j, i] <- mag * edge_sign
    }
  }
  diag(W) <- 0
  worst_kkt <- max(vapply(fits, function(f) f$kkt, numeric(1)))

  # Nodewise representation for net_predict(): center/scale map raw data to the
  # glmnet predictor scale, beta_std holds the (unthresholded) effective
  # coefficients on that scale, intercept the (effective) intercept.
  intercepts <- vapply(fits, function(f) f$b0_eff, numeric(1))
  families   <- vapply(fits, function(f) f$family, character(1))
  # The p-wide beta matrix is one coefficient per predictor VARIABLE, which can
  # only represent a design where every predictor is a single column. With a
  # multi-level node present it cannot, so it is left empty and net_predict()
  # refuses rather than predicting from a truncated model.
  B_eff <- matrix(0, p, p, dimnames = list(labels, labels))
  if (!has_multi) {
    for (i in seq_len(p)) B_eff[i, fits[[i]]$pred_var] <- fits[[i]]$beta_eff
  }

  .new_psychnet(W, labels, method = "mgm", directed = FALSE,
                n_obs = n, data = mat,
                extra = list(types = stats::setNames(types, labels),
                             levels = stats::setNames(as.integer(levels_v), labels),
                             kkt = worst_kkt, threshold = threshold,
                             native = FALSE,
                             nlambda = nlambda,
                             lambda_min_ratio = lambda_min_ratio,
                             nodewise = list(intercept = intercepts,
                                             beta_std = B_eff,
                                             families = families,
                                             exact = !has_multi,
                                             center = resp_center,
                                             scale = resp_scale,
                                             resp_center = resp_center,
                                             resp_scale = resp_scale)))
}
