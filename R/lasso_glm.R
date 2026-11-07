# Clean-room L1-penalized GLM (Friedman, Hastie & Tibshirani 2010, JSS) -- the
# nodewise kernel shared by ising_fit() (logistic link) and mgm_fit() (gaussian
# + logistic). Pure base R. Penalized IRLS: an outer Newton loop forms a
# quadratic (weighted-least-squares) approximation of the log-likelihood, and an
# inner coordinate descent solves the resulting weighted lasso. For the gaussian
# link the working weights are constant, so there is no outer loop: the weighted
# lasso is solved exactly by `.gauss_lasso()` (Gram coordinate descent plus an
# exact active-set KKT solve).
#
# Predictors are standardized to mean 0, unit (population) variance before
# fitting and coefficients are returned on both scales; the penalty therefore
# acts on standardized coordinates, matching glmnet's default.

# Validate optional observation weights against the number of rows used, and
# rescale them to a maximum of 1 -- mgm::mgm's convention (`weights /
# max(weights)`, with the effective sample size sum(weights)). The EBIC
# penalty uses log(sum(weights)) while the deviance scales with the weights, so
# without a fixed scale `weights = rep(1/n, n)` selects the empty network and
# `rep(2, n)` a denser one than `rep(1, n)`. After rescaling the network depends
# only on the relative weights.
.check_weights <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  if (!is.numeric(weights) && !is.logical(weights)) {
    .psn_stop("`weights` must be a numeric vector.")
  }
  weights <- as.numeric(weights)
  if (length(weights) != n) {
    .psn_stop(sprintf("`weights` must have one value per row (%d); got %d.",
                 n, length(weights)), call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) {
    .psn_stop("`weights` must be finite, non-negative, non-missing, and not all zero.",
         call. = FALSE)
  }
  weights / max(weights)
}

# Weighted lasso coordinate descent of z on standardized X with weights w.
# Intercept is unpenalized. Returns updated (b0, beta).
.wls_lasso <- function(X, z, w, lambda, b0, beta, max_inner, tol_inner) {
  n <- nrow(X); p <- ncol(X)
  sw <- sum(w)
  xx <- colSums(w * X^2) / n            # per-column weighted second moment
  resid <- z - b0 - as.numeric(X %*% beta)
  for (inner in seq_len(max_inner)) {
    max_diff <- 0
    # intercept (unpenalized)
    db0 <- sum(w * resid) / sw
    b0 <- b0 + db0
    resid <- resid - db0
    for (j in seq_len(p)) {
      if (xx[j] < 1e-12) next
      rho_j <- sum(w * X[, j] * resid) / n + xx[j] * beta[j]
      new_j <- .soft(rho_j, lambda) / xx[j]
      d <- new_j - beta[j]
      if (d != 0) {
        resid <- resid - X[, j] * d
        if (abs(d) > max_diff) max_diff <- abs(d)
        beta[j] <- new_j
      }
    }
    if (max_diff < tol_inner) break
  }
  list(b0 = b0, beta = beta)
}

# Exact weighted lasso for a gaussian response (the mgm gaussian nodes):
# minimise (2n)^-1 sum_i w_i (z_i - b0 - x_i'beta)^2 + lambda * ||beta||_1.
#
# The quadratic is fixed, so the solve runs on the Gram matrix of the design
# with an intercept column (G = Xa' W Xa / n, cvec = Xa' W z / n): one
# coordinate update costs O(p), not O(n p). Plain coordinate descent with a
# coefficient-change stop is NOT enough on collinear items -- its contraction
# rate approaches 1, so each sweep moves the coefficients by less than the stop
# rule while the iterate is still far from the optimum (KKT ~1e-3 at 100
# sweeps). The descent is therefore only used to find the support and signs:
# every time it settles (and every 25 sweeps) the KKT system on that active set
# is solved exactly, G_AA theta_A = cvec_A - lambda s_A, and accepted when the
# solution keeps the signs and every inactive gradient is within the penalty.
# That is the optimum itself, not an approximation of it. If no exact solve is
# accepted the descent runs to its own stopping rule, and `converged` reports
# whether the stationarity residual reached `kkt_tol`.
#
# `gram` (list(G, cvec)) may be supplied by a caller that fits a whole lambda
# path on the same X, z, w.
#' @noRd
.gauss_gram <- function(X, z, w) {
  n <- nrow(X)
  Xa <- cbind(1, X)
  list(G = crossprod(Xa, w * Xa) / n,
       cvec = as.numeric(crossprod(Xa, w * z)) / n)
}

#' @noRd
.gauss_active_solve <- function(G, cvec, theta, lambda, usable) {
  k <- length(theta)
  active <- usable & (seq_len(k) == 1L | theta != 0)
  A <- which(active)
  s <- sign(theta[A]); s[1L] <- 0                       # intercept unpenalized
  qa <- qr(G[A, A, drop = FALSE])
  if (qa$rank < length(A)) return(NULL)                 # collinear active set
  th <- qr.coef(qa, cvec[A] - lambda * s)
  if (anyNA(th) || any(sign(th[-1L]) != s[-1L])) return(NULL)
  out <- numeric(k); out[A] <- th
  inactive <- which(!active & usable)
  if (length(inactive)) {
    g <- cvec[inactive] - as.numeric(G[inactive, A, drop = FALSE] %*% th)
    if (any(abs(g) > lambda * (1 + 1e-10) + 1e-14)) return(NULL)
  }
  out
}

#' @noRd
.gauss_lasso <- function(X, z, w, lambda, b0, beta, gram = NULL,
                         max_iter = 10000L, tol_cd = 1e-12, kkt_tol = 1e-9) {
  p <- ncol(X)
  gram <- gram %||% .gauss_gram(X, z, w)
  G <- gram$G; cvec <- gram$cvec
  dg <- diag(G)
  usable <- dg > 1e-12
  pen <- c(0, rep(lambda, p))
  theta <- c(b0, beta)
  theta[!usable] <- 0
  grad <- cvec - as.numeric(G %*% theta)
  kkt_of <- function(theta, grad) {
    act <- usable & theta != 0
    act[1L] <- FALSE
    ina <- usable & theta == 0
    ina[1L] <- FALSE
    max(abs(grad[1L]),
        if (any(act)) max(abs(grad[act] - lambda * sign(theta[act]))) else 0,
        if (any(ina)) max(pmax(abs(grad[ina]) - lambda, 0)) else 0)
  }
  coords <- which(usable)
  converged <- FALSE
  iter <- 0L
  # Gauss-Seidel coordinate descent is genuinely sequential (each update reads
  # the gradient the previous one changed), hence the loops.
  for (iter in seq_len(max_iter)) {
    max_diff <- 0
    for (j in coords) {
      new_j <- .soft(grad[j] + dg[j] * theta[j], pen[j]) / dg[j]
      d <- new_j - theta[j]
      if (d != 0) {
        grad <- grad - G[, j] * d
        theta[j] <- new_j
        if (abs(d) > max_diff) max_diff <- abs(d)
      }
    }
    if (max_diff < 1e-6 || iter %% 25L == 0L) {
      exact <- .gauss_active_solve(G, cvec, theta, lambda, usable)
      if (!is.null(exact)) {
        theta <- exact
        converged <- TRUE
        break
      }
    }
    if (max_diff < tol_cd) {
      grad <- cvec - as.numeric(G %*% theta)            # drop accumulated rounding
      converged <- kkt_of(theta, grad) <= kkt_tol
      break
    }
  }
  list(b0 = theta[1L], beta = theta[-1L], converged = converged, iter = iter)
}

# Single penalized-GLM fit at one lambda (standardized X). family in
# {"gaussian","binomial"}. Returns standardized coefficients + intercept, and
# `converged`: for the gaussian family whether the exact solve (or the descent's
# stationarity check) succeeded; for the binomial family whether the outer IRLS
# loop met its coefficient-change rule. A loop that has not converged after
# `max_outer` passes is continued for up to 4 * `max_outer` more: a slow but
# finite fit (extreme weights or probabilities) converges there, while under
# (quasi-)complete separation the coefficients keep growing without bound.
# `separation` is set by a drift test on that continuation (see below).
.glm_lasso_fit <- function(X, y, family, lambda,
                           max_outer = 100L, tol_outer = 1e-7,
                           max_inner = 100L, tol_inner = 1e-7,
                           b0 = NULL, beta = NULL, weights = NULL,
                           gram = NULL) {
  n <- nrow(X); p <- ncol(X)
  wobs <- if (is.null(weights)) rep(1, n) else weights
  if (is.null(beta)) beta <- numeric(p)
  if (is.null(b0))   b0 <- if (family == "binomial")
    stats::qlogis(min(max(stats::weighted.mean(y, wobs), 1e-10), 1 - 1e-10))
    else stats::weighted.mean(y, wobs)

  if (family == "gaussian") {
    fit <- .gauss_lasso(X, y, wobs, lambda, b0, beta, gram = gram)
    return(list(b0 = fit$b0, beta = fit$beta, converged = fit$converged,
                separation = FALSE))
  }

  # binomial: penalized IRLS. Observation weights multiply the working weights.
  # Sequential by nature (each pass linearizes at the previous solution), hence
  # the explicit loop; it stops at the first pass meeting the change rule.
  irls <- function(b0, beta, passes, floor = 1e-5) {
    for (outer in seq_len(passes)) {
      eta <- b0 + as.numeric(X %*% beta)
      pr  <- 1 / (1 + exp(-eta))
      pir <- pmax(pr * (1 - pr), floor)         # bounded IRLS weights (glmnet)
      z   <- eta + (y - pr) / pir               # working response
      v   <- wobs * pir                         # weighted least-squares weights
      b0_old <- b0; beta_old <- beta
      fit <- .wls_lasso(X, z, v, lambda, b0, beta, max_inner, tol_inner)
      b0 <- fit$b0; beta <- fit$beta
      if (max(abs(beta - beta_old)) < tol_outer && abs(b0 - b0_old) < tol_outer)
        return(list(b0 = b0, beta = beta, converged = TRUE))
    }
    list(b0 = b0, beta = beta, converged = FALSE)
  }
  run <- irls(b0, beta, max_outer)
  separation <- FALSE
  eta <- run$b0 + as.numeric(X %*% run$beta)
  # The glmnet floor on the working weights is active once |eta| > ~11.5; there
  # the steps shrink and the change rule can stop a fit that is still drifting
  # (separation) or still far from a slow optimum. Continue with (nearly) exact
  # Newton weights, then test the drift: along the direction the fit moved, a
  # finite optimum's (strictly concave) objective falls, while under
  # (quasi-)complete separation it keeps rising toward its supremum.
  if (!run$converged || any(abs(eta[wobs > 0]) > 11)) {
    th0 <- c(run$b0, run$beta)
    run <- irls(run$b0, run$beta, 4L * max_outer, floor = 1e-14)
    th1 <- c(run$b0, run$beta)
    objective <- function(th) {                  # penalized log-likelihood / n
      e <- th[1L] + as.numeric(X %*% th[-1L])
      log1pexp <- ifelse(e > 0, e + log1p(exp(-e)), log1p(exp(e)))
      sum(wobs * (y * e - log1pexp)) / n - lambda * sum(abs(th[-1L]))
    }
    d <- th1 - th0
    step <- sqrt(sum(d^2))
    at <- objective(th1)
    eta <- th1[1L] + as.numeric(X %*% th1[-1L])
    separation <- if (step > 1e-6) {
      objective(th1 + 10 * d / step) >= at - 1e-12 * max(1, abs(at))
    } else {
      any(abs(eta[wobs > 0]) > 30)               # stalled at saturation
    }
    if (separation) run$converged <- FALSE
  }
  list(b0 = run$b0, beta = run$beta, converged = run$converged,
       separation = separation)
}

#' Stationarity (KKT) residual of an L1-penalized GLM fit
#'
#' Dependency-free correctness certificate for a nodewise lasso, analogous to
#' [glasso_kkt()] for the graphical lasso. With standardized predictors `X` and
#' fitted mean `mu` (identity link for gaussian, logistic for binomial), the
#' subgradient conditions are \eqn{n^{-1} X_j^\top (y - \mu) = \lambda\,
#' \mathrm{sign}(\beta_j)} for active coordinates and
#' \eqn{|n^{-1} X_j^\top (y - \mu)| \le \lambda} otherwise. Near-zero certifies
#' the penalized-likelihood optimum.
#'
#' @param X Standardized predictor matrix (mean 0, unit variance columns).
#' @param y Response.
#' @param b0 Fitted intercept.
#' @param beta Fitted (standardized) coefficients.
#' @param lambda Penalty.
#' @param family `"gaussian"` or `"binomial"`.
#' @param weights Optional observation weights (`NULL` = unweighted).
#' @param active_tol Magnitude above which a coefficient is "active".
#' @return Maximum absolute stationarity violation (scalar). Near-zero certifies
#'   the fit is at the penalized-likelihood optimum.
#' @examples
#' set.seed(1)
#' x <- scale(matrix(stats::rnorm(200 * 3), 200, 3))
#' y <- as.numeric(x %*% c(0.5, 0, -0.3) + stats::rnorm(200))
#' # x has centred columns, so the least-squares intercept is mean(y) and the
#' # slopes are qr.solve(x, y): the unpenalized optimum, residual ~0.
#' glm_lasso_kkt(x, y, b0 = mean(y), beta = qr.solve(x, y), lambda = 0)
#' # The generating coefficients are not the sample optimum.
#' glm_lasso_kkt(x, y, b0 = mean(y), beta = c(0.5, 0, -0.3), lambda = 0)
#' @export
glm_lasso_kkt <- function(X, y, b0, beta, lambda, family = "gaussian",
                          weights = NULL, active_tol = 1e-8) {
  family <- .psn_match_arg(family, c("gaussian", "binomial"))
  if (!(is.matrix(X) && is.numeric(X)) || length(y) != nrow(X) ||
      length(beta) != ncol(X) ||
      !(is.null(weights) || length(weights) == nrow(X))) {
    .psn_stop(paste0("glm_lasso_kkt(): `X` must be a numeric matrix with one ",
                     "row per `y` (and per `weights`) and one column per `beta`."))
  }
  if (!(is.numeric(y) && all(is.finite(y))) ||
      (family == "binomial" && any(y < 0 | y > 1)))
    .psn_stop("glm_lasso_kkt(): `y` must be finite (in [0, 1] for binomial).")
  if (!all(is.finite(X)) || !(is.numeric(beta) && all(is.finite(beta))))
    .psn_stop("glm_lasso_kkt(): `X` and `beta` must be finite.")
  .psn_check_number(b0, "b0")
  .psn_check_number(lambda, "lambda", lower = 0)
  .psn_check_number(active_tol, "active_tol", lower = 0)
  if (!is.null(weights) &&
      !(is.numeric(weights) && all(is.finite(weights)) && all(weights >= 0) &&
        sum(weights) > 0))
    .psn_stop(paste0("glm_lasso_kkt(): `weights` must be finite, non-negative ",
                     "and not all zero (a negative weight turns the objective ",
                     "into one with no minimum, so stationarity certifies ",
                     "nothing)."))
  n <- nrow(X)
  wobs <- if (is.null(weights)) rep(1, n) else weights
  # Normalize by n (the row count), matching .wls_lasso's objective scaling, so
  # the certificate and the solver speak the same lambda scale under weighting.
  eta <- b0 + as.numeric(X %*% beta)
  mu  <- if (family == "binomial") 1 / (1 + exp(-eta)) else eta
  grad <- as.numeric(crossprod(X, wobs * (y - mu))) / n
  v_0 <- abs(sum(wobs * (y - mu)) / n)         # unpenalized intercept stationarity
  active <- abs(beta) > active_tol
  v_a <- if (any(active))  max(abs(grad[active] - lambda * sign(beta[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(grad[!active]) - lambda, 0)) else 0
  max(v_0, v_a, v_i)
}

# Standardize columns to mean 0, unit population sd. Returns standardized matrix
# plus centers/scales for back-transformation.
.standardize <- function(X, w = NULL) {
  if (is.null(w)) {
    ctr <- colMeans(X)
    Xc  <- sweep(X, 2L, ctr, "-")
    scl <- sqrt(colMeans(Xc^2))
  } else {
    sw  <- sum(w)
    ctr <- colSums(w * X) / sw
    Xc  <- sweep(X, 2L, ctr, "-")
    scl <- sqrt(colSums(w * Xc^2) / sw)
  }
  scl[scl < 1e-12] <- 1
  list(X = sweep(Xc, 2L, scl, "/"), center = ctr, scale = scl)
}

# Nodewise lasso path with EBIC selection for one response column.
# Returns the EBIC-selected fit: original-scale coefficients, intercept,
# standardized coefficients (for certification), the lambda, and the KKT value.
.nodewise_ebic <- function(X, y, family, gamma, nlambda, lambda_min_ratio,
                           weights = NULL) {
  n <- nrow(X); p <- ncol(X)
  wobs <- if (is.null(weights)) rep(1, n) else weights
  sw <- sum(wobs)
  std <- .standardize(X, weights)
  Xs <- std$X

  # lambda_max: smallest penalty zeroing all coefficients (KKT at beta = 0).
  # The empty-model fitted mean is mu0 = weighted.mean(y) on the response scale;
  # the intercept on the model's linear-predictor scale is the logit of that for
  # a binomial node, not the probability itself. The penalty grid and EBIC use
  # the effective sample size sw = sum(weights).
  mu0 <- stats::weighted.mean(y, wobs)
  lambda_max <- max(abs(crossprod(Xs, wobs * (y - mu0)))) / n
  # A binary response with one value on every positively weighted row (left
  # by `min_sum` filtering or zero weights) has no finite intercept: there is
  # no optimum to certify, so the node is flagged degenerate, not certified.
  yw <- y[wobs > 0]
  if (family == "binomial" && all(yw == yw[1L])) {
    return(list(beta = numeric(p), b0 = stats::qlogis(min(max(mu0, 1e-10), 1 - 1e-10)),
                beta_std = numeric(p), lambda = 0, kkt = NA_real_,
                converged = FALSE, separation = FALSE, degenerate = TRUE))
  }
  if (lambda_max < 1e-12) {
    b0_empty <- if (family == "binomial") {
      stats::qlogis(min(max(mu0, 1e-10), 1 - 1e-10))
    } else mu0
    return(list(beta = numeric(p), b0 = b0_empty, beta_std = numeric(p),
                lambda = 0, kkt = 0, converged = TRUE, separation = FALSE))
  }
  lambdas <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                     length.out = nlambda))

  best <- NULL; best_ebic <- Inf
  b0 <- NULL; beta <- NULL
  # The gaussian quadratic does not change along the path: build its Gram once.
  gram <- if (family == "gaussian") .gauss_gram(Xs, y, wobs) else NULL
  for (lam in lambdas) {
    fit <- .glm_lasso_fit(Xs, y, family, lam, b0 = b0, beta = beta,
                          weights = weights, gram = gram)
    b0 <- fit$b0; beta <- fit$beta                      # warm start
    eta <- fit$b0 + as.numeric(Xs %*% fit$beta)
    if (family == "binomial") {
      pr <- 1 / (1 + exp(-eta))
      pr <- pmin(pmax(pr, 1e-10), 1 - 1e-10)
      dev <- -2 * sum(wobs * (y * log(pr) + (1 - y) * log(1 - pr)))
    } else {
      dev <- sum(wobs * (y - eta)^2)                    # gaussian deviance = RSS (glmnet/mgm)
    }
    df <- sum(abs(fit$beta) > 1e-10)
    ebic <- dev + df * log(sw) + 2 * gamma * df * log(p)
    if (ebic < best_ebic) {
      best_ebic <- ebic
      best <- list(b0 = fit$b0, beta_std = fit$beta, lambda = lam,
                   converged = fit$converged, separation = fit$separation)
    }
  }
  # Inactive coordinates must be EXACT zeros. Coordinate descent can leave
  # BLAS-rounding crumbs (~1e-17) on coordinates whose gradient sits exactly at
  # the penalty boundary, and whether the crumb is 0.0 or 4e-17 differs between
  # BLAS implementations. The EBIC df above already treats |beta| <= 1e-10 as
  # "not a parameter"; the returned support must agree, or the AND rule counts
  # a platform-dependent crumb as an edge.
  best$beta_std[abs(best$beta_std) <= 1e-10] <- 0
  beta_orig <- best$beta_std / std$scale
  kkt <- glm_lasso_kkt(Xs, y, best$b0, best$beta_std, best$lambda, family,
                       weights = weights)
  list(beta = beta_orig, b0 = best$b0, beta_std = best$beta_std,
       lambda = best$lambda, kkt = kkt, ebic = best_ebic,
       converged = best$converged && kkt <= .nodewise_kkt_tol,
       separation = best$separation)
}

# Stationarity residual above which a nodewise fit is reported as unconverged:
# the default tolerance of certificate(), so a fit the solver accepts is never
# one the certificate then rejects without a warning.
.nodewise_kkt_tol <- 1e-6

# Surface nodewise fits that did not reach their optimum: one classed warning
# per estimator call, naming the nodes. `fits` are .nodewise_ebic() (or
# per-node) results carrying `converged`/`kkt`.
#' @noRd
.warn_nodewise_convergence <- function(fits, labels, caller) {
  degenerate <- vapply(fits, function(f) isTRUE(f$degenerate), logical(1))
  if (any(degenerate)) {
    .psn_warn(sprintf(paste0(
      "%s(): node(s) %s take a single value on every retained (positively ",
      "weighted) row, e.g. after `min_sum` filtering, so their logistic ",
      "regression has no finite intercept; the fit is not certified."),
      caller, paste(labels[degenerate], collapse = ", ")),
      class = "psychnet_degenerate_sample")
  }
  bad <- !vapply(fits, function(f) isTRUE(f$converged), logical(1)) & !degenerate
  if (!any(bad)) return(invisible(any(degenerate)))
  kkt <- vapply(fits, function(f) f$kkt %||% NA_real_, numeric(1))
  .psn_warn(sprintf(paste0(
    "%s(): the nodewise fit(s) for %s did not reach the optimum (worst ",
    "stationarity residual %.3g, tolerance %.0e); the returned weights are ",
    "not certified."), caller, paste(labels[bad], collapse = ", "),
    max(kkt[bad]), .nodewise_kkt_tol), class = "psychnet_no_converge")
  invisible(TRUE)
}
