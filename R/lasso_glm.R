# Clean-room L1-penalized GLM (Friedman, Hastie & Tibshirani 2010, JSS) -- the
# nodewise kernel shared by ising_fit() (logistic link) and mgm_fit() (gaussian
# + logistic). Pure base R. Penalized IRLS: an outer Newton loop forms a
# quadratic (weighted-least-squares) approximation of the log-likelihood, and an
# inner coordinate descent solves the resulting weighted lasso. For the gaussian
# link the working weights are constant, so the outer loop converges in one pass.
#
# Predictors are standardized to mean 0, unit (population) variance before
# fitting and coefficients are returned on both scales; the penalty therefore
# acts on standardized coordinates, matching glmnet's default.

# Validate optional observation weights against the number of rows used.
.check_weights <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  weights <- as.numeric(weights)
  if (length(weights) != n) {
    stop(sprintf("`weights` must have one value per row (%d); got %d.",
                 n, length(weights)), call. = FALSE)
  }
  if (anyNA(weights) || any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be non-negative, non-missing, and not all zero.",
         call. = FALSE)
  }
  weights
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

# Single penalized-GLM fit at one lambda (standardized X). family in
# {"gaussian","binomial"}. Returns standardized coefficients + intercept.
.glm_lasso_fit <- function(X, y, family, lambda,
                           max_outer = 100L, tol_outer = 1e-7,
                           max_inner = 100L, tol_inner = 1e-7,
                           b0 = NULL, beta = NULL, weights = NULL) {
  n <- nrow(X); p <- ncol(X)
  wobs <- if (is.null(weights)) rep(1, n) else weights
  if (is.null(beta)) beta <- numeric(p)
  if (is.null(b0))   b0 <- if (family == "binomial")
    stats::qlogis(min(max(stats::weighted.mean(y, wobs), 1e-10), 1 - 1e-10))
    else stats::weighted.mean(y, wobs)

  if (family == "gaussian") {
    fit <- .wls_lasso(X, y, wobs, lambda, b0, beta, max_inner, tol_inner)
    return(list(b0 = fit$b0, beta = fit$beta))
  }

  # binomial: penalized IRLS. Observation weights multiply the working weights.
  for (outer in seq_len(max_outer)) {
    eta <- b0 + as.numeric(X %*% beta)
    pr  <- 1 / (1 + exp(-eta))
    pir <- pmax(pr * (1 - pr), 1e-5)            # bounded IRLS weights (glmnet)
    z   <- eta + (y - pr) / pir                 # working response
    v   <- wobs * pir                           # weighted least-squares weights
    b0_old <- b0; beta_old <- beta
    fit <- .wls_lasso(X, z, v, lambda, b0, beta, max_inner, tol_inner)
    b0 <- fit$b0; beta <- fit$beta
    if (max(abs(beta - beta_old)) < tol_outer && abs(b0 - b0_old) < tol_outer) break
  }
  list(b0 = b0, beta = beta)
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
#' fit <- stats::lm.fit(cbind(1, x), y)
#' glm_lasso_kkt(x, y, fit$coefficients[1], fit$coefficients[-1], lambda = 0)
#' @export
glm_lasso_kkt <- function(X, y, b0, beta, lambda, family = "gaussian",
                          weights = NULL, active_tol = 1e-8) {
  family <- match.arg(family, c("gaussian", "binomial"))
  stopifnot(is.matrix(X), length(y) == nrow(X), length(beta) == ncol(X),
            is.null(weights) || length(weights) == nrow(X))
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
  if (lambda_max < 1e-12) {
    b0_empty <- if (family == "binomial") {
      stats::qlogis(min(max(mu0, 1e-10), 1 - 1e-10))
    } else mu0
    return(list(beta = numeric(p), b0 = b0_empty, beta_std = numeric(p),
                lambda = 0, kkt = 0))
  }
  lambdas <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
                     length.out = nlambda))

  best <- NULL; best_ebic <- Inf
  b0 <- NULL; beta <- NULL
  for (lam in lambdas) {
    fit <- .glm_lasso_fit(Xs, y, family, lam, b0 = b0, beta = beta,
                          weights = weights)
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
      best <- list(b0 = fit$b0, beta_std = fit$beta, lambda = lam)
    }
  }
  beta_orig <- best$beta_std / std$scale
  kkt <- glm_lasso_kkt(Xs, y, best$b0, best$beta_std, best$lambda, family,
                       weights = weights)
  list(beta = beta_orig, b0 = best$b0, beta_std = best$beta_std,
       lambda = best$lambda, kkt = kkt, ebic = best_ebic)
}
