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
                           b0 = NULL, beta = NULL) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(beta)) beta <- numeric(p)
  if (is.null(b0))   b0 <- if (family == "binomial") log(mean(y) / (1 - mean(y))) else mean(y)

  if (family == "gaussian") {
    w <- rep(1, n)
    fit <- .wls_lasso(X, y, w, lambda, b0, beta, max_inner, tol_inner)
    return(list(b0 = fit$b0, beta = fit$beta))
  }

  # binomial: penalized IRLS
  for (outer in seq_len(max_outer)) {
    eta <- b0 + as.numeric(X %*% beta)
    pr  <- 1 / (1 + exp(-eta))
    w   <- pmax(pr * (1 - pr), 1e-5)            # bounded working weights (glmnet)
    z   <- eta + (y - pr) / w                   # working response
    b0_old <- b0; beta_old <- beta
    fit <- .wls_lasso(X, z, w, lambda, b0, beta, max_inner, tol_inner)
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
#' @param active_tol Magnitude above which a coefficient is "active".
#' @return Maximum absolute stationarity violation (scalar).
#' @export
glm_lasso_kkt <- function(X, y, b0, beta, lambda, family = "gaussian",
                          active_tol = 1e-8) {
  n <- nrow(X)
  eta <- b0 + as.numeric(X %*% beta)
  mu  <- if (family == "binomial") 1 / (1 + exp(-eta)) else eta
  grad <- as.numeric(crossprod(X, y - mu)) / n
  v_0 <- abs(mean(y - mu))                     # unpenalized intercept stationarity
  active <- abs(beta) > active_tol
  v_a <- if (any(active))  max(abs(grad[active] - lambda * sign(beta[active]))) else 0
  v_i <- if (any(!active)) max(pmax(abs(grad[!active]) - lambda, 0)) else 0
  max(v_0, v_a, v_i)
}

# Standardize columns to mean 0, unit population sd. Returns standardized matrix
# plus centers/scales for back-transformation.
.standardize <- function(X) {
  ctr <- colMeans(X)
  Xc  <- sweep(X, 2L, ctr, "-")
  scl <- sqrt(colMeans(Xc^2))
  scl[scl < 1e-12] <- 1
  list(X = sweep(Xc, 2L, scl, "/"), center = ctr, scale = scl)
}

# Nodewise lasso path with EBIC selection for one response column.
# Returns the EBIC-selected fit: original-scale coefficients, intercept,
# standardized coefficients (for certification), the lambda, and the KKT value.
.nodewise_ebic <- function(X, y, family, gamma, nlambda, lambda_min_ratio) {
  n <- nrow(X); p <- ncol(X)
  std <- .standardize(X)
  Xs <- std$X

  # lambda_max: smallest penalty zeroing all coefficients (KKT at beta = 0).
  # The empty-model fitted mean is mu0 = mean(y) on the response scale; the
  # intercept on the model's linear-predictor scale is the logit of that for a
  # binomial node, not the probability itself.
  mu0 <- mean(y)
  lambda_max <- max(abs(crossprod(Xs, y - mu0))) / n
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
    fit <- .glm_lasso_fit(Xs, y, family, lam, b0 = b0, beta = beta)
    b0 <- fit$b0; beta <- fit$beta                      # warm start
    eta <- fit$b0 + as.numeric(Xs %*% fit$beta)
    if (family == "binomial") {
      pr <- 1 / (1 + exp(-eta))
      pr <- pmin(pmax(pr, 1e-10), 1 - 1e-10)
      dev <- -2 * sum(y * log(pr) + (1 - y) * log(1 - pr))
    } else {
      dev <- sum((y - eta)^2)                           # gaussian deviance = RSS (glmnet/mgm)
    }
    df <- sum(abs(fit$beta) > 1e-10)
    ebic <- dev + df * log(n) + 2 * gamma * df * log(p)
    if (ebic < best_ebic) {
      best_ebic <- ebic
      best <- list(b0 = fit$b0, beta_std = fit$beta, lambda = lam)
    }
  }
  beta_orig <- best$beta_std / std$scale
  kkt <- glm_lasso_kkt(Xs, y, best$b0, best$beta_std, best$lambda, family)
  list(beta = beta_orig, b0 = best$b0, beta_std = best$beta_std,
       lambda = best$lambda, kkt = kkt, ebic = best_ebic)
}
