# Stationarity (KKT) residual of an L1-penalized GLM fit

Dependency-free correctness certificate for a nodewise lasso, analogous
to
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
for the graphical lasso. With standardized predictors `X` and fitted
mean `mu` (identity link for gaussian, logistic for binomial), the
subgradient conditions are \\n^{-1} X_j^\top (y - \mu) = \lambda\\
\mathrm{sign}(\beta_j)\\ for active coordinates and \\\|n^{-1} X_j^\top
(y - \mu)\| \le \lambda\\ otherwise. Near-zero certifies the
penalized-likelihood optimum.

## Usage

``` r
glm_lasso_kkt(
  X,
  y,
  b0,
  beta,
  lambda,
  family = "gaussian",
  weights = NULL,
  active_tol = 1e-08
)
```

## Arguments

- X:

  Standardized predictor matrix (mean 0, unit variance columns).

- y:

  Response.

- b0:

  Fitted intercept.

- beta:

  Fitted (standardized) coefficients.

- lambda:

  Penalty.

- family:

  `"gaussian"` or `"binomial"`.

- weights:

  Optional observation weights (`NULL` = unweighted).

- active_tol:

  Magnitude above which a coefficient is "active".

## Value

Maximum absolute stationarity violation (scalar). Near-zero certifies
the fit is at the penalized-likelihood optimum.

## Examples

``` r
set.seed(1)
x <- scale(matrix(stats::rnorm(200 * 3), 200, 3))
y <- as.numeric(x %*% c(0.5, 0, -0.3) + stats::rnorm(200))
# x has centred columns, so the least-squares intercept is mean(y) and the
# slopes are qr.solve(x, y): the unpenalized optimum, residual ~0.
glm_lasso_kkt(x, y, b0 = mean(y), beta = qr.solve(x, y), lambda = 0)
#> [1] 1.854072e-16
# The generating coefficients are not the sample optimum.
glm_lasso_kkt(x, y, b0 = mean(y), beta = c(0.5, 0, -0.3), lambda = 0)
#> [1] 0.07437882
```
