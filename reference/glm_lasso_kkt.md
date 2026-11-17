# Check that coefficients solve a lasso-penalized regression

`glm_lasso_kkt()` checks whether coefficients solve a lasso-penalized
linear or logistic regression, given its predictors, response,
intercept, coefficients and penalty. Such a regression, with a penalty
on the absolute coefficients (the lasso), is a convex problem, and its
coefficients are the optimum if and only if they satisfy the optimality
(Karush-Kuhn-Tucker) conditions.

Ising models and mixed graphical models are estimated node by node: each
variable is regressed on all the others with the lasso, logistic for a
binary variable and linear for a continuous one, and the coefficients
become the edges of the network. Checking the conditions therefore shows
whether a regression was solved, without fitting it a second time. In
words, the residuals must average to zero; for every coefficient that is
not zero, the average product of its predictor with the residuals must
equal the penalty, with the sign of the coefficient; and for every
coefficient that is zero it must not exceed the penalty in absolute
value. With `lambda = 0` the check applies to an unpenalized
maximum-likelihood regression.

`glm_lasso_kkt()` returns the largest violation of the conditions for
one regression; a value near zero shows that the coefficients are the
optimum, up to numerical precision.

[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
store the largest value over their nodewise regressions, and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports it. Graphical lasso networks are checked by
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md).

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

  Numeric predictor matrix, one row per observation, on the scale on
  which the fit was penalized.

- y:

  Numeric response, one value per row of `X`; between 0 and 1 for
  `family = "binomial"`.

- b0:

  Fitted intercept, a single finite number.

- beta:

  Fitted coefficients, one per column of `X`.

- lambda:

  The penalty: a single finite number \>= 0.

- family:

  `"gaussian"` (default, linear regression) or `"binomial"` (logistic
  regression).

- weights:

  Optional observation weights, finite, non-negative and not all zero,
  one per row of `X`. Default `NULL`, equal weights.

- active_tol:

  Absolute value above which a coefficient counts as nonzero. Default
  `1e-8`.

## Value

A single number, the largest absolute violation of the optimality
conditions; 0 at the exact optimum.

## Details

With fitted mean \\\mu = \beta_0 + X\beta\\ (gaussian) or \\\mu = 1 /
(1 + e^{-(\beta_0 + X\beta)})\\ (binomial), observation weights \\w\\
and \\n\\ rows, let \\g_j = n^{-1} \sum_i w_i x\_{ij} (y_i - \mu_i)\\.
The conditions are \\n^{-1} \sum_i w_i (y_i - \mu_i) = 0\\ for the
unpenalized intercept, \\g_j = \lambda\\\mathrm{sign}(\beta_j)\\ for
every active coefficient, and \\\|g_j\| \le \lambda\\ otherwise. The
function returns the largest absolute deviation from these conditions.
The penalty acts on the scale of `X`; the package's nodewise fits
standardize the predictors to mean 0 and unit variance, so a fit is
checked with its standardized predictors and coefficients.

## Conditions

`psychnet_bad_input` when `X` is not a finite numeric matrix matching
`y`, `beta` and `weights` in size, `y` is not finite (or not between 0
and 1 for binomial), `b0`, `lambda` or `active_tol` is not a valid
single number, `family` is not one of the two choices, or `weights`
contain a negative or non-finite value or are all zero.

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
which reports the stored value of a fitted network;
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md);
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
and
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
for the other optimality checks.

## Examples

``` r
set.seed(1)
x <- scale(matrix(stats::rnorm(200 * 3), 200, 3))
y <- as.numeric(x %*% c(0.5, 0, -0.3) + stats::rnorm(200))
# x has centred columns, so the least-squares intercept is mean(y) and the
# slopes are qr.solve(x, y): the unpenalized optimum, violation ~0.
glm_lasso_kkt(x, y, b0 = mean(y), beta = qr.solve(x, y), lambda = 0)
#> [1] 1.854072e-16
# The generating coefficients are not the sample optimum.
glm_lasso_kkt(x, y, b0 = mean(y), beta = c(0.5, 0, -0.3), lambda = 0)
#> [1] 0.07437882
```
