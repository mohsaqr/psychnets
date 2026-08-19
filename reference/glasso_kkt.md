# Graphical-lasso stationarity (KKT) residual

A dependency-free correctness certificate for a fitted Gaussian
graphical model. For the convex objective \$\$\min\_{\Theta \succ 0}
-\log\det\Theta + \mathrm{tr}(S\Theta) + \rho \sum\_{i \neq j}
\|\Theta\_{ij}\|\$\$ (off-diagonal penalty), let \\W = \Theta^{-1}\\.
The subgradient optimality conditions are \\W\_{ii} = S\_{ii}\\;
\\W\_{ij} - S\_{ij} = \rho\\\mathrm{sign}(\Theta\_{ij})\\ where
\\\Theta\_{ij} \neq 0\\; and \\\|W\_{ij} - S\_{ij}\| \le \rho\\
otherwise. By strict convexity, a precision matrix with zero violation
is the unique global optimum, so a near-zero return certifies
correctness independently of any reference solver.

## Usage

``` r
glasso_kkt(
  theta,
  cor_matrix,
  rho,
  active_tol = 1e-08,
  penalize_diagonal = FALSE
)
```

## Arguments

- theta:

  Precision matrix to test.

- cor_matrix:

  Correlation / covariance the model was fit to.

- rho:

  Scalar penalty.

- active_tol:

  Magnitude above which an off-diagonal entry is "active".

- penalize_diagonal:

  Was the diagonal penalized in the fit being graded? Default `FALSE`,
  matching
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)'s
  default.

## Value

Maximum absolute stationarity violation (scalar); 0 = exact optimum.

## Details

When the penalty also applies to the diagonal
(`penalize_diagonal = TRUE`, the objective adds \\\rho \sum_i
\|\Theta\_{ii}\|\\), the diagonal condition becomes \\W\_{ii} =
S\_{ii} + \rho\\. Grading a diagonally-penalized fit with the default
would report a spurious violation of exactly \\\rho\\, so this flag must
match the fit being tested.

## Examples

``` r
S <- 0.5^abs(outer(1:5, 1:5, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 200)
glasso_kkt(fit$precision, S, fit$lambda)
#> [1] 7.168638e-11
```
