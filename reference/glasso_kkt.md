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

Maximum absolute stationarity violation (scalar); 0 = exact optimum,
`Inf` for an asymmetric or non-positive-definite `theta`.

## Details

When the penalty also applies to the diagonal
(`penalize_diagonal = TRUE`, the objective adds \\\rho \sum_i
\|\Theta\_{ii}\|\\), the diagonal condition becomes \\W\_{ii} =
S\_{ii} + \rho\\. Grading a diagonally-penalized fit with the default
would report a spurious violation of exactly \\\rho\\, so this flag must
match the fit being tested.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
which reads the stored residual of a fitted network.

## Examples

``` r
S <- 0.5^abs(outer(1:5, 1:5, "-"))
# At rho = 0 the optimum is the inverse of S, so its residual is ~0 ...
glasso_kkt(solve(S), S, rho = 0)
#> [1] 1.110223e-16
# ... and a perturbed precision matrix violates the conditions.
glasso_kkt(solve(S) + diag(0.1, 5), S, rho = 0)
#> [1] 0.135257
# A fitted network reports the same residual through certificate().
certificate(ebic_glasso(cor_matrix = S, n = 200))
#>   method  certificate kind certified
#> 1 glasso 7.168638e-11  kkt      TRUE
```
