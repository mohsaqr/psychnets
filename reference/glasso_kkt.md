# Check that a precision matrix solves the graphical lasso

`glasso_kkt()` checks whether a precision matrix solves the graphical
lasso for a given correlation or covariance matrix and penalty. The
graphical lasso estimates a Gaussian graphical model as the precision
(inverse covariance) matrix that maximizes the Gaussian likelihood minus
a penalty on the absolute partial associations (Friedman et al., 2008),
and a matrix is its solution if and only if it satisfies the optimality
(Karush-Kuhn-Tucker) conditions of the problem.

The objective is strictly convex, so it has exactly one optimum, and
checking the conditions shows whether a network is the graphical lasso
solution without solving the problem a second time. In words, the
covariance implied by the precision matrix must equal the observed
covariance on the diagonal, differ from it by exactly the penalty for
every pair joined by an edge, and by at most the penalty for every pair
that is not.

`glasso_kkt()` returns the largest violation of the conditions; a value
near zero shows that the matrix is the optimum, up to numerical
precision.

Networks from
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)
store this value, and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports it. Unpenalized networks restricted to a graph are checked by
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md),
and nodewise Ising and mixed graphical models by
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md).

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

  The precision matrix to check: a finite, square numeric matrix.

- cor_matrix:

  The correlation or covariance matrix the model was fitted to, of the
  same size as `theta`.

- rho:

  The penalty: a single finite number \>= 0.

- active_tol:

  Absolute value above which an off-diagonal entry of `theta` counts as
  an edge. Default `1e-8`.

- penalize_diagonal:

  Whether the fit being checked also penalized the diagonal. Default
  `FALSE`, the default of
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

## Value

A single number, the largest absolute violation of the optimality
conditions: 0 at the exact optimum, `Inf` for an asymmetric or
non-positive-definite `theta`.

## Details

The objective, with the penalty on the off-diagonal entries, is
\$\$\min\_{\Theta \succ 0} -\log\det\Theta + \mathrm{tr}(S\Theta) + \rho
\sum\_{i \neq j} \|\Theta\_{ij}\|.\$\$ With \\W = \Theta^{-1}\\, the
optimality conditions are \\W\_{ii} = S\_{ii}\\; \\W\_{ij} - S\_{ij} =
\rho\\\mathrm{sign}(\Theta\_{ij})\\ where \\\Theta\_{ij} \neq 0\\; and
\\\|W\_{ij} - S\_{ij}\| \le \rho\\ otherwise. The function returns the
largest absolute deviation from these conditions. Because the objective
is strictly convex, a precision matrix with zero violation is its unique
global optimum.

When the penalty also applies to the diagonal
(`penalize_diagonal = TRUE`, the objective adds \\\rho \sum_i
\|\Theta\_{ii}\|\\), the diagonal condition becomes \\W\_{ii} =
S\_{ii} + \rho\\. Checking a diagonally penalized fit with the default
would report a violation of exactly \\\rho\\, so this argument must
match the fit being checked.

The objective is defined only for symmetric positive-definite matrices,
so a `theta` that is asymmetric (by more than `1e-8`) or not positive
definite is not a solution and returns `Inf`.

## Conditions

`psychnet_bad_input` when `theta` is not a finite square numeric matrix,
`cor_matrix` is not a finite matrix of the same size, `rho` or
`active_tol` is not a single finite number \>= 0, or `penalize_diagonal`
is not `TRUE` or `FALSE`.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
which reports the stored value of a fitted network;
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md);
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
and
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
for the other optimality checks.

## Examples

``` r
S <- 0.5^abs(outer(1:5, 1:5, "-"))
# At rho = 0 the optimum is the inverse of S, so the violation is ~0 ...
glasso_kkt(solve(S), S, rho = 0)
#> [1] 1.110223e-16
# ... and a perturbed precision matrix violates the conditions.
glasso_kkt(solve(S) + diag(0.1, 5), S, rho = 0)
#> [1] 0.135257
# A fitted network reports the same value through certificate().
certificate(ebic_glasso(cor_matrix = S, n = 200))
#>   method  certificate kind certified
#> 1 glasso 7.168638e-11  kkt      TRUE
```
