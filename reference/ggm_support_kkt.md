# Constrained Gaussian-MRF (graph-restricted MLE) stationarity residual

Certificate for an *unregularized* Gaussian graphical model whose
precision is constrained to a fixed graph (the estimator behind
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
and
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)).
The maximum-likelihood / maximum-entropy conditions for a Gaussian
Markov random field on a graph \\G\\ are exact: \\W\_{ij} = S\_{ij}\\
for every \\(i,j)\\ on the graph and on the diagonal (\\W =
\Theta^{-1}\\), and \\\Theta\_{ij} = 0\\ for every \\(i,j)\\ not on the
graph. A near-zero return certifies the constrained optimum with no
reference solver.

## Usage

``` r
ggm_support_kkt(theta, cor_matrix, support, active_tol = 1e-08)
```

## Arguments

- theta:

  Precision matrix to test.

- cor_matrix:

  Correlation / covariance the model was fit to.

- support:

  Logical p x p matrix; `TRUE` where an edge is allowed. It is
  symmetrized (`support | t(support)`), so an upper- or lower-triangular
  support describes the same undirected graph; the diagonal is ignored.

- active_tol:

  Magnitude above which an off-support entry counts as a nonzero
  violation.

## Value

Maximum absolute stationarity violation (scalar); 0 = exact optimum,
`Inf` for an asymmetric or non-positive-definite `theta`. A `support`
that is not a logical matrix of the right size, or contains `NA`, is an
error of class `psychnet_bad_input`.

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
which reads the stored residual of a fitted network.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
# An AR(1) correlation has a tridiagonal inverse: the chain graph is exact.
chain <- abs(row(S) - col(S)) <= 1
ggm_support_kkt(solve(S), S, support = chain)
#> [1] 2.220446e-16
# The same precision violates an empty-graph (diagonal) constraint.
ggm_support_kkt(solve(S), S, support = diag(6) > 0)
#> [1] 0.4761905
# A fitted network reports its residual through certificate().
certificate(ggm_modselect(cor_matrix = S, n = 250))
#>   method  certificate kind certified
#> 1    ggm 2.220446e-16  kkt      TRUE
```
