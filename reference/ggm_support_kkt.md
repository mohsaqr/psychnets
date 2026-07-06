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

  Logical p x p matrix; `TRUE` where an edge is allowed.

- active_tol:

  Magnitude above which an off-support entry counts as a nonzero
  violation.

## Value

Maximum absolute stationarity violation (scalar); 0 = exact optimum.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
fit <- ggm_modselect(cor_matrix = S, n = 250)
ggm_support_kkt(fit$precision, S, fit$support)
#> [1] 2.220446e-16
```
