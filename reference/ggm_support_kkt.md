# Check that a precision matrix is the maximum-likelihood network on a graph

`ggm_support_kkt()` checks whether a precision matrix is the
maximum-likelihood Gaussian graphical model on a given graph, for the
correlation or covariance matrix it was fitted to. Some Gaussian
graphical models are estimated without a penalty: a graph is chosen
first, and the precision matrix is then the maximum-likelihood estimate
under the constraint that every pair of variables not joined by an edge
has a partial correlation of zero.

This constrained estimate is unique and is characterized by two
conditions: the covariance implied by the model reproduces the observed
covariance for every pair joined by an edge and for every variance, and
the precision matrix is zero for every pair that is not joined. A matrix
meeting both conditions is the maximum-likelihood network on that graph.

`ggm_support_kkt()` returns the largest violation of these conditions; a
value near zero shows that the matrix is the constrained optimum, up to
numerical precision.

Networks from
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md),
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)
and `ebic_glasso(refit = "unregularized")` store this value, and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports it. Penalized graphical lasso networks are checked by
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md).

## Usage

``` r
ggm_support_kkt(theta, cor_matrix, support, active_tol = 1e-08)
```

## Arguments

- theta:

  The precision matrix to check: a finite, square numeric matrix.

- cor_matrix:

  The correlation or covariance matrix the model was fitted to, of the
  same size as `theta`.

- support:

  A logical matrix of the same size, `TRUE` where an edge is allowed,
  with no missing values. It is made symmetric, so an upper- or
  lower-triangular matrix describes the same graph; the diagonal is
  ignored.

- active_tol:

  A single finite number \>= 0 (default `1e-8`). It is checked but does
  not affect the result: every entry of `theta` outside the graph counts
  as a violation of its absolute size.

## Value

A single number, the largest absolute violation of the conditions: 0 at
the exact constrained optimum, `Inf` for an asymmetric or
non-positive-definite `theta`.

## Details

For a Gaussian graphical model on a graph \\G\\, with \\W =
\Theta^{-1}\\, the maximum-likelihood conditions are \\W\_{ij} =
S\_{ij}\\ for every \\(i,j)\\ in \\G\\ and on the diagonal, and
\\\Theta\_{ij} = 0\\ for every \\(i,j)\\ not in \\G\\. The function
returns the largest absolute deviation from these conditions; a `theta`
that is asymmetric (by more than `1e-8`) or not positive definite
returns `Inf`.

## Conditions

`psychnet_bad_input` when `theta` is not a finite square numeric matrix,
`cor_matrix` is not a finite matrix of the same size, `support` is not a
logical matrix of that size without missing values, or `active_tol` is
not a single finite number \>= 0.

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
which reports the stored value of a fitted network;
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
and
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md);
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
for penalized networks.

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
# A fitted network reports its value through certificate().
certificate(ggm_modselect(cor_matrix = S, n = 250))
#>   method  certificate kind certified
#> 1    ggm 2.220446e-16  kkt      TRUE
```
