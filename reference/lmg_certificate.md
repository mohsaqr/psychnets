# Relative-importance (LMG / Shapley) certificate

By Shapley efficiency, the importance shares a node receives from its
predictors sum exactly to that node's full-model R-squared. Returns the
maximum absolute deviation from that identity; near zero certifies the
decomposition.

## Usage

``` r
lmg_certificate(x)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object produced by
  [`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md).

## Value

Maximum absolute deviation of incoming-share sums from the per-node
R-squared (scalar); 0 = exact decomposition.

## Examples

``` r
S <- 0.4^abs(outer(1:5, 1:5, "-"))
lmg_certificate(relimp_network(cor_matrix = S))
#> [1] 5.551115e-17
```
