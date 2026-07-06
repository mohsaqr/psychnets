# Bridge centrality

Computes bridge centrality for an undirected weighted network: how
strongly each node connects to communities other than its own. You
supply the community membership; psychnets does not detect it.

## Usage

``` r
net_bridge(x, communities, normalize = FALSE, labels = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix.

- communities:

  Community membership, one entry per node: a vector aligned to the node
  order, or a named vector / list keyed by node label.

- normalize:

  If `TRUE`, divide each metric by the number of available
  other-community nodes (comparable across differently sized networks).
  Default `FALSE`.

- labels:

  Optional node labels (used when `x` is a bare matrix).

## Value

A tidy `data.frame` (class `psychnet_bridge`), one row per node, with
columns `node`, `community`, `bridge_strength`, `bridge_betweenness`,
`bridge_closeness`, `bridge_ei1`, `bridge_ei2`. Visualise with
[`plot.psychnet_bridge()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_bridge.md).

## References

Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality.
*Multivariate Behavioral Research*, 56(2), 353-367.

## Examples

``` r
S <- 0.3^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 400)
net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2))
#>   node community bridge_strength bridge_betweenness bridge_closeness bridge_ei1
#> 1   V1         1       0.0000000                  0       0.06918582  0.0000000
#> 2   V2         1       0.0000000                  3       0.09139545  0.0000000
#> 3   V3         1       0.2729255                  6       0.13741057  0.2729255
#> 4   V4         2       0.2729255                  6       0.13741057  0.2729255
#> 5   V5         2       0.0000000                  3       0.09139545  0.0000000
#> 6   V6         2       0.0000000                  0       0.06918582  0.0000000
#>   bridge_ei2
#> 1 0.00000000
#> 2 0.07448834
#> 3 0.34741385
#> 4 0.34741385
#> 5 0.07448834
#> 6 0.00000000
```
