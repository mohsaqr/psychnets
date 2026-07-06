# Weighted clustering coefficients

Per-node clustering coefficients for a weighted (optionally signed)
network: Watts-Strogatz, Zhang-Horvath, Onnela, and Barrat. For signed
networks the signed variants are also returned (a closed triangle of
three negative or one negative edge lowers the signed coefficient).
Ported from
[`qgraph::clustcoef_auto`](https://rdrr.io/pkg/qgraph/man/clustcoef_auto.html).

## Usage

``` r
net_clustering(x, labels = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix.

- labels:

  Optional node labels (used when `x` is a bare matrix).

## Value

A tidy `data.frame` (class `psychnet_clustering`), one row per node:
`node`, `clustWS`, `signed_clustWS`, `clustZhang`, `signed_clustZhang`,
`clustOnnela`, `signed_clustOnnela`, `clustBarrat`.

## References

Costantini, G., & Perugini, M. (2014); Watts & Strogatz (1998); Zhang &
Horvath (2005); Onnela et al. (2005); Barrat et al. (2004).

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_clustering(ebic_glasso(cor_matrix = S, n = 400))
#>   node clustWS signed_clustWS clustZhang signed_clustZhang clustOnnela
#> 1   V1       0              0          0                 0           0
#> 2   V2       0              0          0                 0           0
#> 3   V3       0              0          0                 0           0
#> 4   V4       0              0          0                 0           0
#> 5   V5       0              0          0                 0           0
#> 6   V6       0              0          0                 0           0
#>   signed_clustOnnela clustBarrat
#> 1                  0           0
#> 2                  0           0
#> 3                  0           0
#> 4                  0           0
#> 5                  0           0
#> 6                  0           0
```
