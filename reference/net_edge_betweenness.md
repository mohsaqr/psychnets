# Edge betweenness centrality

For each edge, the share of weighted shortest paths (across all node
pairs) that pass through it - a high value marks an edge that bridges
otherwise distant parts of the network. Geodesics are computed on
inverse absolute weights, so strong edges count as short, matching
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)'s
node betweenness/closeness.

## Usage

``` r
net_edge_betweenness(x, invert = TRUE, labels = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix.

- invert:

  If `TRUE` (default) edge weights are inverted to distances (strong
  association = short path). Set `FALSE` to treat weights as distances.

- labels:

  Optional node labels (used when `x` is a bare matrix).

## Value

A tidy `data.frame`, one row per edge: `from`, `to`, `edge_betweenness`.
Undirected networks give one row per unordered edge.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_edge_betweenness(ebic_glasso(cor_matrix = S, n = 400))
#>   from to edge_betweenness
#> 1   V1 V2                5
#> 2   V2 V3                8
#> 3   V3 V4                9
#> 4   V4 V5                8
#> 5   V5 V6                5
```
