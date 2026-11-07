# Edge betweenness centrality

For each edge, the number of weighted shortest paths (across all node
pairs) that pass through it - a high value marks an edge that bridges
otherwise distant parts of the network. When a pair is joined by several
equally short paths each counts as its share
(`1 / number of geodesics`), so the value is the raw, unnormalized edge
betweenness (igraph's `edge_betweenness()`), not a proportion: on the
path `a-b-c-d-e` the four edges score 4, 6, 6, 4. Geodesics are computed
on inverse absolute weights, so strong edges count as short, matching
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)'s
node betweenness/closeness.

## Usage

``` r
net_edge_betweenness(x, invert = TRUE, labels = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix. A
  `psychnet` is directed when its `directed` field says so; a bare
  matrix is treated as directed when it is not symmetric
  ([`isSymmetric()`](https://rdrr.io/r/base/isSymmetric.html)), with
  `x[i, j]` the edge `i -> j`.

- invert:

  If `TRUE` (default) edge weights are inverted to distances (strong
  association = short path). Set `FALSE` to treat weights as distances.

- labels:

  Optional node labels (used when `x` is a bare matrix): one unique
  label per node.

## Value

A tidy `data.frame`, one row per edge: `from`, `to`, `edge_betweenness`
(count of shortest paths through the edge). Undirected networks give one
row per unordered edge (counting each unordered node pair once);
directed networks one row per arc, counting ordered pairs.

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
