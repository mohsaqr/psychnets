# Compute the edge betweenness of a network

`net_edge_betweenness()` computes the betweenness of each edge in a
network. Edge betweenness describes the role of an edge in connecting
the network: it counts how many of the shortest paths between pairs of
nodes run through the edge. An edge with high betweenness links parts of
the network that would otherwise be far apart, for example two groups of
closely related symptoms.

The function finds shortest paths on the absolute edge weights, inverted
to distances so that a strong association is a short step (with
`invert = FALSE` the weights are used as distances). When two nodes are
joined by several equally short paths, each path counts as its share, 1
/ (number of such paths).

It returns one row per edge. The value is a count, as returned by
[`igraph::edge_betweenness()`](https://r.igraph.org/reference/betweenness.html),
not a proportion: on the path a-b-c-d-e the four edges score 4, 6, 6 and
4.

For the betweenness of nodes, use
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md);
for the connections of nodes to predefined communities, use
[`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md).

## Usage

``` r
net_edge_betweenness(x, invert = TRUE, labels = NULL)
```

## Arguments

- x:

  A fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
  or a square weighted adjacency matrix. A bare matrix that is not
  symmetric is treated as directed, the entry in row i and column j
  being the edge from node i to node j.

- invert:

  If `TRUE` (default), edge weights are inverted to distances (strong
  association = short path); if `FALSE`, weights are distances.

- labels:

  Optional node labels for a bare matrix, one unique label per node.

## Value

A data frame with one row per edge and columns `from`, `to` and
`edge_betweenness` (the number of shortest paths through the edge). An
undirected network has one row per edge, counting each pair of nodes
once; a directed network has one row per arc, counting ordered pairs.

## Conditions

`psychnet_bad_input` when `x` is not a fitted network or a square
numeric matrix, contains non-finite weights, or `labels` does not give
one unique label per node.

## See also

[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
for node betweenness and other measures;
[`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md)
for bridge centrality.

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
