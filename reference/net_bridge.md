# Compute bridge centrality between communities of nodes

`net_bridge()` computes the bridge centrality of each node in a network
whose nodes belong to known communities, such as the symptoms of two
disorders or the items of different scales. Bridge centrality (Jones, Ma
& McNally, 2021) describes how strongly a node is connected to the
communities other than its own; nodes with high bridge centrality are
those through which the communities are linked, for example symptoms
that may connect two disorders in comorbidity.

The function computes five measures for each node of an undirected
network: bridge strength, the sum of the absolute weights of its edges
to other communities; one-step bridge expected influence (`bridge_ei1`),
the sum of the signed weights of those edges; two-step bridge expected
influence (`bridge_ei2`), which adds the signed weights of the paths of
two edges that end in another community; bridge betweenness, the number
of shortest paths between nodes of different communities that pass
through the node; and bridge closeness, the inverse of the mean distance
from the node to the nodes of other communities. Shortest paths use the
positive edges, with weights inverted to distances. The communities are
given by the user; `normalize` divides the measures by the number of
nodes they can involve.

It returns one row per node, with its community and the five measures.

For centrality without regard to communities, use
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).
To merge each community into a single node and estimate the network
between communities, use
[`net_aggregate()`](https://pak.dynasite.org/psychnets/reference/net_aggregate.md).

## Usage

``` r
net_bridge(x, communities, normalize = FALSE, labels = NULL)

# S3 method for class 'psychnet_bridge'
plot(x, measures = NULL, ...)
```

## Arguments

- x:

  An undirected fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
  or a symmetric weighted adjacency matrix. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the table
  returned by `net_bridge()`.

- communities:

  The community of each node, in one of three forms: an unnamed vector
  with one entry per node, in node order; a vector named by node label;
  or a list of node-label vectors, one per community, whose names are
  the community names, e.g.
  `list(anxiety = c("A1", "A2"), mood = c("M1", "M2"))`. Every node must
  belong to exactly one community, and there must be at least two.

- normalize:

  If `TRUE`, divide bridge strength and one-step bridge expected
  influence by the number `q` of nodes in other communities, bridge
  betweenness by `(p - 1) q`, and two-step bridge expected influence by
  `q + (p - 1)(q - q / (p - 1))`, where `p` is the number of nodes, as
  `networktools::bridge(normalize = TRUE)` does. Bridge closeness is
  unchanged. Default `FALSE`.

- labels:

  Optional node labels for a bare matrix, one unique label per node.

- measures:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  bridge columns to draw. Default: all five.

- ...:

  Unused.

## Value

A data frame of class `psychnet_bridge`, one row per node, with columns
`node`, `community`, `bridge_strength`, `bridge_betweenness`,
`bridge_closeness`, `bridge_ei1` and `bridge_ei2`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
sorted panel per measure, each node's point coloured and shaped by its
community, with a legend of the communities, and returns `x` invisibly.

## Functions

- `plot(psychnet_bridge)`: One sorted panel per bridge measure, with
  points coloured and shaped by community.

## Conditions

`psychnet_bad_input` when `x` is directed, is not a finite symmetric
square matrix, or `labels` is invalid; when `communities` lists a node
twice, names a node that is not in the network, leaves a node without a
community, has the wrong length, or defines fewer than two communities.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) raises
`psychnet_bad_input` for a `measures` entry that is not a bridge column,
`psychnet_plot_skipped` (warning) when it skips a measure without finite
values, and `psychnet_nothing_to_plot` (an error) when no requested
measure has one.

## References

Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality.
*Multivariate Behavioral Research*, 56(2), 353-367.

## See also

[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
for centrality without communities;
[`net_aggregate()`](https://pak.dynasite.org/psychnets/reference/net_aggregate.md)
for the network between communities.

## Examples

``` r
S <- 0.3^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 400)
br <- net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2))
br
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
plot(br)
```
