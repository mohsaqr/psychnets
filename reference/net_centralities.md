# Compute the centrality of each node in a network

`net_centralities()` computes the centrality of each node in a network,
by default its strength and expected influence. Centrality indices
summarize the position of each node in a network; for psychometric
networks, whose edges are signed and weighted, the most informative are
strength, the sum of the absolute weights of a node's edges, and
expected influence, the sum of their signed weights, in which negative
edges offset positive ones (Robinaugh, Millner & McNally, 2016). A node
with high strength is strongly connected to the other variables;
expected influence also shows whether these connections are, on balance,
positive or negative.

On request the function adds two-step expected influence, which also
counts the expected influence a node passes on through its neighbours,
and betweenness and closeness, which are computed on shortest paths
through the absolute edge weights, inverted to distances (Opsahl,
Agneessens & Skvoretz, 2010). Because betweenness and closeness ignore
the sign of an edge, they are less suited to signed networks. Further
measures can be supplied as a function (`centrality_fn`). For a directed
network, such as one from
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md),
strength and expected influence sum a node's outgoing edges.

It returns one row per node, with one column per requested measure.

Whether differences in centrality can be interpreted depends on their
accuracy and stability; see
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
and
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md).
Centrality with respect to communities is computed by
[`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md),
and the betweenness of edges by
[`net_edge_betweenness()`](https://pak.dynasite.org/psychnets/reference/net_edge_betweenness.md).

## Usage

``` r
net_centralities(
  x,
  measures = c("strength", "expected_influence"),
  centrality_fn = NULL,
  ...
)

# S3 method for class 'psychnet_centrality'
plot(
  x,
  type = c("bar", "line"),
  scale = c("raw", "z", "relative"),
  measures = NULL,
  ...
)
```

## Arguments

- x:

  A fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)),
  a weighted adjacency matrix, or a
  [psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  (one table per group). For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), the table
  returned by `net_centralities()`.

- measures:

  Measures to compute: any of `"strength"`, `"expected_influence"` (the
  defaults), `"expected_influence_2step"`, `"betweenness"`,
  `"closeness"`, and names supplied by `centrality_fn`. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the measure
  columns to draw (default: all).

- centrality_fn:

  Optional function that takes the weighted adjacency matrix (with node
  labels as dimnames) and returns a named list of centrality vectors,
  one numeric value per node, for measures that are not built in. A
  named vector is matched to the nodes by label; an unnamed one is taken
  in node order.

- ...:

  Unused.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `"bar"`
  (default; one sorted panel per measure) or `"line"` (one panel per
  measure, with the nodes in a shared order).

- scale:

  For `plot(type = "line")`: `"raw"` (default), `"z"` (z-scores per
  measure) or `"relative"` (rescaled to 0-1 per measure). `type = "bar"`
  always shows raw values.

## Value

A data frame of class `psychnet_centrality`, one row per node, with a
`node` column and one column per requested measure, named as in
`measures`. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws it and returns it invisibly. For a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md),
one table per group (see
[psychnet_result_group](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md)).

## Details

Two-step expected influence is the node's expected influence plus the
sum, over its neighbours, of each edge weight times that neighbour's
expected influence.

Betweenness and closeness are normalized, as in igraph's
`betweenness(normalized = TRUE)` and `closeness(normalized = TRUE)`:
betweenness is divided by the number of pairs of other nodes,
`(p - 1)(p - 2) / 2`, and closeness is the number of reachable nodes
divided by the sum of the distances to them (the inverse mean distance).
They therefore differ from
[`qgraph::centrality_auto()`](https://rdrr.io/pkg/qgraph/man/centrality_auto.html),
which reports betweenness as a count of paths and closeness as
`1 / sum(distances)`.

## Functions

- `plot(psychnet_centrality)`: One panel per measure, nodes sorted by
  value (`type = "bar"`) or in a shared order (`type = "line"`).

## Conditions

`psychnet_bad_input` when `x` is not a fitted network or a square
numeric matrix, a measure is not built in and `centrality_fn` does not
supply it, `centrality_fn` returns something other than a named list of
one numeric value per node (or names that are not the node labels), or a
group holds within- and between-person networks.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) raises
`psychnet_bad_input` for a `measures` entry that is not a measure
column, `psychnet_plot_skipped` (warning) when it skips a measure
without finite values, and `psychnet_nothing_to_plot` (an error) when no
requested measure has one.

## References

Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016). Identifying
highly influential nodes in the complicated grief network. *Journal of
Abnormal Psychology*, 125(6), 747-757.
[doi:10.1037/abn0000181](https://doi.org/10.1037/abn0000181)

Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in
weighted networks: Generalizing degree and shortest paths. *Social
Networks*, 32(3), 245-251.
[doi:10.1016/j.socnet.2010.03.006](https://doi.org/10.1016/j.socnet.2010.03.006)

## See also

[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
and
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
for the accuracy and stability of centralities;
[`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md);
[`net_edge_betweenness()`](https://pak.dynasite.org/psychnets/reference/net_edge_betweenness.md);
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
for node predictability.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
cent <- net_centralities(ebic_glasso(cor_matrix = S, n = 250))
cent
#>   node  strength expected_influence
#> 1   V1 0.3681824          0.3681824
#> 2   V2 0.7105013          0.7105013
#> 3   V3 0.6846378          0.6846378
#> 4   V4 0.6846378          0.6846378
#> 5   V5 0.7105013          0.7105013
#> 6   V6 0.3681824          0.3681824
plot(cent)

plot(cent, type = "line")
```
