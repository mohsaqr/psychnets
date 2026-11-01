# Node centrality

Summarizes each node's position in a fitted network as a tidy table. For
psychometric networks the recommended measures are strength (the sum of
absolute edge weights) and expected influence (the signed sum, so
negative edges offset positive ones; Robinaugh, Millner & McNally 2016).
Betweenness and closeness follow Opsahl, Agneessens and Skvoretz (2010)
on the inverted absolute weights. For a directed network (e.g.
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md))
the measures sum the outgoing edges.

## Usage

``` r
net_centralities(
  x,
  measures = c("strength", "expected_influence"),
  centrality_fn = NULL,
  ...
)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object or a weighted adjacency matrix.

- measures:

  Character vector of measures to return. Any of `"strength"`,
  `"expected_influence"` (the defaults, recommended for psychometric
  networks), `"expected_influence_2step"` (a node's one-step expected
  influence plus the one-step expected influence it transmits through
  its neighbours), `"betweenness"`, `"closeness"`, plus any names
  supplied via `centrality_fn`. Betweenness and closeness are computed
  on the absolute, inverted-weight graph and are not generally
  meaningful on signed networks – request them only when a downstream
  comparison needs them.

- centrality_fn:

  Optional function taking the weighted adjacency matrix and returning a
  named list of node-centrality vectors, used to supply any `measures`
  not built in.

- ...:

  Unused.

## Value

A tidy `data.frame` of class `c("psychnet_centrality", "data.frame")`,
one row per node, with a `node` column and one column per requested
measure (`strength` = sum of absolute edge weights, `expected_influence`
= sum of signed edge weights, by default). Plot it with
[`plot.psychnet_centrality()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_centrality.md).
For a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/print.psychnet_group.md)
the result is a per-group list (see
[`print.psychnet_result_group()`](https://pak.dynasite.org/psychnets/reference/print.psychnet_result_group.md)).

## References

Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016). Identifying
highly influential nodes in the complicated grief network. *Journal of
Abnormal Psychology*, 125(6), 747-757.
[doi:10.1037/abn0000181](https://doi.org/10.1037/abn0000181)

Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in
weighted networks: Generalizing degree and shortest paths. *Social
Networks*, 32(3), 245-251.
[doi:10.1016/j.socnet.2010.03.006](https://doi.org/10.1016/j.socnet.2010.03.006)

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_centralities(ebic_glasso(cor_matrix = S, n = 250))
#>   node  strength expected_influence
#> 1   V1 0.3681824          0.3681824
#> 2   V2 0.7105013          0.7105013
#> 3   V3 0.6846378          0.6846378
#> 4   V4 0.6846378          0.6846378
#> 5   V5 0.7105013          0.7105013
#> 6   V6 0.3681824          0.3681824
```
