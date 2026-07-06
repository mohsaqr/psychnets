# Node centrality

Node centrality

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
  networks), `"betweenness"`, `"closeness"`, plus any names supplied via
  `centrality_fn`. Betweenness and closeness are computed on the
  absolute, inverted-weight graph and are not generally meaningful on
  signed networks – request them only when a downstream comparison needs
  them.

- centrality_fn:

  Optional function taking the weighted adjacency matrix and returning a
  named list of node-centrality vectors, used to supply any `measures`
  not built in.

- ...:

  Unused.

## Value

A tidy `data.frame`, one row per node, with a `node` column and one
column per requested measure (`strength` = sum of absolute edge weights,
`expected_influence` = sum of signed edge weights, by default).

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
