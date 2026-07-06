# Triangulated Maximally Filtered Graph (TMFG)

Builds a sparse, planar, chordal association network by greedily
retaining the `3(p - 2)` most informative edges (Massara et al. 2016).
Equivalent in purpose to `NetworkToolbox::TMFG()` / `bootnet`'s `"TMFG"`
default, pure base R; correctness is certified structurally by
[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md).

## Usage

``` r
tmfg_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations). Optional if
  `cor_matrix` is supplied.

- cor_matrix:

  Optional correlation matrix.

- n:

  Accepted and ignored. TMFG is a structural filter and needs no sample
  size; the argument exists only so a uniform `(cor_matrix=, n=)` call
  shared with the other estimators does not partial-match `na_method`.

- cor_method:

  Correlation when `data` is supplied: `"pearson"` (default),
  `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default)
  or `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the filtered (signed)
correlation matrix on the retained edges, with `$adjacency`, `$cliques`,
`$separators` (the chordal decomposition used by
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)),
and `$cor_matrix`.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 6), 200, 6)
tmfg_network(x)
#> <psychnet> tmfg network
#>   nodes: 6   edges: 12   (undirected)
```
