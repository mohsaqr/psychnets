# Relative-importance network (LMG / Shapley)

Builds a directed network in which the edge predictor -\> outcome is the
predictor's LMG (Shapley) share of the outcome node's regression
R-squared. Equivalent in purpose to
`relaimpo::calc.relimp(type = "lmg")` applied nodewise / `bootnet`'s
`"relimp"` default, pure base R and self-certified via
[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md).

## Usage

``` r
relimp_network(
  data = NULL,
  cor_matrix = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  max_nodes = 21L,
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

- cor_method:

  Correlation when `data` is supplied: `"pearson"` (default),
  `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- max_nodes:

  Refuse to run above this many nodes (the cost grows as `2^(p-1)` per
  node). Default 21.

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default)
  or `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the directed importance matrix
(`weights[k, j]` = importance of `k` for outcome `j`), with `$r2`
(per-node full-model R-squared), `$cor_matrix`, and `$kkt` (the
decomposition residual).

## Examples

``` r
S <- 0.4^abs(outer(1:5, 1:5, "-"))
relimp_network(cor_matrix = S)
#> <psychnet> relimp network
#>   nodes: 5   edges: 20   (directed)
#>   optimality (KKT residual): 5.55e-17
```
