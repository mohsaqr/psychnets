# Correlation network

Marginal (zero-order) association network: the Pearson correlation
matrix with the diagonal removed. Equivalent to `bootnet`'s `"cor"`
default.

## Usage

``` r
cor_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  threshold = 0,
  alpha = NULL,
  adjust = "none",
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations). Optional if
  `cor_matrix` is supplied.

- cor_matrix:

  Optional precomputed correlation matrix; if given, `data` is ignored
  and `n` is required when `alpha` is used.

- n:

  Sample size (needed for significance testing when `cor_matrix` is
  supplied).

- cor_method:

  Correlation method: `"pearson"` (default), `"spearman"`, `"kendall"`,
  or `"auto"` (polychoric/polyserial for ordinal items, the
  [`qgraph::cor_auto`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
  default; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- threshold:

  Correlations with absolute value below this are set to zero. Default
  0.

- alpha:

  Significance level; if set, correlations not significant at `alpha`
  are zeroed. `NULL` (default) keeps every edge.

- adjust:

  Multiple-comparison adjustment for the edge p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

- na_method:

  Missing-data handling: `"pairwise"` (default) uses pairwise-complete
  correlations projected to the nearest positive-definite matrix;
  `"listwise"` drops rows with any `NA`. Identical when data is
  complete.

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the thresholded correlation
matrix, with `$cor_matrix`, `$n_eff`, `$na_method` (and `$p_values` when
`alpha` is used).

## Examples

``` r
x <- matrix(stats::rnorm(200 * 4), 200, 4)
cor_network(x)
#> <psychnet> cor network
#>   nodes: 4   edges: 6   (undirected)
cor_network(x, alpha = 0.05, adjust = "BH")
#> <psychnet> cor network
#>   nodes: 4   edges: 0   (undirected)
```
