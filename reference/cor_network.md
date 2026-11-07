# Correlation network

Marginal (zero-order) association network: the correlation matrix of the
chosen `cor_method` (Pearson by default) with the diagonal removed.
Equivalent to `bootnet`'s `"cor"` default.

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
  `cor_matrix` is supplied. Ordered-factor columns enter as their
  integer level codes; other non-numeric columns are dropped with a
  `psychnet_dropped_columns` message, and constant / all-missing columns
  with a `psychnet_constant_variable` warning.

- cor_matrix:

  Optional precomputed correlation matrix; if given, `data` is ignored
  and `n` is required when `alpha` is used.

- n:

  Sample size (needed for significance testing when `cor_matrix` is
  supplied): a single number \>= 1, and a whole number of at least 3
  when `alpha` is set.

- cor_method:

  Correlation method: `"pearson"` (default), `"spearman"`, `"kendall"`,
  or `"auto"` (polychoric/polyserial for ordinal items, the
  [`qgraph::cor_auto`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
  default; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- threshold:

  Correlations with absolute value below this are set to zero. A single
  number \>= 0; default 0. (qgraph's logical `threshold = TRUE` is
  refused – use `alpha` for a significance rule.)

- alpha:

  Significance level in (0, 1); if set, correlations not significant at
  `alpha` are zeroed. `NULL` (default) keeps every edge. Kendall
  marginal correlations use a tie-corrected asymptotic test and require
  raw data. Significance testing of full-order partial Kendall
  correlations is not supported. A sample too small for the test (df =
  n - 2 - k \<= 0, with k the number of variables partialled out) is an
  error of class `psychnet_degenerate_sample`, never an empty network.

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

  Optional node labels, one unique label per network node.

## Value

A `psychnet` object whose `$weights` is the thresholded correlation
matrix, with `$cor_matrix`, `$n_eff`, `$na_method` (and `$p_values` when
`alpha` is used). `$n_pair` (the pairwise co-observation counts) is
present only when the network is estimated from `data`. Node order in
`$weights` and `$nodes` is always the column order of the input `data` /
`cor_matrix`, never sorted. Invalid arguments raise errors of class
`psychnet_bad_input`.

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
