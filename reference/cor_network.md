# Estimate a correlation network

`cor_network()` estimates a correlation network, in which each pair of
variables is joined by their marginal (zero-order) correlation. An edge
shows that two variables covary without controlling for any other
variable, so it also reflects associations that run through other
variables in the network.

The correlations are computed from raw data, or taken from a supplied
correlation matrix. `cor_method` chooses Pearson, Spearman, Kendall or
polychoric and polyserial correlations; `alpha` keeps only correlations
that differ significantly from zero, with an optional
multiple-comparison adjustment (`adjust`); `threshold` removes weak
correlations.

`cor_network()` returns the correlation matrix as a network; the network
involves no optimization and has no certificate.

The network gives a descriptive overview of the association structure,
such as the clustering of items or symptoms, but it does not separate
direct from indirect associations. For associations that control for the
other variables use
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)
(all partial correlations) or
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
(a regularized partial-correlation network);
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md)
filters the correlations to a sparse backbone. The network is the
`"cor"` estimator of `bootnet`.

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

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given. Ordered factors enter as their
  level codes; other non-numeric columns are dropped with a message, and
  constant or all-missing columns with a warning.

- cor_matrix:

  A correlation matrix to use instead of `data`; a covariance matrix is
  rescaled to a correlation matrix. With `alpha`, `n` is then required.

- n:

  Sample size of `cor_matrix`, needed only for the significance test: a
  single number \>= 1, and a whole number \>= 3 when `alpha` is set.
  Ignored when `data` is given.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- threshold:

  Absolute edge weight below which an edge is set to zero: a single
  number \>= 0. Default `0`. A logical value (the `threshold = TRUE`
  rule of `qgraph`) is refused; use `alpha` for a significance rule.

- alpha:

  Significance level in (0, 1): edges not significant at this level are
  set to zero. Default `NULL`, no test. A Kendall test needs raw `data`;
  partial Kendall correlations cannot be tested.

- adjust:

  Multiple-comparison adjustment of the edge p-values, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`.

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) uses
  pairwise-complete correlations adjusted to the nearest
  positive-definite matrix; `"listwise"` drops incomplete rows. The two
  agree for complete data.

- labels:

  Optional node labels, one unique label per node of the network.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
whose edge weights are correlations.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports that the network has no certificate. The object also stores the
correlation matrix (`cor_matrix`), the sample size (`n_eff`; with
pairwise data, the median number of observations per pair), the
missing-data handling used (`na_method`; `"none"` with `cor_matrix`),
the observations per pair when estimated from `data` (`n_pair`) and,
with `alpha`, the unadjusted p-values (`p_values`). Nodes keep the
column order of the input.

## Details

Each correlation is tested against zero with the statistic \\t = r
\sqrt{df / (1 - r^2)}\\ on \\df = n - 2\\ degrees of freedom for a
correlation and \\df = n - p\\ for a full-order partial correlation of
\\p\\ variables. With `na_method = "pairwise"` and incomplete data, a
correlation is tested with the number of observations available for that
pair. Kendall correlations use a tie-corrected asymptotic test of the
raw data. The p-values are adjusted over the distinct pairs and are
computed from the correlations before `threshold` is applied; an edge is
kept only if it passes both rules.

## Conditions

`psychnet_degenerate_sample` (error, also of class `psychnet_bad_input`)
when the sample is too small for the test (\\df \le 0\\);
`psychnet_bad_input` when an argument is out of range, `threshold` is
logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
matrix with a positive diagonal, `n` is missing with `cor_matrix` and
`alpha`, a Kendall test is requested with `cor_matrix`, `data` has fewer
than two usable variables, a pair of variables is never observed
together, or `labels` do not match the nodes;
`psychnet_constant_variable` (warning) when a constant or all-missing
column is dropped. A `psychnet_dropped_columns` message names
non-numeric columns that were dropped.

## See also

[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md);
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md);
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md);
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "cor"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

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
