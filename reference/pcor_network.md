# Estimate a partial correlation network

`pcor_network()` estimates a partial correlation network, in which each
pair of variables is joined by their correlation after controlling for
all other variables in the network. An edge is therefore an association
that the remaining variables do not account for, and for normally
distributed data a partial correlation of zero means that the two
variables are conditionally independent.

This is the unregularized Gaussian graphical model, in which every
partial correlation is estimated and none is set to zero by the model
itself. The partial correlations are computed from the inverse of the
correlation matrix of raw data or of a supplied correlation matrix.
`alpha` keeps only partial correlations that differ significantly from
zero, with an optional multiple-comparison adjustment (`adjust`), and
`threshold` removes weak ones.

`pcor_network()` returns the partial correlations as a network; the
network involves no optimization and has no certificate.

With many variables relative to the sample size, partial correlations
are estimated with large error, and with fewer observations than
variables they are not identified.
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
then estimate a sparse network instead;
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md)
gives the marginal correlations. The network is the `"pcor"` estimator
of `bootnet`.

## Usage

``` r
pcor_network(
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
whose edge weights are partial correlations.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports that the network has no certificate. The object also stores the
precision matrix (`precision`), the correlation matrix (`cor_matrix`),
the sample size (`n_eff`), the missing-data handling used (`na_method`),
the observations per pair when estimated from `data` (`n_pair`) and,
with `alpha`, the unadjusted p-values (`p_values`). Nodes keep the
column order of the input.

## Details

With \\\Theta\\ the inverse correlation matrix, the partial correlation
of variables \\i\\ and \\j\\ is \\-\Theta\_{ij} /
\sqrt{\Theta\_{ii}\Theta\_{jj}}\\. The significance test is described in
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md);
for \\p\\ variables it has \\n - p\\ degrees of freedom. It needs a
single sample size, so with pairwise data whose pairs have different
numbers of observations `alpha` is refused.

When the correlation matrix is singular, as with fewer observations than
variables or with collinear or duplicate columns, the nearest
positive-definite matrix is inverted so that a network is still
returned, and a warning says that its edges are not interpretable.

## Conditions

`psychnet_singular` (warning) when the correlation matrix is singular or
nearly so; `psychnet_degenerate_sample` (error, also of class
`psychnet_bad_input`) when the sample is too small for the test (\\n \le
p\\); `psychnet_bad_input` when an argument is out of range, `threshold`
is logical, `alpha` is combined with Kendall correlations or with
unequal pairwise sample sizes, `cor_matrix` is not a finite, symmetric,
positive semi-definite matrix with a positive diagonal, `n` is missing
with `cor_matrix` and `alpha`, `data` has fewer than two usable
variables, or `labels` do not match the nodes;
`psychnet_constant_variable` (warning) when a constant or all-missing
column is dropped.

## See also

[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md);
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
for sparse networks;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "pcor"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

## Examples

``` r
x <- matrix(stats::rnorm(200 * 4), 200, 4)
pcor_network(x)
#> <psychnet> pcor network
#>   nodes: 4   edges: 6   (undirected)
pcor_network(x, alpha = 0.05, adjust = "holm")
#> <psychnet> pcor network
#>   nodes: 4   edges: 0   (undirected)
```
