# Detect redundant node pairs ("goldbricker")

Flags pairs of items that behave redundantly in a network: their
correlations with all other items are mostly statistically
indistinguishable (a small proportion of significantly different
correlations) and the two items are themselves strongly correlated.
Ported from
[`networktools::goldbricker`](https://rdrr.io/pkg/networktools/man/goldbricker.html).

## Usage

``` r
redundancy(
  data,
  p = 0.05,
  threshold = 0.25,
  cor_min = 0.5,
  cor_method = c("auto", "pearson", "spearman", "kendall")
)
```

## Arguments

- data:

  A numeric data frame or matrix (rows = observations).

- p:

  Significance level for each pairwise correlation-difference test.
  Default 0.05.

- threshold:

  Maximum proportion of significantly different correlations for a pair
  to be flagged redundant. Default 0.25.

- cor_min:

  Minimum correlation between the two items themselves. Default 0.5.

- cor_method:

  Correlation type: `"auto"` (default,
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md) -
  polychoric/polyserial as appropriate, matching goldbricker),
  `"pearson"`, `"spearman"`, or `"kendall"`.

## Value

A tidy `data.frame` (class `psychnet_redundancy`), one row per flagged
pair (sorted most-redundant first), with columns `item1`, `item2`,
`proportion` (share of significantly different correlations) and
`correlation`. Zero rows when nothing is flagged. The full proportion
matrix is in `attr(x, "proportion_matrix")`.

## References

Hallquist, M. N., Wright, A. G. C., & Molenaar, P. C. M. (2021).
Problems with centrality measures in psychopathology networks.
*Multivariate Behavioral Research*, 56(2), 199-223.

## Examples

``` r
redundancy(SRL_Claude)
#> # redundant pairs (proportion < 0.25, r > 0.50): 0 found
#>   none
```
