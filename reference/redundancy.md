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

# S3 method for class 'psychnet_redundancy'
print(x, ...)
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

- x:

  A `psychnet_redundancy` object.

- ...:

  Unused.

## Value

A tidy `data.frame` (class `psychnet_redundancy`), one row per flagged
pair (sorted most-redundant first), with columns `item1`, `item2`,
`proportion` (share of significantly different correlations) and
`correlation`. Zero rows when nothing is flagged. The full proportion
matrix is stored with the result (attribute `proportion_matrix`).
[`print()`](https://rdrr.io/r/base/print.html) shows the flagged pairs,
or a note that none were flagged, and returns `x` invisibly.

## Details

**Data handling.** Non-numeric columns are dropped, and so is any column
that is constant once incomplete rows are removed. Rows with a missing
value in any numeric column are dropped (listwise deletion, as
`goldbricker` does), with a warning of class `psychnet_rows_dropped`
giving the count; the remaining number of rows is the `n` of every
correlation-difference test.

## Functions

- `print(psychnet_redundancy)`: The flagged pairs with the thresholds
  used.

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
