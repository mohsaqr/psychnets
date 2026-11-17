# Detect redundant pairs of items

`redundancy()` detects pairs of items that measure nearly the same
thing, which appear in a network as two nodes with nearly the same
connections. Such topological overlap complicates the interpretation of
the network, and of centrality in particular (Hallquist, Wright &
Molenaar, 2021), so redundant pairs are commonly identified before a
network is estimated, and then combined or reduced to one item.

The function compares, for every pair of items, their correlations with
each of the other items, using a significance test for two dependent
correlations that share a variable. A pair is flagged as redundant when
the proportion of significantly different correlations is below
`threshold` and the two items correlate above `cor_min` with each other.
This is the procedure of `goldbricker()` in the `networktools` package.

It returns the flagged pairs, most redundant first.

To combine flagged items into one variable before estimating the
network,
[`net_aggregate()`](https://pak.dynasite.org/psychnets/reference/net_aggregate.md)
forms a composite score per group of items.

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

  A numeric data frame or matrix with one row per observation and one
  column per item.

- p:

  Significance level of each test of two correlations. Default `0.05`.

- threshold:

  A pair is flagged when its proportion of significantly different
  correlations is below this value. Default `0.25`.

- cor_min:

  Minimum correlation between the two items: a pair is flagged only when
  its correlation exceeds this value. Default `0.5`.

- cor_method:

  Correlation type: `"auto"` (default; polychoric or polyserial
  correlations where items are ordinal, see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)),
  `"pearson"`, `"spearman"` or `"kendall"`.

- x:

  A `psychnet_redundancy` object.

- ...:

  Unused.

## Value

A data frame of class `psychnet_redundancy`, one row per flagged pair,
most redundant first, with columns `item1`, `item2`, `proportion` (share
of significantly different correlations) and `correlation` (the
correlation of the two items); it has no rows when no pair is flagged.
The proportions for all pairs are kept with the result as the attribute
`proportion_matrix`. [`print()`](https://rdrr.io/r/base/print.html)
shows the flagged pairs, or a note that none were flagged, and returns
`x` invisibly.

## Details

Non-numeric columns are dropped, and so is any column that is constant
once incomplete rows are removed. Rows with a missing value in any
numeric column are dropped (listwise deletion, as in `goldbricker()`),
and the number of remaining rows is the sample size of every test.

## Functions

- `print(psychnet_redundancy)`: The thresholds used and the flagged
  pairs.

## Conditions

`psychnet_bad_input` for an invalid argument, for data that are not a
data frame or matrix, or for fewer than three usable items;
`psychnet_degenerate_sample` when fewer than four complete rows remain;
`psychnet_rows_dropped` (warning) when incomplete rows are dropped;
`psychnet_constant_variable` (warning) when a constant column is
dropped; `psychnet_dropped_columns` (a message) when non-numeric columns
are dropped.

## References

Hallquist, M. N., Wright, A. G. C., & Molenaar, P. C. M. (2021).
Problems with centrality measures in psychopathology networks.
*Multivariate Behavioral Research*, 56(2), 199-223.

## See also

[`net_aggregate()`](https://pak.dynasite.org/psychnets/reference/net_aggregate.md)
to combine items;
[`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)
for the correlations.

## Examples

``` r
redundancy(SRL_Claude)
#> # redundant pairs (proportion < 0.25, r > 0.50): 0 found
#>   none
```
