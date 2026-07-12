# Split-half reliability of the network edge structure

Repeatedly splits the sample into two halves, estimates a network on
each, and compares their edge-weight vectors. Reports, across splits,
the edge-weight correlation between halves plus the mean/median/maximum
absolute edge deviation - a psychometric reliability view of the
estimated structure.

## Usage

``` r
network_reliability(
  data,
  method = "glasso",
  iter = 100L,
  split = 0.5,
  cor_method = c("pearson", "spearman", "kendall"),
  labels = NULL,
  ...
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations), or a
  `psychnet_group` (split-half per level).

- method:

  Estimator (see
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)).
  Default `"glasso"`.

- iter:

  Number of split-half iterations. Default 100.

- split:

  Fraction of rows in the first half. Default 0.5.

- cor_method:

  Correlation method for the between-halves edge comparison: `"pearson"`
  (default), `"spearman"`, or `"kendall"`.

- labels:

  Optional node labels.

- ...:

  Passed to the estimator.

## Value

A tidy `data.frame` (class `psychnet_reliability`), one row per metric
with columns `metric`, `mean`, `sd`, `lower`, `upper`. The per-split
draws are carried in `attr(x, "iterations")` for
[`plot.psychnet_reliability()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_reliability.md).

## Examples

``` r
# `iter` is kept small here so the example runs quickly; the default
# (iter = 100) is what a real reliability assessment should use.
network_reliability(SRL_Claude, iter = 10)
#> # split-half reliability: glasso | 10 iterations (50/50 split)
#>           metric       mean         sd      lower      upper
#> 1   mean_abs_dev 0.06505893 0.02720464 0.02471090 0.10065454
#> 2 median_abs_dev 0.05980563 0.02943647 0.01557087 0.09904948
#> 3    correlation 0.98131883 0.01347642 0.95790513 0.99665236
#> 4    max_abs_dev 0.15438705 0.06804919 0.06582753 0.26122348
```
