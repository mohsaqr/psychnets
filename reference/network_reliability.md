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
# \donttest{
network_reliability(SRL_Claude)
#> # split-half reliability: glasso | 100 iterations (50/50 split)
#>           metric       mean         sd      lower     upper
#> 1   mean_abs_dev 0.06693026 0.02270286 0.03359414 0.1253135
#> 2 median_abs_dev 0.05931370 0.02360581 0.02751320 0.1156776
#> 3    correlation 0.98010564 0.01314679 0.94539317 0.9944821
#> 4    max_abs_dev 0.15698277 0.04968053 0.07990015 0.2594862
# }
```
