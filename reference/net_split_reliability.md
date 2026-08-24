# Split-half reliability of the network edge structure

Repeatedly splits the sample into two halves, estimates a network on
each, and compares their edge-weight vectors. Reports, across splits,
the edge-weight correlation between halves plus the mean/median/maximum
absolute edge deviation - a psychometric reliability view of the
estimated structure.

## Usage

``` r
net_split_reliability(
  data,
  method = "glasso",
  iter = 100L,
  split = 0.5,
  cor_method = c("pearson", "spearman", "kendall"),
  labels = NULL,
  estimator_args = list(),
  ...
)
```

## Arguments

- data:

  Data frame or matrix (rows = observations), resampled exactly as given
  (see
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)),
  or a `psychnet_group` (split-half per level).

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

- estimator_args:

  Named list of estimator arguments. Use this for names consumed by the
  diagnostic itself, such as estimator `threshold`.

- ...:

  Passed to the estimator.

## Value

A tidy `data.frame` (class `psychnet_reliability`), one row per metric
with columns `metric`, `mean`, `sd`, `lower`, `upper`. The per-split
draws are carried in `attr(x, "iterations")` for
[`plot.psychnet_reliability()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_reliability.md).

## Examples

``` r
# `method` and `iter` are chosen so the example runs quickly; the defaults
# (method = "glasso", iter = 100) are what a real assessment should use.
net_split_reliability(SRL_Claude, method = "pcor", iter = 10)
#> # split-half reliability: pcor | 10 iterations (50/50 split)
#>           metric       mean         sd      lower      upper
#> 1   mean_abs_dev 0.07808955 0.02457258 0.03846060 0.10963050
#> 2 median_abs_dev 0.06588545 0.02220565 0.02991325 0.09716984
#> 3    correlation 0.97498013 0.01592382 0.94570362 0.99516211
#> 4    max_abs_dev 0.19047875 0.07735301 0.07527392 0.32270329
```
