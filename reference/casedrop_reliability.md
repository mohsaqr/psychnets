# Edge-weight stability coefficient (case-dropping subset bootstrap)

The edge-vector complement to
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md).
For each drop proportion the network is re-estimated on random
case-dropped subsets and the subset edge-weight vector is compared with
the full-sample one. The edge-weight CS-coefficient is the largest drop
proportion at which the edge-vector correlation stays `>= threshold`
with probability `>= certainty` (Epskamp, Borsboom & Fried 2018).

## Usage

``` r
casedrop_reliability(
  data,
  method = "glasso",
  drop_prop = seq(0.1, 0.9, by = 0.1),
  iter = 100L,
  threshold = 0.7,
  certainty = 0.95,
  cor_method = c("spearman", "pearson", "kendall"),
  labels = NULL,
  ...
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations), or a
  `psychnet_group` (case-dropped per level).

- method:

  Estimator (see
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)).
  Default `"glasso"`.

- drop_prop:

  Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.

- iter:

  Subsets per proportion. Default 100.

- threshold:

  Minimum acceptable edge-vector correlation. Default 0.7.

- certainty:

  Probability the correlation must exceed `threshold`. Default 0.95.

- cor_method:

  Correlation method for the edge-vector comparison: `"spearman"`
  (default, robust to the wide range of edge weights), `"pearson"`, or
  `"kendall"`.

- labels:

  Optional node labels.

- ...:

  Passed to the estimator.

## Value

A tidy `data.frame` (class `psychnet_casedrop`), one row per metric per
drop proportion, with columns `metric`, `drop_prop`, `mean`, `sd`. The
edge-weight CS-coefficient is carried in `attr(x, "cs")` and shown when
the result is printed. Visualise it with
[`plot.psychnet_casedrop()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_casedrop.md).

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy. *Behavior Research Methods*,
50(1), 195-212.

## Examples

``` r
# `iter` and `drop_prop` are kept small here so the example runs quickly;
# the defaults (iter = 100, drop_prop = seq(0.1, 0.9, 0.1)) are what a real
# reliability assessment should use.
casedrop_reliability(SRL_Claude, iter = 5, drop_prop = c(0.25, 0.5))
#> # edge-weight stability: glasso | CS = 0.50 (spearman cor >= 0.70 at 95%)
#>           metric drop_prop       mean          sd
#> 1   mean_abs_dev      0.25 0.01887120 0.009120429
#> 2   mean_abs_dev      0.50 0.02206496 0.005387005
#> 3 median_abs_dev      0.25 0.01604977 0.008175395
#> 4 median_abs_dev      0.50 0.01955155 0.006003321
#> 5    correlation      0.25 0.97513583 0.008683012
#> 6    correlation      0.50 0.97512845 0.012960209
#> 7    max_abs_dev      0.25 0.04745924 0.024045972
#> 8    max_abs_dev      0.50 0.05437319 0.015481959
```
