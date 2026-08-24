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
net_casedrop_reliability(
  data,
  method = "glasso",
  drop_prop = seq(0.1, 0.9, by = 0.1),
  iter = 100L,
  threshold = 0.7,
  certainty = 0.95,
  cor_method = c("spearman", "pearson", "kendall"),
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
  or a `psychnet_group` (case-dropped per level).

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

- estimator_args:

  Named list of estimator arguments. Use this for names consumed by the
  diagnostic itself, such as estimator `threshold`.

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
# `method`, `iter`, and `drop_prop` are chosen so the example runs quickly;
# the defaults (method = "glasso", iter = 100, drop_prop = seq(0.1, 0.9,
# 0.1)) are what a real reliability assessment should use.
net_casedrop_reliability(SRL_Claude, method = "pcor", iter = 5,
                         drop_prop = c(0.25, 0.5))
#> # edge-weight stability: pcor | CS = 0.50 (spearman cor >= 0.70 at 95%)
#>           metric drop_prop       mean          sd
#> 1   mean_abs_dev      0.25 0.02624533 0.005282477
#> 2   mean_abs_dev      0.50 0.03307525 0.011005165
#> 3 median_abs_dev      0.25 0.02232825 0.006143106
#> 4 median_abs_dev      0.50 0.03082581 0.013881714
#> 5    correlation      0.25 0.97818182 0.010141334
#> 6    correlation      0.50 0.96363636 0.030903149
#> 7    max_abs_dev      0.25 0.05938727 0.016592220
#> 8    max_abs_dev      0.50 0.07034004 0.023615060
```
