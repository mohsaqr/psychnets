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
# \donttest{
casedrop_reliability(SRL_Claude)
#> # edge-weight stability: glasso | CS = 0.90 (spearman cor >= 0.70 at 95%)
#>            metric drop_prop       mean          sd
#> 1    mean_abs_dev       0.1 0.01210425 0.006842117
#> 2    mean_abs_dev       0.2 0.01890183 0.009083151
#> 3    mean_abs_dev       0.3 0.02149010 0.011285685
#> 4    mean_abs_dev       0.4 0.02698512 0.011699679
#> 5    mean_abs_dev       0.5 0.03276844 0.015423245
#> 6    mean_abs_dev       0.6 0.04346120 0.019907911
#> 7    mean_abs_dev       0.7 0.05458419 0.019680164
#> 8    mean_abs_dev       0.8 0.07213005 0.028998661
#> 9    mean_abs_dev       0.9 0.11435088 0.041262261
#> 10 median_abs_dev       0.1 0.01087812 0.007526945
#> 11 median_abs_dev       0.2 0.01662288 0.009894813
#> 12 median_abs_dev       0.3 0.01908751 0.012849561
#> 13 median_abs_dev       0.4 0.02328488 0.012107749
#> 14 median_abs_dev       0.5 0.02834539 0.016470417
#> 15 median_abs_dev       0.6 0.03893692 0.021357956
#> 16 median_abs_dev       0.7 0.04928713 0.022370871
#> 17 median_abs_dev       0.8 0.06351304 0.032443227
#> 18 median_abs_dev       0.9 0.10422689 0.046986840
#> 19    correlation       0.1 0.97948093 0.018021434
#> 20    correlation       0.2 0.97726341 0.017660241
#> 21    correlation       0.3 0.97759688 0.018398626
#> 22    correlation       0.4 0.97386763 0.017045752
#> 23    correlation       0.5 0.96892625 0.019271040
#> 24    correlation       0.6 0.96668355 0.019234762
#> 25    correlation       0.7 0.95693658 0.020374142
#> 26    correlation       0.8 0.95491323 0.024149980
#> 27    correlation       0.9 0.91966867 0.047363394
#> 28    max_abs_dev       0.1 0.02782703 0.013178779
#> 29    max_abs_dev       0.2 0.04431631 0.019171457
#> 30    max_abs_dev       0.3 0.05063819 0.022210782
#> 31    max_abs_dev       0.4 0.06527754 0.027760050
#> 32    max_abs_dev       0.5 0.07769319 0.031758263
#> 33    max_abs_dev       0.6 0.09759017 0.038810847
#> 34    max_abs_dev       0.7 0.12589880 0.043891748
#> 35    max_abs_dev       0.8 0.16921224 0.061437261
#> 36    max_abs_dev       0.9 0.26500688 0.104004806
# }
```
