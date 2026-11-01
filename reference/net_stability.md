# Centrality-stability coefficient (case-dropping subset bootstrap)

Asks how far the sample can shrink before the ranking of node
centralities changes. For each drop proportion the network is
re-estimated on `iter` random subsets without replacement, and each
subset's centralities are rank-correlated (Spearman) with the
full-sample centralities. The CS-coefficient of a measure is the largest
drop proportion at which that correlation stays at or above `threshold`
with probability at least `certainty`; Epskamp, Borsboom and Fried
(2018) suggest that it should not fall below 0.25 and preferably exceed
0.5 before centrality differences are interpreted.

## Usage

``` r
net_stability(
  data,
  method = "glasso",
  measures = c("strength", "expected_influence"),
  centrality_fn = NULL,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  iter = 100L,
  threshold = 0.7,
  certainty = 0.95,
  labels = NULL,
  estimator_args = list(),
  ...
)
```

## Arguments

- data:

  Data frame or matrix (rows = observations), resampled exactly as given
  (see
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)).

- method:

  Estimator (see
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)).
  Default `"glasso"`.

- measures:

  Centrality measures to assess. Defaults to the two recommended for
  psychometric networks (`c("strength", "expected_influence")`);
  `"betweenness"`/`"closeness"` and custom measures (via
  `centrality_fn`) are also accepted. See
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).

- centrality_fn:

  Optional function supplying any non-built-in `measures` (see
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)).

- drop_prop:

  Proportions of cases to drop. Default `seq(0.1, 0.9, 0.1)`.

- iter:

  Subsets per proportion. Default 100.

- threshold:

  Minimum acceptable rank correlation. Default 0.7.

- certainty:

  Probability the correlation must exceed `threshold`. Default 0.95.

- labels:

  Optional node labels.

- estimator_args:

  Named list of estimator arguments. Use this for names consumed by the
  stability diagnostic itself, such as estimator `threshold`.

- ...:

  Passed to the estimator.

## Value

An object of class `psychnet_stability` with `$cs` (CS-coefficient per
measure) and a tidy `$table` (columns `measure`, `drop_prop`,
`mean_cor`, `sd_cor`, `prop_above`) of the case-dropping correlations by
drop proportion. Visualise it with
[`plot.psychnet_stability()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_stability.md).
Undefined correlations count as failing the threshold in `prop_above`
and the CS-coefficient; they are omitted from `mean_cor` and `sd_cor`.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods*, 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
# method = "pcor" and small iter keep this example fast; the "glasso"
# default and iter = 100 are what a real analysis should use.
cs <- net_stability(x, method = "pcor", drop_prop = c(0.3, 0.6), iter = 10)
cs          # the CS-coefficient of each measure
#> <psychnet_stability> pcor, 10 subsets/proportion
#>   CS-coefficient (cor >= 0.70 with 95% certainty):
#>     strength             0.00
#>     expected_influence   0.00
plot(cs)    # mean correlation by drop proportion
```
