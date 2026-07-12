# Centrality-stability coefficient (case-dropping subset bootstrap)

Centrality-stability coefficient (case-dropping subset bootstrap)

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
  ...
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations).

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

- ...:

  Passed to the estimator.

## Value

An object of class `psychnet_stability` with `$cs` (CS-coefficient per
measure) and a tidy `$table` (columns `measure`, `drop_prop`,
`mean_cor`, `sd_cor`, `prop_above`) of the case-dropping correlations by
drop proportion. Visualise it with
[`plot.psychnet_stability()`](https://pak.dynasite.org/psychnets/reference/plot.psychnet_stability.md).

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
cs <- net_stability(x, drop_prop = c(0.3, 0.6), iter = 10)
cs$cs
#>           strength expected_influence 
#>                0.3                0.3 
```
