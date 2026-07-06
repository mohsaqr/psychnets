# Unregularized Ising network for binary data

Estimates an Ising model by *unpenalized* nodewise logistic regression,
with optional Wald p-value edge pruning, combined by the AND (default)
or OR rule. The unregularized counterpart of
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md);
self-certified by the maximum-likelihood score residual (see
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
at `lambda = 0`).

## Usage

``` r
ising_sampler(
  data,
  rule = c("AND", "OR"),
  alpha = NULL,
  adjust = "none",
  min_sum = NULL,
  weights = NULL,
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  Binary (0/1) data frame or matrix (rows = observations).

- rule:

  Edge-combination rule: `"AND"` (default) or `"OR"`.

- alpha:

  Significance level for Wald edge pruning; `NULL` (default) keeps every
  edge.

- adjust:

  Multiple-comparison adjustment for the edge p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

- min_sum:

  Minimum row sum-score; rows below it are dropped before fitting.
  `NULL` (default) keeps every row.

- weights:

  Optional non-negative observation weights, one per retained row.
  `NULL` (default) is unweighted.

- na_method:

  Missing-data handling: `"pairwise"` (default, mode-impute) or
  `"listwise"`. See
  [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md).

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the symmetric weight matrix,
with `$thresholds` (node intercepts), `$rule`, `$p_values`, `$nodewise`
(for
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)),
and `$kkt` (worst nodewise score residual).

## Examples

``` r
set.seed(1)
z <- matrix(stats::rnorm(500 * 2), 500, 2)
x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(500 * 4), 500)
b <- (x > 0) * 1L
colnames(b) <- paste0("V", 1:4)
ising_sampler(b)
#> <psychnet> ising_sampler network
#>   nodes: 4   edges: 6   (undirected)
#>   optimality (KKT residual): 7.03e-10
```
