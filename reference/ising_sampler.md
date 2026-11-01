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

  Optional non-negative observation weights, one per input row after
  missing-data preparation. When `min_sum` is used, weights are filtered
  alongside their corresponding rows. `NULL` (default) is unweighted.

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
and `$kkt` (worst nodewise score residual). If Wald pruning changes any
coefficient, `$kkt` is `NA` for the altered fit; `$fit_kkt` retains the
unpruned residual and `$certificate_target` identifies it as
`"pre_pruning_fit"`.

## Examples

``` r
b <- dichotomize(SRL_GPT)
ising_sampler(b)                 # every edge, certified ML fit
#> <psychnet> ising_sampler network
#>   nodes: 5   edges: 10   (undirected)
#>   optimality (KKT residual): 1.07e-08
pruned <- ising_sampler(b, alpha = 0.05)
as.data.frame(pruned)
#>   from to    weight
#> 1  CSU IV 2.4336970
#> 2  CSU SE 1.9776719
#> 3   IV SE 1.5285878
#> 4  CSU SR 1.7964616
#> 5   IV SR 2.0903313
#> 6   SE SR 0.9333181
certificate(pruned)              # NA: pruning altered the ML solution
#>          method certificate kind certified
#> 1 ising_sampler          NA  kkt        NA
```
