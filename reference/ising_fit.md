# Ising network for binary data

Estimates an Ising model by nodewise L1-penalized logistic regression
with EBIC selection, combined by the AND (default) or OR rule.
Equivalent in purpose to
[`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html),
but pure base R and self-certified: each node's regression reports its
stationarity (KKT) residual (see
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)).

## Usage

``` r
ising_fit(
  data,
  gamma = 0.25,
  rule = c("AND", "OR"),
  nlambda = 100L,
  lambda_min_ratio = 0.01,
  min_sum = NULL,
  weights = NULL,
  na_method = c("pairwise", "listwise"),
  native = TRUE,
  labels = NULL
)
```

## Arguments

- data:

  Binary data frame or matrix (rows = observations). Columns may be
  numeric 0/1, logical, or two-level factors (the first level is coded
  0). Any other column – a numeric column not coded 0/1 (e.g. 1/2), a
  factor with more than two levels, a character column – is an error
  naming it, never a silently dropped node. A constant column carries no
  information and is dropped.

- gamma:

  EBIC hyperparameter. Default 0.25.

- rule:

  Edge-combination rule: `"AND"` (default) or `"OR"`.

- nlambda:

  Number of penalties per nodewise path. Default 100.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest. With `native = FALSE`,
  omitting the argument retains glmnet's reference default for
  compatibility; an explicitly supplied value is honored.

- min_sum:

  Minimum row sum-score (number of endorsed items); rows below it are
  dropped before fitting. `NULL` (default) keeps every row.

- weights:

  Optional non-negative observation weights, one per input row after
  missing-data preparation. When `min_sum` is used, weights are filtered
  alongside their corresponding rows. Weights are rescaled to a maximum
  of 1 (the [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html)
  convention), so only their relative sizes matter: `weights = w` and
  `weights = 7 * w` give the same network. `NULL` (default) is
  unweighted.

- na_method:

  Missing-data handling: `"pairwise"` (default) single-imputes each
  column over its observed values (mode for binary), keeping the full
  sample; `"listwise"` drops incomplete rows. Identical for complete
  data.

- native:

  Solver switch. `TRUE` (default) uses psychnet's own pure-R,
  dependency-free, self-certified L1 logistic path (KKT ~1e-9). `FALSE`
  delegates each per-node fit to the `glmnet` package with the IsingFit
  EBIC path, so the returned `$weights`/`$thresholds` byte-match
  [`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html)
  (to ~1e-16) at the cost of glmnet's looser self-certificate.
  `native = FALSE` needs the optional `glmnet` package (Suggests).
  `min_sum` is supported by both engines; observation `weights` require
  `native = TRUE`.

- labels:

  Optional node labels: one per node, or one per input column (labels of
  dropped constant columns are then discarded).

## Value

A `psychnet` object whose `$weights` is the symmetric weight matrix,
with `$thresholds` (node intercepts), `$rule`, `$nodewise` (the per-node
fits used by
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md))
and `$kkt` (the worst nodewise stationarity residual, read with
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)).

Errors of class `psychnet_bad_input` for non-binary columns, invalid
`gamma`/`nlambda`/`lambda_min_ratio`/`weights`, or mislengthed `labels`.
A nodewise fit that does not reach its optimum raises a
`psychnet_no_converge` warning naming the node.

## References

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

## Examples

``` r
# Ising models need 0/1 data: split the five SRL construct scores at the median.
b <- dichotomize(SRL_GPT)
fit <- ising_fit(b)
fit
#> <psychnet> ising network
#>   nodes: 5   edges: 6   (undirected)
#>   optimality (KKT residual): 1.45e-09
as.data.frame(fit)
#>   from to    weight
#> 1  CSU IV 2.3063761
#> 2  CSU SE 1.8399656
#> 3   IV SE 1.2860533
#> 4  CSU SR 1.7394177
#> 5   IV SR 1.9836011
#> 6   SE SR 0.7570295
certificate(fit)
#>   method  certificate kind certified
#> 1  ising 1.451526e-09  kkt      TRUE
```
