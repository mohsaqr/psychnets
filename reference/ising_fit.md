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

  Binary (0/1) data frame or matrix (rows = observations).

- gamma:

  EBIC hyperparameter. Default 0.25.

- rule:

  Edge-combination rule: `"AND"` (default) or `"OR"`.

- nlambda:

  Number of penalties per nodewise path. Default 100.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest.

- min_sum:

  Minimum row sum-score (number of endorsed items); rows below it are
  dropped before fitting. `NULL` (default) keeps every row.

- weights:

  Optional non-negative observation weights, one per retained row.
  `NULL` (default) is unweighted.

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
  `native = FALSE` needs the optional `glmnet` package (Suggests);
  `weights`/`min_sum` are supported with `native = TRUE` only.

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the symmetric weight matrix,
with `$thresholds` (node intercepts) and `$kkt` (the worst nodewise
stationarity residual).

## Examples

``` r
set.seed(1)
z <- matrix(stats::rnorm(400 * 2), 400, 2)
x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(400 * 4), 400)
b <- (x > 0) * 1L
colnames(b) <- paste0("V", 1:4)
ising_fit(b)
#> <psychnet> ising network
#>   nodes: 4   edges: 2   (undirected)
#>   optimality (KKT residual): 3.77e-14
```
