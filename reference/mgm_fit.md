# Mixed graphical model

Estimates a mixed graphical model by nodewise L1-penalized regression –
a gaussian (linear) lasso for continuous nodes and a logistic lasso for
binary nodes – with per-node EBIC selection, combined by the AND rule.
Equivalent in purpose to
[`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html), but pure base R
and self-certified: each node's regression reports its stationarity
(KKT) residual (see
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)).

## Usage

``` r
mgm_fit(
  data,
  gamma = 0.25,
  types = NULL,
  nlambda = 100L,
  lambda_min_ratio = 0.01,
  threshold = c("LW", "HW", "none"),
  rule = c("AND", "OR"),
  moderators = NULL,
  weights = NULL,
  na_method = c("pairwise", "listwise"),
  native = TRUE,
  labels = NULL
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations); columns are
  continuous or binary (0/1).

- gamma:

  EBIC hyperparameter. Default 0.25.

- types:

  Optional character vector of node types (`"g"` gaussian, `"c"`
  binary); auto-detected if `NULL`.

- nlambda:

  Number of penalties per nodewise path. Default 100.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest.

- threshold:

  Post-selection coefficient threshold: `"LW"` (default), `"HW"`, or
  `"none"`, matching
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html).

- rule:

  Edge-combination rule: `"AND"` (default) or `"OR"`.

- moderators:

  Optional single column index of a moderator variable. When supplied,
  fits a *moderated* MGM (that variable moderates every pairwise edge)
  and returns a `psychnet_moderated` object to be read with
  [`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md).
  Honours `native` like the unmoderated path: the default base kernel is
  pure R and KKT-certified, and covers gaussian and binary nodes;
  `native = FALSE` uses `glmnet` and additionally allows multi-level
  categorical nodes. `weights` are not supported in this mode.

- weights:

  Optional non-negative observation weights, one per row of the
  (NA-prepared) data. `NULL` (default) is unweighted.

- na_method:

  Missing-data handling: `"pairwise"` (default) single-imputes each
  column over its observed values (mean for continuous, mode for
  binary), keeping the full sample; `"listwise"` drops incomplete rows.

- native:

  Solver switch. `TRUE` (default) uses psychnet's own pure-R,
  dependency-free, self-certified L1 path (KKT ~1e-9). `FALSE` delegates
  each per-node fit to the `glmnet` package with mgm's exact EBIC/LW
  path (gaussian lasso for continuous nodes, 2-class multinomial lasso
  for binary nodes), so the returned edge magnitudes byte-match
  `abs(mgm::mgm()$pairwise$wadj)` (to ~1e-6) at the cost of glmnet's
  looser self-certificate. `native = FALSE` needs the optional `glmnet`
  package (Suggests); `weights` are supported with `native = TRUE` only.

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the symmetric standardized
weight matrix, with `$types` and `$kkt` (the worst nodewise residual). A
binary-binary edge carries the sign of its nodewise-logistic
coefficient; [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html)
reports the same edge as a magnitude only (its sign is undefined for a
categorical-categorical interaction), so compare such edges on
[`abs()`](https://rdrr.io/r/base/MathFun.html). Continuous columns are
standardized internally, binary predictors enter the graph on their 0/1
dummy scale, and binary-response logit coefficients are converted to
`mgm`'s two-class multinomial scale before edge aggregation. With these
conventions the edge magnitudes match
[`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html) closely for
gaussian-gaussian, gaussian-binary, and binary-binary edges alike; weak
edges near the EBIC/threshold boundary can still differ in support
because the penalty is selected on an independent base-R path.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(400)
g1 <- f + stats::rnorm(400); g2 <- f + stats::rnorm(400)
b1 <- (f + stats::rnorm(400) > 0) * 1L
d <- data.frame(g1 = g1, g2 = g2, b1 = b1, n = stats::rnorm(400))
mgm_fit(d)
#> <psychnet> mgm network
#>   nodes: 4   edges: 3   (undirected)
#>   optimality (KKT residual): 6.84e-09
```
