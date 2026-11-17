# Estimate a mixed graphical model for continuous and binary variables

`mgm_fit()` estimates a mixed graphical model: one network for variables
of different types, such as continuous symptom ratings and binary
indicators. An edge connects two variables that remain associated after
conditioning on all other variables, and an absent edge indicates
conditional independence (Haslbeck & Waldorp, 2020), so the network
shows which variables are directly related when they cannot all be
treated as Gaussian or all as binary.

The function regresses each variable on all others with the L1-penalized
(lasso) regression suited to its type, linear for continuous and
logistic for binary variables, and selects each variable's penalty by
the extended Bayesian information criterion (EBIC). Small coefficients
are then removed by a threshold (`threshold`, default `"LW"`), and the
two estimates of each pair are combined by the AND rule (default) or the
OR rule into an undirected edge whose weight is their mean. Variable
types are detected from the data or given in `types`. Each regression is
checked against its optimality (Karush–Kuhn–Tucker) conditions;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the largest violation. A threshold that removes a coefficient
changes the network after fitting, so its certificate is then reported
as missing; `threshold = "none"` keeps it. With `moderators`, the
function instead fits a moderated model in which one variable can change
the strength of every edge; see
[`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md).

The result is an undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)),
or a `psychnet_moderated` fit with `moderators`.

For binary data only, use
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md);
for continuous data only,
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).
Categorical variables with more than two levels need `native = FALSE`.
The model is the one fitted by
[`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html).

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

  A data frame or matrix with one row per observation. Columns are
  continuous, binary, or (with `native = FALSE`) categorical with more
  than two levels. With `native = TRUE` every column must be numeric.
  Columns with fewer than two distinct observed values are dropped. A
  character column with more than 10 distinct values is refused as an
  identifier; convert it with
  [`factor()`](https://rdrr.io/r/base/factor.html) if it is a
  categorical variable.

- gamma:

  EBIC hyperparameter, a single number \>= 0; larger values give sparser
  networks. Default `0.25`.

- types:

  Node types, one per column of `data`: `"g"` (Gaussian) or `"c"`
  (categorical). Default `NULL` detects them as
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) does: factor and
  character columns are categorical, and so are numeric columns with at
  most 10 distinct integer values. Such a numeric column with more than
  two values is named in a warning, because it changes the model for
  that node; Likert items often meet this rule. A non-factor column
  declared `"c"` with more than 10 distinct values is refused.

- nlambda:

  Number of penalty values on each nodewise path, a whole number \>= 2.
  Default `100`.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest, in (0, 1). Default
  `0.01`. With `native = FALSE`, leaving it unset uses glmnet's own
  default, as [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) does.

- threshold:

  Threshold on the selected coefficients: `"LW"` (default), `"HW"` or
  `"none"` (see Details), as in
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html).

- rule:

  `"AND"` (default) keeps an edge when both nodewise estimates are
  nonzero; `"OR"` when either is.

- moderators:

  Optional moderator: one column of `data`, by name or by position. The
  moderated model lets this variable change the strength of every
  pairwise edge and returns a `psychnet_moderated` fit, read with
  [`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md)
  and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
  Needs at least three variables and does not accept `weights`. Default
  `NULL`.

- weights:

  Optional non-negative observation weights, one per row after
  missing-data handling. Only their relative sizes matter, as in
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html). Requires
  `native = TRUE`. Default `NULL` (unweighted).

- na_method:

  `"pairwise"` (default) fills each missing value with the variable's
  mean (continuous) or most frequent level (categorical) and keeps every
  row; `"listwise"` drops rows with a missing value. The two agree for
  complete data.

- native:

  `TRUE` (default) uses the package's own solver, for continuous and
  binary variables. `FALSE` fits each regression with the `glmnet`
  package (which must be installed), following
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html); it also allows
  categorical variables with more than two levels, but not `weights`,
  and its certificate reflects glmnet's convergence tolerance. A network
  with a multi-level categorical node cannot be passed to
  [`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md).
  One-hot encoding such a variable is not equivalent: it turns one node
  into several.

- labels:

  Optional node labels: one per node, or one per input column (labels of
  dropped constant columns are then discarded).

## Value

An undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md));
with `moderators`, a `psychnet_moderated` fit (see
[`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md)).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per edge (`from`, `to`, `weight`),
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
the largest optimality violation of the nodewise regressions (missing
when the threshold removed a coefficient), and
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
the predictability of each node. The object also stores the node types
(`types`), the number of levels per node (`levels`, 1 for continuous)
and the certificate of the regressions before thresholding (`fit_kkt`).

## Details

Continuous variables are standardized before fitting, binary predictors
enter as 0/1 indicators, and the coefficient of a binary response is
expressed on the two-class multinomial scale of
[`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) (half the logistic
coefficient). Edge magnitudes are therefore comparable with those of
[`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html). A binary-binary
edge keeps the sign of its logistic coefficient, whereas
[`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) reports such an
edge without a sign; compare those edges on their absolute value. Weak
edges near the selection or threshold boundary can differ between the
two engines.

With \\\beta_i\\ the selected coefficients of variable \\i\\, \\q\\ the
number of its predictor columns (\\p - 1\\ for \\p\\ continuous or
binary variables) and \\n\\ the sample size, the `"LW"` threshold
removes coefficients smaller in absolute value than \\\\\beta_i\\\_2
\sqrt{\log(q) / n}\\, and the `"HW"` threshold those smaller than
\\\sqrt{\log(q) / n}\\.

Node order is the column order of `data`, and a categorical variable
with several levels remains one node.

## Conditions

`psychnet_bad_input` when `data` is not a data frame or matrix, `types`
is not one `"g"`/`"c"` per column, a character column looks like an
identifier, a non-factor column declared `"c"` has more than 10 values,
a categorical variable has more than two levels or a column is not
numeric with `native = TRUE`, `moderators` does not name one column or
names a dropped constant column, a moderated model has fewer than three
variables or is given `weights`, `weights` is combined with
`native = FALSE`, `native = FALSE` is requested without the `glmnet`
package, or `gamma`, `nlambda`, `lambda_min_ratio`, `weights` or
`labels` is invalid; `psychnet_type_promotion` (warning) when a numeric
column with more than two values is detected as categorical;
`psychnet_no_converge` (warning) when a nodewise regression does not
reach its optimum; `psychnet_degenerate_sample` (warning) when a binary
variable takes a single value on every row with positive weight. The
last two are raised with `native = TRUE`.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
time-varying mixed graphical models in high-dimensional data. *Journal
of Statistical Software*, 93(8), 1-46.
[doi:10.18637/jss.v093.i08](https://doi.org/10.18637/jss.v093.i08)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

## See also

[`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md)
for the moderated model;
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
for binary data;
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
for continuous data;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
for the optimality check;
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md);
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `method = "mgm"`.

## Examples

``` r
set.seed(1)
f <- stats::rnorm(400)
g1 <- f + stats::rnorm(400); g2 <- f + stats::rnorm(400)
b1 <- (f + stats::rnorm(400) > 0) * 1L
d <- data.frame(g1 = g1, g2 = g2, b1 = b1, n = stats::rnorm(400))
fit <- mgm_fit(d)
fit
#> <psychnet> mgm network
#>   nodes: 4   edges: 3   (undirected)
#>   optimality (KKT residual): 1.01e-10
as.data.frame(fit)
#>   from to    weight
#> 1   g1 g2 0.3421962
#> 2   g1 b1 0.3288350
#> 3   g2 b1 0.3589604
```
