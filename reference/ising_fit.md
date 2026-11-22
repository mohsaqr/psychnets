# Estimate a regularized Ising network for binary items

`ising_fit()` estimates a regularized Ising network from binary items,
such as symptoms that are present or absent or items that are endorsed
or not. The Ising model describes the items through a threshold for each
item and an interaction for each pair, the association between two items
after conditioning on all other items, and is estimated here by the
nodewise lasso method of van Borkulo et al. (2014).

A positive edge means that the two items tend to be endorsed together, a
negative edge that endorsing one makes the other less likely, and an
absent edge that the two are conditionally independent; the network
shows which items are directly associated once the rest are taken into
account. The function regresses each item on all other items by
L1-penalized (lasso) logistic regression and selects each item's penalty
by the extended Bayesian information criterion (EBIC). Every pair is
estimated twice, once in each item's regression. Under the AND rule
(default) an edge is kept when both estimates are nonzero, under the OR
rule when either is; its weight is the mean of the two. Rows can be
filtered by their sum score (`min_sum`) and weighted (`weights`). Each
regression is checked against its optimality (Karush–Kuhn–Tucker)
conditions, and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the largest violation across the items.

The result is an undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

For unpenalized estimates with significance-based pruning, use
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md);
for data that combine continuous and binary variables,
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md).
Ordinal or continuous items can be split with
[`dichotomize()`](https://pak.dynasite.org/psychnets/reference/dichotomize.md)
first. The model is the one fitted by
[`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html),
which `native = FALSE` reproduces.

## Usage

``` r
ising_fit(
  data,
  gamma = 0.25,
  rule = c("AND", "OR"),
  nlambda = 100L,
  lambda_min_ratio = NULL,
  min_sum = NULL,
  weights = NULL,
  na_method = c("pairwise", "listwise"),
  native = TRUE,
  labels = NULL
)
```

## Arguments

- data:

  Binary data: a data frame or matrix with one row per observation.
  Columns may be numeric 0/1, logical, or two-level factors (the first
  level is coded 0). Any other column, such as a numeric column coded
  1/2 or a factor with more than two levels, is an error naming it. A
  constant column is dropped.

- gamma:

  EBIC hyperparameter, a single number \>= 0; larger values give sparser
  networks. Default `0.25`.

- rule:

  `"AND"` (default) keeps an edge when both nodewise estimates are
  nonzero; `"OR"` when either is.

- nlambda:

  Number of penalty values on each nodewise path, a whole number \>= 2.
  Default `100`.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest, in (0, 1). Default
  `NULL` uses glmnet's default, as
  [`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html)
  does: `1e-4`, or `0.01` when there are fewer observations than
  predictors. Both engines select from glmnet's path (see Details).

- min_sum:

  Minimum number of endorsed items per row; rows below it are dropped
  before fitting. Default `NULL` keeps every row.

- weights:

  Optional non-negative observation weights, one per row after
  missing-data handling (and filtered with `min_sum`). Only their
  relative sizes matter: `weights = w` and `weights = 7 * w` give the
  same network. Requires `native = TRUE`. Default `NULL` (unweighted).

- na_method:

  `"pairwise"` (default) fills each missing value with the item's more
  frequent value and keeps every row; `"listwise"` drops rows with a
  missing value. The two agree for complete data.

- native:

  `TRUE` (default) uses the package's own solver. `FALSE` fits each
  regression with the `glmnet` package (which must be installed) and
  reproduces
  [`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html);
  its certificate reflects glmnet's convergence tolerance and is usually
  larger. `weights` are not accepted with `native = FALSE`.

- labels:

  Optional node labels: one per node, or one per input column (labels of
  dropped constant columns are then discarded).

## Value

An undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per edge (`from`, `to`, `weight`),
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
the largest optimality violation of the nodewise regressions, and
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
the predictability of each node. The object also stores the node
thresholds (`thresholds`), the combination rule (`rule`) and the
nodewise regressions (`nodewise`).

## Details

Each regression's penalty is chosen by the EBIC from glmnet's default
path (Friedman, Hastie & Tibshirani 2010), as in
[`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html):
`nlambda` penalties spaced evenly on the log scale from the smallest
penalty that sets every coefficient to zero down to `lambda_min_ratio`
times it, stopped early by glmnet's rule (from the fifth penalty on,
when the explained deviance grows by less than 1e-5, relative for a
continuous response, or exceeds 0.999). The native solver fits each
penalty to its optimality conditions, so its weights agree with
[`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html)
up to glmnet's convergence tolerance.

Edge weights and thresholds are on the log-odds scale of the 0/1 coding.
The weight of an edge is the change in the log-odds of one item when the
other is endorsed, averaged over the two regressions; an item's
threshold is its log-odds when every other item is 0.

## Conditions

`psychnet_bad_input` when a column is not binary, fewer than two usable
columns remain, `gamma`, `nlambda`, `lambda_min_ratio`, `min_sum`,
`weights`, `native` or `labels` is invalid, `min_sum` removes every row
(or every row with positive weight), `weights` is combined with
`native = FALSE`, or `native = FALSE` is requested without the `glmnet`
package; `psychnet_constant_variable` (warning) when a constant column
is dropped; `psychnet_no_converge` (warning) when a nodewise regression
does not reach its optimum; `psychnet_degenerate_sample` (warning) when
an item takes a single value on every retained row, so its regression
has no finite intercept. The last two are raised with `native = TRUE`.

## References

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths
for generalized linear models via coordinate descent. *Journal of
Statistical Software*, 33(1), 1-22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

## See also

[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
for the unregularized model;
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
for mixed continuous and binary data;
[`dichotomize()`](https://pak.dynasite.org/psychnets/reference/dichotomize.md)
to prepare items;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
for the optimality check;
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md);
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `method = "ising"`.

## Examples

``` r
# Ising models need 0/1 data: split the five SRL construct scores at the median.
b <- dichotomize(SRL_GPT)
fit <- ising_fit(b)
fit
#> <psychnet> ising network
#>   nodes: 5   edges: 6   (undirected)
#>   optimality (KKT residual): 8.65e-10
as.data.frame(fit)
#>   from to    weight
#> 1  CSU IV 2.3226325
#> 2  CSU SE 1.8563260
#> 3   IV SE 1.2819940
#> 4  CSU SR 1.7682981
#> 5   IV SR 1.9940481
#> 6   SE SR 0.7660948
certificate(fit)
#>   method  certificate kind certified
#> 1  ising 8.650803e-10  kkt      TRUE
```
