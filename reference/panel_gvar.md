# Panel graphical VAR network (panel-GVAR)

Estimates the three networks of the panel graphical
vector-autoregressive model (Epskamp 2020) from few-wave, many-person
panel data: a directed **temporal** network (within-person lag-1
effects), an undirected **contemporaneous** network (partial
correlations of the within-person innovations), and an undirected
**between**-person network (partial correlations of the stable person
means). The within-person process is a stationary VAR(1) around a random
intercept; means are stationary across waves. This is the model fitted
by
[`psychonetrics::panelgvar()`](https://rdrr.io/pkg/psychonetrics/man/panelvar_family.html),
estimated by maximum likelihood in base R.

## Usage

``` r
panel_gvar(
  data,
  vars,
  id = NULL,
  wave = NULL,
  select = c("none", "prune", "ebic"),
  alpha = 0.01,
  adjust = "none",
  recursive = FALSE,
  stepup = TRUE,
  criterion = c("bic", "none"),
  gamma = 0.5,
  nlambda = 50L,
  lambda_min_ratio = 1e-04,
  refit = FALSE,
  standardize = c("none", "z", "z_per_wave"),
  labels = NULL,
  na_method = "listwise"
)
```

## Arguments

- data:

  A data.frame (or matrix). Either **long** (one row per person per
  wave; give `vars`, `id`, `wave`) or **wide** (one row per person; give
  `vars` as a design matrix).

- vars:

  Long data: character vector of the variable columns. Wide data: a
  character matrix of column names, one row per variable and one column
  per wave (the `psychonetrics` `vars` design); its row names become the
  labels.

- id, wave:

  Long data only: names of the person-identifier and wave columns.
  Person identifiers must be non-missing. Waves are ordered by
  [`sort()`](https://rdrr.io/r/base/sort.html) (or factor level order).

- select:

  Edge selection. `"none"` (default): every edge is estimated and
  reported with its standard error and p-value – significance is shown,
  not imposed (the `psychonetrics` default). `"prune"`: remove edges
  with p \> `alpha` and refit, then step up (see `stepup`) –
  psychonetrics' `prune() |> stepup()`. `"ebic"`: penalized maximum
  likelihood (lasso on every temporal, contemporaneous and between
  parameter) with the penalty chosen by EBIC along a path –
  psychonetrics' `panelgvar(estimator = "PML")`, the panel counterpart
  of
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- alpha:

  Significance level for the `significant` flag, the solid/dashed edges
  of [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  `select = "prune"`. Default `0.01`
  ([`psychonetrics::prune()`](https://rdrr.io/pkg/psychonetrics/man/prune.html)'s
  default).

- adjust:

  Multiplicity correction applied to the edge p-values before pruning,
  one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).
  Default `"none"`
  ([`psychonetrics::prune()`](https://rdrr.io/pkg/psychonetrics/man/prune.html)'s
  default).

- recursive:

  Logical; repeat pruning until no edge has `p > alpha`. Default
  `FALSE`.

- stepup:

  Logical; after pruning, run the
  [`psychonetrics::stepup()`](https://rdrr.io/pkg/psychonetrics/man/stepup.html)
  search: re-add the pruned edge with the largest modification index
  while it is significant at `alpha` and the refit passes a chi-square
  difference test at `alpha` (and lowers BIC, with `criterion = "bic"`).
  Default `TRUE`, because one-pass pruning can remove parameters that
  substitute for each other (e.g. an autoregression and correlated
  cross-lagged effects), each non-significant alone but not jointly.
  `FALSE` reproduces
  [`psychonetrics::prune()`](https://rdrr.io/pkg/psychonetrics/man/prune.html)
  alone. Used only with `select = "prune"`.

- criterion:

  Step-up acceptance criterion besides the chi-square test: `"bic"`
  (default, as
  [`psychonetrics::stepup()`](https://rdrr.io/pkg/psychonetrics/man/stepup.html))
  or `"none"`.

- gamma, nlambda, lambda_min_ratio:

  For `select = "ebic"`: the EBIC hyperparameter (default 0.5), the
  number of penalty values on the path (50; a finite integer \>= 2) and
  the smallest as a fraction of the largest (1e-4) – psychonetrics'
  `find_penalized_lambda()` defaults.

- refit:

  For `select = "ebic"`: refit the selected network by unpenalized
  maximum likelihood (standard errors and p-values then come from that
  refit). Default `FALSE`: the penalized estimates are returned, and
  standard errors, p-values and `significant` are `NA` (they do not
  apply to lasso estimates).

- standardize:

  Rescale the data before fitting, as psychonetrics'
  `panelvar(standardize =)`: `"none"` (default), `"z"` (one mean and SD
  per variable, pooled over waves), or `"z_per_wave"` (z-score every
  variable separately within every wave). `"z_per_wave"` removes
  wave-to-wave drift in means and variances – useful when the
  stationary-means assumption fails, at the price of no longer modelling
  that drift – and makes the edges standardized. Standardization uses
  the available values (SD with denominator n - 1) before listwise
  deletion.

- labels:

  Optional node labels (default: `vars` for long data, the design row
  names for wide data).

- na_method:

  Missing-data handling. Only `"listwise"` is supported: a person
  missing any wave/variable is dropped (with a warning giving the
  count). A pairwise covariance is not a likelihood for this model and
  full- information ML is not yet implemented, so `"pairwise"` is
  refused.

## Value

An object of class `c("psychnet_panelgvar", "netobject_group")`: a named
list of three
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
networks, `temporal` (directed; edge `from -> to` weight
`beta[to, from]`, autoregressive effects kept in `$beta`),
`contemporaneous`, and `between` (undirected partial correlations; the
innovation / between covariance in `$sigma`). Being a `netobject_group`,
it plots as a grid with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
Read it with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
`what = "edges"` (one row per edge with standard error, p-value and 95%
Wald interval), `what = "parameters"` (every model parameter), or
`what = "fit"` (one row of fit indices: chi-square, CFI, TLI, RMSEA,
AIC, BIC), or `what = "diagnostics"` (one row: the score certificate –
largest absolute gradient of the ML fit function on the free parameters,
or the lasso KKT residual for `select = "ebic"` – boundary flag, SE
availability, and how many of the three optimizer starts reached the
best optimum). Fitting, inference and the certificate run on
standardized variables (pooled SD per variable), so they are unit-free
and the certificate does not depend on the variables' units; the model
is scale-equivariant, so estimates and standard errors map back exactly
and are returned on the data scale.

Errors are classed: `psychnet_bad_input` for unusable input (fewer than
three waves, too few persons, non-numeric or missing columns, a
requested `"pairwise"`), `psychnet_no_converge` if the optimizer ends at
an infeasible point. Warnings are classed too: `psychnet_no_converge`
(score above `1e-6` at an interior solution), `psychnet_boundary` (the
final contemporaneous or between covariance is near-singular, i.e.
partial correlations near +/-1 – common with three waves; those
estimates are unreliable, and the score, which need not vanish there, is
quoted in the message), `psychnet_singular_information` (some parameters
are not identified at the returned estimate: their standard errors are
NA, the others exact), `psychnet_approximate_se` (the fit being pruned
had a singular information; the parameters it does not identify were
treated as having infinite standard errors, so pruning removes them),
and `psychnet_listwise_drop` (persons removed for missing values).

## Details

By default the full (saturated) model is returned: every edge with its
standard error, p-value, 95% Wald interval and a `significant` flag at
`alpha`, drawn solid (significant) or dashed (not) by
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

`select = "ebic"` instead minimizes F + lambda \* sum\|theta\| over the
temporal (autoregressive diagonal included), contemporaneous and between
parameters, for 50 log-spaced lambdas from the smallest penalty that
empties every network, keeping the lambda with the lowest EBIC = -2
logL + k log n

- 4 k gamma log V (k free parameters, V wave x variable columns). The
  solver (accelerated proximal gradient, then an active-set Newton
  polish) stays inside the model's domain and is certified by the lasso
  KKT conditions; the selected lambda is re-solved from three starts.
  The likelihood is not convex, so the certificate proves a stationary
  point, and [`print()`](https://rdrr.io/r/base/print.html) reports how
  many starts agree.

With `select = "prune"` the saturated model is fitted, every temporal,
contemporaneous, and between edge whose Wald p-value exceeds `alpha` is
fixed to zero in a single pass, and the model is refitted – the
[`psychonetrics::prune()`](https://rdrr.io/pkg/psychonetrics/man/prune.html)
default workflow (`recursive = TRUE` repeats until every remaining edge
is significant). Then, with `stepup = TRUE` (default), the
[`psychonetrics::stepup()`](https://rdrr.io/pkg/psychonetrics/man/stepup.html)
search adds pruned edges back one at a time while modification indices,
chi-square difference tests and BIC support them, so `select = "prune"`
equals psychonetrics'
`panelgvar() |> runmodel() |> prune() |> stepup()`.

## References

Epskamp, S. (2020). Psychometric network models from time-series and
panel data. *Psychometrika*, 85(1), 206-231.
[doi:10.1007/s11336-020-09697-3](https://doi.org/10.1007/s11336-020-09697-3)

## Examples

``` r
# 300 persons, 3 variables, 4 waves: correlated stable traits plus a
# within-person lag-1 process in which X1 drives X2.
set.seed(1)
n <- 300
beta <- matrix(c(0.4, 0.3, 0, 0, 0.4, 0, 0, 0, 0.4), 3)
trait <- matrix(rnorm(n * 3), n) + rnorm(n)
states <- Reduce(function(w, t) w %*% t(beta) + matrix(rnorm(n * 3), n),
                 seq_len(5), matrix(0, n, 3), accumulate = TRUE)[-(1:2)]
long <- do.call(rbind, Map(function(w, t)
  data.frame(id = seq_len(n), wave = t, trait + w), states, seq_len(4)))
fit <- panel_gvar(long, vars = c("X1", "X2", "X3"), id = "id", wave = "wave")
fit
#> <psychnet_panelgvar> panel graphical VAR (temporal + contemporaneous + between)
#>   persons: 300   waves: 4   nodes: 3
#>   edges (significant/estimated, p < 0.01): temporal 1/6 (+3/3 autoregressive)   contemporaneous 0/3   between 3/3
#>   selection: none (full model)
#>   fit: chisq(66) = 68.78, p = 0.383   CFI = 0.992   TLI = 0.990   RMSEA = 0.012
#>   optimality (ML score): 3.84e-13
#>   starts: 3 of 3 reached the best saturated optimum
as.data.frame(fit)
#>            network from to       weight         se            p      ci_lower
#> 1         temporal   X1 X1  0.482747594 0.05930052 3.930081e-16  0.3665207066
#> 2         temporal   X1 X2  0.235271004 0.04498760 1.697990e-07  0.1470969198
#> 3         temporal   X1 X3  0.005564114 0.04885089 9.093171e-01 -0.0901818643
#> 4         temporal   X2 X1 -0.042891438 0.05061982 3.968140e-01 -0.1421044541
#> 5         temporal   X2 X2  0.361971693 0.05546525 6.750464e-11  0.2532618093
#> 6         temporal   X2 X3 -0.056120844 0.05030812 2.646180e-01 -0.1547229413
#> 7         temporal   X3 X1  0.025718975 0.04893278 5.991675e-01 -0.0701875069
#> 8         temporal   X3 X2 -0.075423378 0.04537402 9.646086e-02 -0.1643548288
#> 9         temporal   X3 X3  0.371808150 0.05975427 4.899623e-10  0.2546919384
#> 10 contemporaneous   X1 X2 -0.017149677 0.03937724 6.631831e-01 -0.0943276448
#> 11 contemporaneous   X1 X3  0.076329060 0.03855186 4.771409e-02  0.0007688039
#> 12 contemporaneous   X2 X3 -0.031210847 0.03866245 4.195137e-01 -0.1069878542
#> 13         between   X1 X2  0.496938649 0.07628681 7.313489e-11  0.3474192455
#> 14         between   X1 X3  0.252362241 0.09541489 8.171661e-03  0.0653524861
#> 15         between   X2 X3  0.262919665 0.08703576 2.520884e-03  0.0923327019
#>      ci_upper significant
#> 1  0.59897448        TRUE
#> 2  0.32344509        TRUE
#> 3  0.10131009       FALSE
#> 4  0.05632158       FALSE
#> 5  0.47068158        TRUE
#> 6  0.04248125       FALSE
#> 7  0.12162546       FALSE
#> 8  0.01350807       FALSE
#> 9  0.48892436        TRUE
#> 10 0.06002829       FALSE
#> 11 0.15188932       FALSE
#> 12 0.04456616       FALSE
#> 13 0.64645805        TRUE
#> 14 0.43937200        TRUE
#> 15 0.43350663        TRUE
as.data.frame(fit, what = "fit")
#>        logl npar df    chisq    pvalue baseline_chisq baseline_df       cfi
#> 1 -5847.958   24 66 68.78458 0.3832688       436.5634          81 0.9921685
#>         tli      rmsea      aic      bic
#> 1 0.9903886 0.01185899 11743.92 11832.81
```
