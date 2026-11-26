# Estimate temporal, contemporaneous and between-person networks from panel data

`panel_gvar()` estimates temporal, contemporaneous and between-person
networks from panel data, which follow many persons over a few
measurement waves. It fits the panel graphical vector-autoregressive
(panel-GVAR) model, which separates three kinds of association among the
measured variables (Epskamp, 2020).

The **temporal** network is directed: an edge from A to B is the
within-person effect of A at one wave on B at the next, given all other
variables, and a self-loop is an autoregressive effect. The
**contemporaneous** network holds the partial correlations among a
person's deviations at the same wave that the temporal effects leave
unexplained. The **between** network holds the partial correlations
among the persons' stable means. Together they distinguish variables
that predict one another over time within persons, that vary together
within persons at the same wave, and that differ together between
persons.

The model is fitted by maximum likelihood to long data (one row per
person and wave) or wide data (one row per person). By default every
edge is estimated and reported with its standard error, p-value, 95%
Wald interval and a `significant` flag at `alpha`. `select = "prune"`
removes the edges that are not significant and refits the model, then
adds back removed edges that improve the fit (`stepup`).
`select = "ebic"` estimates the three networks with a lasso penalty
chosen by the extended Bayesian information criterion (EBIC), the panel
counterpart of
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).
The fitted model is checked against its optimality conditions (see
Details), and `as.data.frame(x, what = "diagnostics")` reports the
result.

The result is a `psychnet_panelgvar` object holding the three networks
([psychnet_panelgvar](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)).

The model needs at least three waves. For a network from cross-sectional
data, use
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md);
for networks of within-person deviations without a temporal model,
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `within`. The model is the one fitted by
[`psychonetrics::panelgvar()`](https://rdrr.io/pkg/psychonetrics/man/panelvar_family.html).

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

  A data frame or matrix: **long** (one row per person per wave; give
  `vars`, `id` and `wave`) or **wide** (one row per person; give `vars`
  as a design matrix).

- vars:

  Long data: the names of the variable columns. Wide data: a character
  matrix of column names with one row per variable and one column per
  wave (the `psychonetrics` design); its row names become the labels.

- id, wave:

  Long data only: the names of the person and wave columns. Neither may
  contain missing values or be one of `vars`, and each person has at
  most one row per wave. Waves are ordered by factor level order, by
  value for numbers, dates and times, numerically for character values
  that are all numbers, and alphabetically (with a warning) for any
  other character values; pass a factor to set the order.

- select:

  Edge selection: `"none"` (default) estimates every edge and reports
  its significance; `"prune"` removes non-significant edges and refits
  (see Details); `"ebic"` estimates penalized networks with the penalty
  chosen by EBIC.

- alpha:

  Significance level, in (0, 1), for the `significant` flag, the solid
  and dashed edges of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and pruning.
  Default `0.01`, as in
  [`psychonetrics::prune()`](https://rdrr.io/pkg/psychonetrics/man/prune.html).

- adjust:

  Multiplicity correction for the p-values tested in pruning, one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).
  Default `"none"`. Used only with `select = "prune"`; the reported
  p-values and `significant` flags are always unadjusted.

- recursive:

  With `select = "prune"`: repeat pruning until no edge has p \>
  `alpha`. Default `FALSE`.

- stepup:

  With `select = "prune"`: after pruning, add back removed edges that
  improve the fit (see Details). Default `TRUE`, because one-pass
  pruning can remove parameters that are each non-significant but
  jointly significant, such as an autoregression and correlated
  cross-lagged effects. `FALSE` gives pruning alone.

- criterion:

  With `stepup`: `"bic"` (default) also requires BIC to decrease when an
  edge is added; `"none"` uses the significance tests only.

- gamma, nlambda, lambda_min_ratio:

  With `select = "ebic"`: the EBIC hyperparameter (default `0.5`), the
  number of penalty values (default `50`, a whole number \>= 2), and the
  smallest penalty as a fraction of the largest (default `1e-4`).

- refit:

  With `select = "ebic"`: refit the selected networks by unpenalized
  maximum likelihood, so that standard errors and p-values are
  available. Default `FALSE` returns the penalized estimates, for which
  standard errors, p-values and `significant` are `NA`.

- standardize:

  Rescaling before fitting: `"none"` (default), `"z"` (one mean and
  standard deviation per variable, pooled over waves), or `"z_per_wave"`
  (each variable z-scored within each wave). `"z_per_wave"` removes
  changes in means and variances across waves, which the model would
  otherwise assume absent, and makes the edges standardized.
  Standardization uses the available values before listwise deletion.

- labels:

  Optional node labels, unique and non-empty (default: `vars` for long
  data, the design row names for wide data).

- na_method:

  Only `"listwise"` (default): a person with any missing value is
  dropped. `"pairwise"` is refused, because a pairwise covariance is not
  a likelihood for this model.

## Value

An object of class `c("psychnet_panelgvar", "netobject_group")`
([psychnet_panelgvar](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)):
a named list of three networks of class `psychnet`, `temporal`
(directed), `contemporaneous` and `between` (undirected).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per edge with its standard error, p-value and Wald interval;
`as.data.frame(x, what = ...)` gives the `"parameters"`, `"fit"`
indices, `"variance"` decomposition, `"covariance"` matrices, EBIC
`"path"` or `"diagnostics"`, whose `score` column is the certificate.
[`summary()`](https://rdrr.io/r/base/summary.html) gives one row per
network and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the three networks.

## Details

Write \\y\_{it}\\ for the scores of person \\i\\ at wave \\t\\. The
model is \\y\_{it} = \mu + b_i + w\_{it}\\, with means \\\mu\\ shared by
all waves, stable person deviations \\b_i\\ with covariance
\\\Sigma_b\\, and a stationary within-person process \\w\_{it} = B
w\_{i,t-1} + \zeta\_{it}\\ with innovation covariance \\\Sigma\_\zeta\\.
The temporal network is \\B\\ (the edge from \\j\\ to \\i\\ is
\\B\_{ij}\\); the contemporaneous and between networks are the partial
correlations implied by \\\Sigma\_\zeta\\ and \\\Sigma_b\\. When the
means drift across waves, `standardize = "z_per_wave"` removes the
drift.

Pruning fixes to zero, in one pass, every temporal, contemporaneous and
between edge whose Wald p-value (after `adjust`) exceeds `alpha`, and
refits; `recursive = TRUE` repeats until every remaining edge is
significant. The step-up search then frees the removed edge with the
largest modification index while that index is significant at `alpha`,
and keeps it when the likelihood-ratio test of the refit is significant
at `alpha` and, with `criterion = "bic"`, BIC decreases; the first
rejected edge ends the search. This is psychonetrics' `prune()` followed
by `stepup()`.

With `select = "ebic"`, the lasso penalty applies to the temporal
(autoregressive effects included), contemporaneous and between
parameters on the data scale, over `nlambda` log-spaced values from the
smallest penalty that empties every network down to `lambda_min_ratio`
times that value. The selected penalty minimizes \\\mathrm{EBIC} = -2
\log L + k \log n + 4 \gamma k \log V\\, with \\k\\ the number of
nonzero free parameters, \\n\\ the number of persons and \\V\\ the
number of wave-by-variable columns. Penalized estimates below `1e-5` in
absolute value on the data scale are set to zero before the certificate
is computed.

The certificate is the largest absolute gradient of the likelihood over
the free parameters or, for `select = "ebic"` without `refit`, the
largest violation of the lasso optimality (Karush–Kuhn–Tucker)
conditions. Pruned and penalized models are refitted or re-solved after
selection, so the certificate refers to the networks returned. It is
computed on standardized variables and does not depend on the variables'
units. The likelihood is not convex: a value near zero shows a
stationary point, and the fit is started from several points whose
agreement is reported.

## Conditions

`psychnet_bad_input` when `data`, `vars`, `id` or `wave` cannot be
resolved (columns missing or duplicated, `id` and `wave` the same or
among `vars`, missing persons or waves, a person with two rows in one
wave), a variable is not numeric, there are fewer than two variables or
fewer than three waves, there are no more complete persons than
wave-by-variable columns, the stacked covariance is singular, a variable
does not vary as `standardize` requires, `labels` is invalid, an
argument is outside its allowed values, or `na_method` is not
`"listwise"`; `psychnet_no_converge` (error) when no start reaches a
feasible optimum, the fitted model is worse than its nested baseline, or
no penalized model is feasible. Warnings: `psychnet_no_converge` when
the score exceeds `1e-6` at an interior solution; `psychnet_boundary`
when the contemporaneous or between covariance is near singular (partial
correlations near +/-1, common with three waves), so estimates involving
it are unreliable; `psychnet_singular_information` when some parameters
are not identified at the estimate (their standard errors are `NA`);
`psychnet_approximate_se` when a fit being pruned had a singular
information and its unidentified parameters were pruned;
`psychnet_listwise_drop` when persons are removed for missing values;
`psychnet_wave_order` when character waves are ordered alphabetically;
`psychnet_ignored_argument` when `adjust` is set without pruning.

## References

Epskamp, S. (2020). Psychometric network models from time-series and
panel data. *Psychometrika*, 85(1), 206-231.
[doi:10.1007/s11336-020-09697-3](https://doi.org/10.1007/s11336-020-09697-3)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

## See also

[psychnet_panelgvar](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)
for reading the fit;
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
for cross-sectional data;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `within`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
for the certificates of the other estimators.

## Examples

``` r
# srl_panel: 100 students x 4 weeks, simulated from a panel-GVAR in which
# effort in one week drives efficacy in the next.
fit <- panel_gvar(srl_panel,
                  vars = c("efficacy", "planning", "effort", "enjoyment",
                           "anxiety"),
                  id = "student", wave = "week")
fit
#> <psychnet_panelgvar> panel graphical VAR (temporal + contemporaneous + between)
#>   persons: 100   waves: 4   nodes: 5
#>   edges (significant/estimated, p < 0.01): temporal 1/20 (+4/5 autoregressive)   contemporaneous 5/10   between 2/10
#>   selection: none (full model)
#>   fit: chisq(170) = 207.19, p = 0.0273   CFI = 0.957   TLI = 0.946   RMSEA = 0.047
#>   optimality (ML score): 5.91e-13
#>   starts: 3 of 3 reached the best saturated optimum
as.data.frame(fit)
#>            network      from        to       weight         se            p
#> 1         temporal   anxiety   anxiety  0.296628278 0.10027278 3.094278e-03
#> 2         temporal   anxiety  efficacy  0.026009263 0.08031053 7.460450e-01
#> 3         temporal   anxiety    effort -0.014936851 0.08549642 8.613096e-01
#> 4         temporal   anxiety enjoyment  0.036480078 0.08218780 6.571420e-01
#> 5         temporal   anxiety  planning  0.003375992 0.07857850 9.657308e-01
#> 6         temporal  efficacy   anxiety -0.122338773 0.13119129 3.510668e-01
#> 7         temporal  efficacy  efficacy -0.096672635 0.13467200 4.728574e-01
#> 8         temporal  efficacy    effort -0.185033790 0.13554948 1.722328e-01
#> 9         temporal  efficacy enjoyment -0.161261754 0.13020900 2.155367e-01
#> 10        temporal  efficacy  planning -0.127621701 0.12388997 3.029530e-01
#> 11        temporal    effort   anxiety  0.087931626 0.12119415 4.681187e-01
#> 12        temporal    effort  efficacy  0.352722389 0.11519665 2.199257e-03
#> 13        temporal    effort    effort  0.483841150 0.13198874 2.465897e-04
#> 14        temporal    effort enjoyment  0.196809315 0.11886939 9.778725e-02
#> 15        temporal    effort  planning  0.168996084 0.11367355 1.370996e-01
#> 16        temporal enjoyment   anxiety  0.036854988 0.10677727 7.299759e-01
#> 17        temporal enjoyment  efficacy  0.125084640 0.10252259 2.224387e-01
#> 18        temporal enjoyment    effort  0.035363362 0.10913195 7.459061e-01
#> 19        temporal enjoyment enjoyment  0.325408489 0.12000699 6.696350e-03
#> 20        temporal enjoyment  planning -0.015083147 0.10027960 8.804404e-01
#> 21        temporal  planning   anxiety  0.024730130 0.15150380 8.703365e-01
#> 22        temporal  planning  efficacy  0.173937878 0.14529357 2.312490e-01
#> 23        temporal  planning    effort  0.084407115 0.15504123 5.861544e-01
#> 24        temporal  planning enjoyment  0.114456276 0.14927613 4.432349e-01
#> 25        temporal  planning  planning  0.396614715 0.15188718 9.021313e-03
#> 26 contemporaneous  efficacy   anxiety  0.028205118 0.06800858 6.783404e-01
#> 27 contemporaneous  efficacy    effort  0.342448405 0.06025581 1.321643e-08
#> 28 contemporaneous  efficacy enjoyment -0.038177031 0.06832147 5.763083e-01
#> 29 contemporaneous  efficacy  planning  0.591106125 0.04431346 1.369209e-40
#> 30 contemporaneous    effort   anxiety  0.002383863 0.06667728 9.714799e-01
#> 31 contemporaneous    effort enjoyment  0.317013534 0.06007371 1.312612e-07
#> 32 contemporaneous enjoyment   anxiety  0.051935419 0.06629439 4.333888e-01
#> 33 contemporaneous  planning   anxiety -0.062045823 0.06668056 3.521157e-01
#> 34 contemporaneous  planning    effort  0.188697279 0.06478954 3.585880e-03
#> 35 contemporaneous  planning enjoyment  0.255668116 0.06276252 4.629603e-05
#> 36         between  efficacy   anxiety -0.148011468 0.12809365 2.478881e-01
#> 37         between  efficacy    effort  0.069437587 0.13041944 5.944369e-01
#> 38         between  efficacy enjoyment  0.083371182 0.14572743 5.672518e-01
#> 39         between  efficacy  planning  0.313900573 0.11991975 8.855349e-03
#> 40         between    effort   anxiety -0.023216660 0.13180138 8.601770e-01
#> 41         between    effort enjoyment  0.344615578 0.12950101 7.788502e-03
#> 42         between enjoyment   anxiety  0.033435880 0.14994682 8.235472e-01
#> 43         between  planning   anxiety  0.196830672 0.12719722 1.217560e-01
#> 44         between  planning    effort  0.302598713 0.11877742 1.084632e-02
#> 45         between  planning enjoyment  0.316230825 0.13044624 1.534117e-02
#>       ci_lower   ci_upper significant
#> 1   0.10009724 0.49315931        TRUE
#> 2  -0.13139648 0.18341501       FALSE
#> 3  -0.18250676 0.15263306       FALSE
#> 4  -0.12460505 0.19756521       FALSE
#> 5  -0.15063504 0.15738703       FALSE
#> 6  -0.37946898 0.13479144       FALSE
#> 7  -0.36062490 0.16727963       FALSE
#> 8  -0.45070588 0.08063830       FALSE
#> 9  -0.41646671 0.09394320       FALSE
#> 10 -0.37044157 0.11519817       FALSE
#> 11 -0.14960455 0.32546780       FALSE
#> 12  0.12694110 0.57850368        TRUE
#> 13  0.22514798 0.74253432        TRUE
#> 14 -0.03617041 0.42978904       FALSE
#> 15 -0.05379998 0.39179215       FALSE
#> 16 -0.17242462 0.24613460       FALSE
#> 17 -0.07585594 0.32602522       FALSE
#> 18 -0.17853134 0.24925806       FALSE
#> 19  0.09019911 0.56061787        TRUE
#> 20 -0.21162754 0.18146125       FALSE
#> 21 -0.27221186 0.32167212       FALSE
#> 22 -0.11083229 0.45870805       FALSE
#> 23 -0.21946811 0.38828234       FALSE
#> 24 -0.17811957 0.40703212       FALSE
#> 25  0.09892130 0.69430813        TRUE
#> 26 -0.10508924 0.16149948       FALSE
#> 27  0.22434919 0.46054762        TRUE
#> 28 -0.17208466 0.09573060       FALSE
#> 29  0.50425334 0.67795891        TRUE
#> 30 -0.12830121 0.13306894       FALSE
#> 31  0.19927123 0.43475583        TRUE
#> 32 -0.07799920 0.18187004       FALSE
#> 33 -0.19273733 0.06864568       FALSE
#> 34  0.06171212 0.31568244        TRUE
#> 35  0.13265584 0.37868039        TRUE
#> 36 -0.39907040 0.10304747       FALSE
#> 37 -0.18617981 0.32505498       FALSE
#> 38 -0.20224933 0.36899169       FALSE
#> 39  0.07886218 0.54893897        TRUE
#> 40 -0.28154261 0.23510929       FALSE
#> 41  0.09079827 0.59843289        TRUE
#> 42 -0.26045448 0.32732624       FALSE
#> 43 -0.05247130 0.44613265       FALSE
#> 44  0.06979925 0.53539818       FALSE
#> 45  0.06056089 0.57190076       FALSE
as.data.frame(fit, what = "fit")
#>        logl npar  df    chisq     pvalue baseline_chisq baseline_df       cfi
#> 1 -6802.528   60 170 207.1872 0.02733013       1084.025         215 0.9572082
#>         tli      rmsea      aic      bic
#> 1 0.9458809 0.04677052 13725.06 13881.37
as.data.frame(fit, what = "variance")
#>    variable  between   within share_between
#> 1  efficacy 65.39817 66.42409     0.4961087
#> 2  planning 72.92447 57.68286     0.5583490
#> 3    effort 82.63922 69.22678     0.5441588
#> 4 enjoyment 53.83746 65.44361     0.4513496
#> 5   anxiety 95.00513 59.25030     0.6158949
```
