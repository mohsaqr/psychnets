# Read a panel graphical VAR fit

Methods for the object returned by
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md).

## Usage

``` r
# S3 method for class 'psychnet_panelgvar'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("edges", "parameters", "fit", "diagnostics", "variance"),
  include_zero = FALSE
)

# S3 method for class 'psychnet_panelgvar'
print(x, ...)

# S3 method for class 'psychnet_panelgvar'
summary(object, ...)

# S3 method for class 'psychnet_panelgvar'
plot(x, edge_labels = TRUE, ...)
```

## Arguments

- x, object:

  A `psychnet_panelgvar` object.

- row.names, optional:

  Ignored (S3 consistency).

- ...:

  Passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  by [`plot()`](https://rdrr.io/r/graphics/plot.default.html); named
  arguments override the plotting defaults, with unique argument
  abbreviations accepted. Otherwise ignored.

- what:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
  `"edges"` (default; one row per estimated edge across the three
  networks, with standard error, p-value and 95% Wald interval –
  temporal rows include the autoregressive `from == to` effects),
  `"parameters"` (every model parameter, including means, scaling
  `delta`s, and fixed/pruned entries), `"fit"` (one row of fit indices),
  or `"diagnostics"` (one row: persons `n`, `waves`, the ML `score`
  certificate, `boundary` / `boundary_part`, `se_available`, and
  `starts_agree` of `starts` optimizer starts reaching the best
  saturated optimum), or `"variance"` (one row per variable: the
  model-implied `between`-person variance, the stationary
  `within`-person variance, and `share_between`, their ratio to the
  total – an intraclass correlation).

- include_zero:

  For `what = "edges"`: also return pruned (fixed-zero) edges. Default
  `FALSE`.

- edge_labels:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): print
  each edge's weight on it. Default `TRUE`.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a
data.frame as described under `what`. The edge table has columns
`network`, `from`, `to`, `weight`, `se`, `p`, `ci_lower`, `ci_upper`,
`significant` (p \< the fit's `alpha`); the parameter table has
`matrix`, `row`, `col`, `estimate`, `se`, `p`, `ci_lower`, `ci_upper`,
`fixed`; the diagnostics row has `n`, `waves`, `score`, `boundary`,
`boundary_part`, `se_available`, `starts_agree`, `starts`; the variance
table has `variable`, `between`, `within`, `share_between`.
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) return `x`
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
data.frame with one row per network (`network`, `directed`, `nodes`,
`edges`, `significant`, `mean_abs_weight`).

## Details

By default, [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws the three networks side by side on one shared circle layout (so a
node sits in the same place in every panel), with the temporal
autoregressive effects as self-loops. Line type shows significance:
solid edges have p \< `alpha` (the fit's `alpha`), dashed edges do not
(an edge without a p-value counts as not significant). Colour shows sign
– blue positive, vermillion negative (Okabe-Ito) – and the signed weight
is printed on each edge, so sign never rests on colour alone. A network
with nothing estimated shows "no edges". It needs the `cograph` package.
If `edge_style` is overridden, the legend no longer assigns significance
to solid or dashed lines.

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
summary(fit)
#>           network directed nodes edges significant mean_abs_weight
#> 1        temporal     TRUE     3     6           1      0.07349829
#> 2 contemporaneous    FALSE     3     3           0      0.04156320
#> 3         between    FALSE     3     3           3      0.33740685
as.data.frame(fit, what = "parameters")
#>                matrix row  col      estimate         se             p
#> 1                  mu  X1 <NA>  0.0004526872 0.08934372  9.959573e-01
#> 2                  mu  X2 <NA> -0.0486728484 0.09680339  6.151041e-01
#> 3                  mu  X3 <NA> -0.1261852806 0.09721925  1.943061e-01
#> 4   omega_zeta_within  X2   X1 -0.0171496775 0.03937724  6.631831e-01
#> 5   omega_zeta_within  X3   X1  0.0763290603 0.03855186  4.771409e-02
#> 6   omega_zeta_within  X3   X2 -0.0312108475 0.03866245  4.195137e-01
#> 7   delta_zeta_within  X1   X1  1.0148239963 0.02782219 2.728321e-291
#> 8   delta_zeta_within  X2   X2  0.9507560159 0.02627600 1.109169e-286
#> 9   delta_zeta_within  X3   X3  1.0224638983 0.02775088 3.597988e-297
#> 10               beta  X1   X1  0.4827475939 0.05930052  3.930081e-16
#> 11               beta  X2   X1  0.2352710035 0.04498760  1.697990e-07
#> 12               beta  X3   X1  0.0055641139 0.04885089  9.093171e-01
#> 13               beta  X1   X2 -0.0428914379 0.05061982  3.968140e-01
#> 14               beta  X2   X2  0.3619716930 0.05546525  6.750464e-11
#> 15               beta  X3   X2 -0.0561208438 0.05030812  2.646180e-01
#> 16               beta  X1   X3  0.0257189751 0.04893278  5.991675e-01
#> 17               beta  X2   X3 -0.0754233776 0.04537402  9.646086e-02
#> 18               beta  X3   X3  0.3718081497 0.05975427  4.899623e-10
#> 19 omega_zeta_between  X2   X1  0.4969386490 0.07628681  7.313489e-11
#> 20 omega_zeta_between  X3   X1  0.2523622409 0.09541489  8.171661e-03
#> 21 omega_zeta_between  X3   X2  0.2629196654 0.08703576  2.520884e-03
#> 22 delta_zeta_between  X1   X1  1.0231153290 0.08177032  6.413739e-36
#> 23 delta_zeta_between  X2   X2  1.1729995519 0.08042754  3.522347e-48
#> 24 delta_zeta_between  X3   X3  1.3047968699 0.08623385  1.012794e-51
#>         ci_lower   ci_upper fixed
#> 1  -0.1746577954 0.17556317 FALSE
#> 2  -0.2384040038 0.14105831 FALSE
#> 3  -0.3167315094 0.06436095 FALSE
#> 4  -0.0943276448 0.06002829 FALSE
#> 5   0.0007688039 0.15188932 FALSE
#> 6  -0.1069878542 0.04456616 FALSE
#> 7   0.9602935050 1.06935449 FALSE
#> 8   0.8992560014 1.00225603 FALSE
#> 9   0.9680731767 1.07685462 FALSE
#> 10  0.3665207066 0.59897448 FALSE
#> 11  0.1470969198 0.32344509 FALSE
#> 12 -0.0901818643 0.10131009 FALSE
#> 13 -0.1421044541 0.05632158 FALSE
#> 14  0.2532618093 0.47068158 FALSE
#> 15 -0.1547229413 0.04248125 FALSE
#> 16 -0.0701875069 0.12162546 FALSE
#> 17 -0.1643548288 0.01350807 FALSE
#> 18  0.2546919384 0.48892436 FALSE
#> 19  0.3474192455 0.64645805 FALSE
#> 20  0.0653524861 0.43937200 FALSE
#> 21  0.0923327019 0.43350663 FALSE
#> 22  0.8628484444 1.18338221 FALSE
#> 23  1.0153644650 1.33063464 FALSE
#> 24  1.1357816391 1.47381210 FALSE
```
