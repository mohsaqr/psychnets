# Bootstrapped difference test for edges or centralities

Tests, within a single network, whether two edge weights or two node
centralities differ. For every pair it forms the per-resample difference
from the stored bootstrap draws, takes the percentile interval of that
difference, and flags the pair `significant` when the interval excludes
zero; it also reports the two-sided bootstrap p-value (Epskamp, Borsboom
& Fried 2018). This is the within-network counterpart to the edge
accuracy intervals reported by
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md).

## Usage

``` r
difference_test(boot, type = "edge", ci = NULL, p_adjust = "none")

# S3 method for class 'psychnet_difference'
plot(x, style = c("box", "forest"), ...)
```

## Arguments

- boot:

  A `psychnet_bootstrap` object from
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md).

- type:

  Quantity to compare: `"edge"` (default), or any one centrality measure
  bootstrapped by
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  (e.g. `"strength"`, `"expected_influence"`).

- ci:

  Confidence level for the difference interval: a single number in (0,
  1). Defaults to the level used by the bootstrap object.

- p_adjust:

  Multiple-comparison adjustment for the pairwise p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

- x:

  A `psychnet_difference` data frame from `difference_test()`.

- style:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `"box"`
  (default) or `"forest"`.

- ...:

  Unused.

## Value

A tidy data frame, one row per pair, with `item1`, `item2`, the two
observed values, their observed difference, the percentile interval of
the bootstrap difference (`lower`, `upper`), the two-sided `p_value`,
and a logical `significant`; class `psychnet_difference`, so
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it. Edges
are labelled `A--B`, or `A->B` for a directed estimator. An invalid
`boot`, `type` or `ci` is a `psychnet_bad_input` error.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns `x`
invisibly.

## Functions

- `plot(psychnet_difference)`: The bootnet "significance box" matrix:
  items on both axes ordered by observed value, the diagonal showing
  that value, and an off-diagonal cell filled when the pair differs
  (diverging blue when the row item is larger, red when smaller; every
  cell also prints the sign as `+` / `-`, followed by significance
  stars). `style = "forest"` draws one row per pair instead: the
  difference with its interval, a reference line at zero, significant
  pairs as filled points on solid segments and the rest as open points
  on dashed segments.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods*, 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
# method = "pcor" keeps this example fast; the "glasso" default (and
# n_boot >= 1000) is what a real analysis should use.
bs <- net_boot(x, method = "pcor", n_boot = 50, cores = 1)
dt <- difference_test(bs, type = "strength")
dt
#>    item1 item2    value1    value2     obs_diff       lower       upper p_value
#> 1     V1    V2 0.4706062 0.7016822 -0.231076055 -0.48433623  0.06865718    0.12
#> 2     V1    V3 0.4706062 0.7057560 -0.235149802 -0.62747879  0.12286757    0.16
#> 3     V2    V3 0.7016822 0.7057560 -0.004073746 -0.39594719  0.36665932    0.76
#> 4     V1    V4 0.4706062 0.7978913 -0.327285134 -0.52087682 -0.08094752    0.00
#> 5     V2    V4 0.7016822 0.7978913 -0.096209079 -0.40837527  0.24606446    0.44
#> 6     V3    V4 0.7057560 0.7978913 -0.092135333 -0.39241669  0.24270074    0.72
#> 7     V1    V5 0.4706062 0.3862460  0.084360127 -0.29984389  0.31054446    0.84
#> 8     V2    V5 0.7016822 0.3862460  0.315436182 -0.10336275  0.51140973    0.24
#> 9     V3    V5 0.7057560 0.3862460  0.319509929 -0.09902091  0.55579333    0.12
#> 10    V4    V5 0.7978913 0.3862460  0.411645262 -0.02614931  0.57157415    0.08
#>    significant
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4         TRUE
#> 5        FALSE
#> 6        FALSE
#> 7        FALSE
#> 8        FALSE
#> 9        FALSE
#> 10       FALSE
plot(dt)                   # significance-box matrix

plot(dt, style = "forest") # forest plot
```
