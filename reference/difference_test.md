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
```

## Arguments

- boot:

  A `psychnet_bootstrap` object from
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md).

- type:

  Quantity to compare: `"edge"` (default), or any centrality measure
  bootstrapped by
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  (e.g. `"strength"`, `"expected_influence"`).

- ci:

  Confidence level for the difference interval. Defaults to the level
  used by the bootstrap object.

- p_adjust:

  Multiple-comparison adjustment for the pairwise p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

## Value

A tidy data frame, one row per pair, with `item1`, `item2`, the two
observed values, their observed difference, the percentile interval of
the bootstrap difference (`lower`, `upper`), the two-sided `p_value`,
and a logical `significant`.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
bs <- net_boot(x, n_boot = 100)
difference_test(bs, type = "strength")
#>    item1 item2    value1    value2     obs_diff       lower       upper p_value
#> 1     V1    V2 0.2647959 0.5762835 -0.311487595 -0.41317357 -0.05972764    0.02
#> 2     V1    V3 0.2647959 0.6051193 -0.340323382 -0.56020707  0.03436308    0.06
#> 3     V2    V3 0.5762835 0.6051193 -0.028835788 -0.27130909  0.24372879    0.92
#> 4     V1    V4 0.2647959 0.5744056 -0.309609637 -0.46080253 -0.07327107    0.04
#> 5     V2    V4 0.5762835 0.5744056  0.001877958 -0.32302967  0.23716978    0.92
#> 6     V3    V4 0.6051193 0.5744056  0.030713745 -0.29444015  0.25212538    0.98
#> 7     V1    V5 0.2647959 0.2807738 -0.015977936 -0.30937799  0.26101233    0.90
#> 8     V2    V5 0.5762835 0.2807738  0.295509659 -0.13797478  0.52488790    0.20
#> 9     V3    V5 0.6051193 0.2807738  0.324345447 -0.06206998  0.46596541    0.10
#> 10    V4    V5 0.5744056 0.2807738  0.293631701  0.00744765  0.51300860    0.02
#>    significant
#> 1         TRUE
#> 2        FALSE
#> 3        FALSE
#> 4         TRUE
#> 5        FALSE
#> 6        FALSE
#> 7        FALSE
#> 8        FALSE
#> 9        FALSE
#> 10        TRUE
```
