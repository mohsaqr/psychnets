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
bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
difference_test(bs, type = "strength")
#>    item1 item2    value1    value2     obs_diff       lower       upper p_value
#> 1     V1    V2 0.2647959 0.5762835 -0.311487595 -0.39115437 -0.10270856    0.04
#> 2     V1    V3 0.2647959 0.6051193 -0.340323382 -0.58191788 -0.06545031    0.04
#> 3     V2    V3 0.5762835 0.6051193 -0.028835788 -0.28366276  0.18863071    0.60
#> 4     V1    V4 0.2647959 0.5744056 -0.309609637 -0.46072124 -0.07140963    0.04
#> 5     V2    V4 0.5762835 0.5744056  0.001877958 -0.30046842  0.23289349    0.72
#> 6     V3    V4 0.6051193 0.5744056  0.030713745 -0.28055334  0.27012035    0.96
#> 7     V1    V5 0.2647959 0.2807738 -0.015977936 -0.31186376  0.27056245    0.84
#> 8     V2    V5 0.5762835 0.2807738  0.295509659 -0.10794019  0.47896542    0.16
#> 9     V3    V5 0.6051193 0.2807738  0.324345447  0.05518897  0.49586776    0.04
#> 10    V4    V5 0.5744056 0.2807738  0.293631701  0.03449230  0.55549046    0.00
#>    significant
#> 1         TRUE
#> 2         TRUE
#> 3        FALSE
#> 4         TRUE
#> 5        FALSE
#> 6        FALSE
#> 7        FALSE
#> 8        FALSE
#> 9         TRUE
#> 10        TRUE
```
