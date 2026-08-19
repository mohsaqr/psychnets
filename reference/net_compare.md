# Network Comparison Test

Permutation test for whether two groups' Gaussian graphical models
differ, on three invariants: global strength (`M`), maximum edge
difference (`S`), and per-edge differences (`E`). Networks are EBIC
graphical lassos (clean-room pure R). Equivalent in purpose to
`NetworkComparisonTest::NCT()`.

## Usage

``` r
net_compare(
  data1,
  data2 = NULL,
  iter = 1000L,
  gamma = 0.5,
  paired = FALSE,
  abs = TRUE,
  weighted = TRUE,
  p_adjust = "none"
)
```

## Arguments

- data1, data2:

  Numeric data frames/matrices with the same columns.

- iter:

  Number of permutations. Default 1000.

- gamma:

  EBIC hyperparameter. Default 0.5.

- paired:

  Logical; within-row swapping for paired designs. Default FALSE.

- abs:

  Logical; compare absolute edge weights. Default TRUE.

- weighted:

  Logical; if FALSE, binarize networks first. Default TRUE.

- p_adjust:

  Multiple-comparison adjustment for per-edge p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

## Value

An object of class `psychnet_nct` with `$nw1`, `$nw2`, and `$M`, `$S`,
`$E` (each `observed`, `perm`, `p_value`); `$E` also carries
`edge_names`, a `from`/`to` data frame aligned to the per-edge vector.

## Examples

``` r
set.seed(1)
a <- matrix(stats::rnorm(150 * 5), 150, 5)
b <- matrix(stats::rnorm(150 * 5), 150, 5)
colnames(a) <- colnames(b) <- paste0("V", 1:5)
fit <- net_compare(a, b, iter = 25)
fit
#> Network Comparison Test (25 permutations)
#>   Global strength (M): observed 0.000, p = 1.000
#>   Network structure (S): observed 0.000, p = 1.000
```
