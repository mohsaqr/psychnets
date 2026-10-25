# Network Comparison Test

Permutation test for whether two groups' Gaussian graphical models
differ, on three invariants (van Borkulo et al. 2022): network structure
(`M`, the largest absolute edge difference), global strength (`S`, the
difference in summed absolute edge weights; with `abs = FALSE`, summed
signed weights, i.e. global expected influence) and individual edges
(`E`). Networks are EBIC graphical lassos (clean-room pure R).
Equivalent in purpose to `NetworkComparisonTest::NCT()`, and equal to it
in its permutation draws: under the same seed both permute the same
rows, so the two differ only by the glasso solver's convergence
tolerance.

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
  p_adjust = "none",
  na_method = c("pairwise", "listwise")
)
```

## Arguments

- data1, data2:

  Numeric data frames/matrices with the same variables. When `data1` is
  a
  [psychnet_group](https://pak.dynasite.org/psychnets/reference/print.psychnet_group.md),
  `data2` names the two levels to compare (defaulting to its two levels
  when there are exactly two).

- iter:

  Number of permutations. Default 1000.

- gamma:

  EBIC hyperparameter. Default 0.5.

- paired:

  Logical; within-row swapping for paired designs (row `i` of `data1`
  and of `data2` are the same unit). Default FALSE.

- abs:

  Logical; compare absolute edge weights. Default TRUE.

- weighted:

  Logical; if FALSE, binarize networks first. Default TRUE.

- p_adjust:

  Multiple-comparison adjustment for per-edge p-values (any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method).
  Default `"none"`.

- na_method:

  `"pairwise"` (default) or `"listwise"` handling of missing values; see
  Details.

## Value

An object of class `psychnet_nct`, a list with `$nw1`, `$nw2` (the two
estimated weighted adjacency matrices); the three tests, each with
`observed`, `perm` and `p_value`: `$structure` (statistic M, the largest
absolute edge difference), `$strength` (statistic S, the absolute
difference in global strength; also `by_group`, each network's global
strength) and `$E` (per-edge absolute differences; also `edge_names`, a
`from`/`to` data frame aligned to the per-edge vector); plus `n_iter`,
`n` (rows used per group), `groups` (display names: the argument names,
or the group levels), `paired` and `params`. The letters follow van
Borkulo et al. (2022) and `NetworkComparisonTest`; `$M` and `$S` (which
before psychnets 0.6.0 held the two global tests the other way round)
raise a `psychnet_renamed` error. Printing the object shows every
result; as tables, [`summary()`](https://rdrr.io/r/base/summary.html)
gives the M and S tests and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) one row
per edge (see
[`summary.psychnet_nct()`](https://pak.dynasite.org/psychnets/reference/summary.psychnet_nct.md)).

Errors of class `psychnet_bad_input` for malformed arguments or groups
whose columns differ (in name or order); of class
`psychnet_degenerate_sample` when (pairwise) a variable pair shares
fewer than two observed rows in a group or a permuted split. Warnings of
class `psychnet_constant_variable` (a constant variable in an observed
group) and `psychnet_estimator_mismatch` (a group object built with an
estimator or settings other than EBIC glasso).

## Details

The two groups must hold the same variables in the same order
(non-numeric columns such as an id or a group label are dropped first).

Missing data follow `na_method`. `"pairwise"` (default) keeps every row,
permutes rows as usual, and estimates each (sub)sample from its
pairwise-complete correlation matrix (projected to the nearest positive
definite correlation matrix), with the median pairwise co-observation
count as the sample size – the same treatment as
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).
`"listwise"` drops incomplete rows first (in a paired design, a pair is
dropped when either member is incomplete). Complete data are unaffected
by the choice.

A variable without variation in a group, or in a permuted split (a rare
binary item can be constant in a split), enters that network as an
isolated node – the convention of every psychnets estimator for a
constant column. The statistic is then defined on every split, so the
permutation test stays exact. A constant variable in an observed group
raises a `psychnet_constant_variable` warning.

## References

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods*. [doi:10.1037/met0000476](https://doi.org/10.1037/met0000476)

## Examples

``` r
set.seed(1)
a <- matrix(stats::rnorm(100 * 4), 100, 4)
b <- matrix(stats::rnorm(100 * 4), 100, 4)
colnames(a) <- colnames(b) <- paste0("V", 1:4)
fit <- net_compare(a, b, iter = 10)  # iter >= 1000 for real use
fit
#> Network Comparison Test: a vs b
#> 10 permutations, independent groups, n = 100 and 100
#> 
#>   Test                  Statistic       p
#>   M  network structure      0.000   1.000 
#>   S  global strength        0.000   1.000 
#>        a 0.000 vs b 0.000
#>   E  edges differing       0 of 6
#> 
#>   Edges, largest difference first
#>   Edge          a       b  |diff|       p
#>   V1 - V2   0.000   0.000   0.000   1.000 
#>   V1 - V3   0.000   0.000   0.000   1.000 
#>   V1 - V4   0.000   0.000   0.000   1.000 
#>   V2 - V3   0.000   0.000   0.000   1.000 
#>   V2 - V4   0.000   0.000   0.000   1.000 
#>   V3 - V4   0.000   0.000   0.000   1.000 
#> 
#>   * p < 0.05
summary(fit)
#>                test statistic observed p_value significant
#> 1 network structure         M        0       1       FALSE
#> 2   global strength         S        0       1       FALSE
as.data.frame(fit)
#>   from to weight_1 weight_2 difference p_value significant
#> 1   V1 V2        0        0          0       1       FALSE
#> 2   V1 V3        0        0          0       1       FALSE
#> 3   V2 V3        0        0          0       1       FALSE
#> 4   V1 V4        0        0          0       1       FALSE
#> 5   V2 V4        0        0          0       1       FALSE
#> 6   V3 V4        0        0          0       1       FALSE
```
