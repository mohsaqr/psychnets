# Plot a Network Comparison Test

Visualises a
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
result. `type = "strength"` (default) and `type = "structure"` draw the
permutation null distribution for the global strength invariance (M) and
the maximum edge-difference (S) statistics, with the observed value and
p-value marked; `type = "edges"` draws the observed per-edge absolute
differences, coloured by whether each edge differs significantly.

## Usage

``` r
# S3 method for class 'psychnet_nct'
plot(x, type = c("strength", "structure", "edges"), alpha = 0.05, ...)
```

## Arguments

- x:

  A `psychnet_nct` object from
  [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md).

- type:

  One of `"strength"`, `"structure"`, `"edges"`.

- alpha:

  Significance level for colouring per-edge differences. Default `0.05`.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
set.seed(1)
mk <- function(s) { set.seed(s)
  matrix(stats::rnorm(120 * 4), 120, 4) %*% chol(0.3^abs(outer(1:4, 1:4, "-"))) }
cmp <- net_compare(mk(1), mk(2), iter = 50)
plot(cmp)                  # global strength permutation null

plot(cmp, type = "edges")  # per-edge differences
```
