# Plot a bootstrapped difference test

Draws the bootnet-style "significance box" matrix for the pairwise
difference test returned by
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md):
items on both axes ordered by their observed value, the diagonal showing
the observed value, and each off-diagonal cell filled when that pair
differs significantly (red when the row item is larger, blue when
smaller). `style = "forest"` instead draws a forest plot: one row per
pair, the bootstrapped difference as a point with its confidence
interval, a reference line at zero, and significant pairs (interval
excluding zero) emphasised.

## Usage

``` r
# S3 method for class 'psychnet_difference'
plot(x, style = c("box", "forest"), ...)
```

## Arguments

- x:

  A `psychnet_difference` data frame from
  [`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md).

- style:

  `"box"` (default) for the significance-box matrix, or `"forest"` for a
  forest plot of the pairwise differences with their CIs.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
plot(difference_test(bs, type = "strength"))                   # box matrix

plot(difference_test(bs, type = "strength"), style = "forest") # forest plot
```
