# Plot centrality stability (case-dropping)

Draws the case-dropping stability curves from
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md):
mean rank correlation with the full-sample centrality against the
proportion of cases dropped, one line per measure, with a +/- 1 SD band,
the acceptance threshold, and the CS-coefficient annotated in the
legend.

## Usage

``` r
# S3 method for class 'psychnet_stability'
plot(x, ...)
```

## Arguments

- x:

  A `psychnet_stability` object from
  [`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md).

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
set.seed(1)
d <- matrix(stats::rnorm(200 * 5), 200, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
s <- net_stability(d, drop_prop = c(0.3, 0.6), iter = 10)
plot(s)
```
