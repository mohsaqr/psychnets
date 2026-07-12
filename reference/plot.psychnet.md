# Plot a psychnet network

Renders the estimated network with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
(a Suggested package); `psychnet` objects inherit from `cograph_network`
for exactly this purpose. For the bootstrap / centrality / difference /
stability diagnostics use the dedicated
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
those result objects instead.

## Usage

``` r
# S3 method for class 'psychnet'
plot(x, ...)
```

## Arguments

- x:

  A `psychnet` object.

- ...:

  Passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## Value

The value of
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
invisibly.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 300)
if (requireNamespace("cograph", quietly = TRUE)) {
  plot(fit)
}
```
