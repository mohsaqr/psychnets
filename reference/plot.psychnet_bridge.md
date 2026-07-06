# Plot bridge centrality

One sorted horizontal panel per bridge measure (nodes coloured by
community).

## Usage

``` r
# S3 method for class 'psychnet_bridge'
plot(x, measures = NULL, ...)
```

## Arguments

- x:

  A `psychnet_bridge` data frame from
  [`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md).

- measures:

  Which bridge columns to draw. Default: all five.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
S <- 0.3^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 400)
plot(net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2)))
```
