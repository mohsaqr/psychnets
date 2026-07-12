# Plot edge-weight case-dropping stability

Draws the four similarity metrics (correlation, mean/median/max absolute
deviation) against the proportion of cases dropped, each with a +/- 1 SD
band; the correlation panel carries the acceptance threshold and
CS-coefficient.

## Usage

``` r
# S3 method for class 'psychnet_casedrop'
plot(x, ...)
```

## Arguments

- x:

  A `psychnet_casedrop` object.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
# Small `iter` / `drop_prop` for a fast example; see [casedrop_reliability()]
# for the defaults a real assessment should use.
plot(casedrop_reliability(SRL_Claude, iter = 5, drop_prop = c(0.25, 0.5)))
```
