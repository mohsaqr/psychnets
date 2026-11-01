# Plot a Network Comparison Test

Visualises a
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
result. `type = "strength"` (default) and `type = "structure"` draw the
permutation null distribution for the global strength invariance (S) and
the network-structure (M, maximum edge difference) statistics – the
notation of van Borkulo et al. (2022) – with the observed value and
p-value marked; `type = "edges"` draws the observed per-edge absolute
differences, with significantly differing edges both coloured and marked
`*`.

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
a <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
b <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.2^abs(outer(1:4, 1:4, "-")))
colnames(a) <- colnames(b) <- paste0("V", 1:4)
cmp <- net_compare(a, b, iter = 10)  # iter >= 1000 for real use
plot(cmp)                  # global strength permutation null

plot(cmp, type = "edges")  # per-edge differences
```
