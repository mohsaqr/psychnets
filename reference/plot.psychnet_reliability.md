# Plot split-half reliability

Histograms of the four between-halves edge metrics across split-half
iterations, with each observed mean marked.

## Usage

``` r
# S3 method for class 'psychnet_reliability'
plot(x, ...)
```

## Arguments

- x:

  A `psychnet_reliability` object.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
# \donttest{
plot(network_reliability(SRL_Claude))

# }
```
