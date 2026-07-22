# Print a psychnet network

Print a psychnet network

## Usage

``` r
# S3 method for class 'psychnet'
print(x, ...)
```

## Arguments

- x:

  A `psychnet` object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
fit <- pcor_network(SRL_GPT)
print(fit)
#> <psychnet> pcor network
#>   nodes: 5   edges: 10   (undirected)
```
