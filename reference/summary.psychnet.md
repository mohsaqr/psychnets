# Summarize a psychnet network

Summarize a psychnet network

## Usage

``` r
# S3 method for class 'psychnet'
summary(object, ...)
```

## Arguments

- object:

  A `psychnet` object.

- ...:

  Unused.

## Value

The tidy edge list (invisibly); prints a summary as a side effect.

## Examples

``` r
fit <- pcor_network(SRL_GPT)
summary(fit)
#> <psychnet> pcor network
#>   nodes: 5   edges: 10   (undirected)
#>   edge weight: range [-0.390, 0.414], mean 0.176
```
