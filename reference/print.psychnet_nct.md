# Print a Network Comparison Test

Shows every result: the two global tests (M, network structure; S,
global strength with each group's value), the number of differing edges,
and every edge – largest difference first – with each group's edge
weight.

## Usage

``` r
# S3 method for class 'psychnet_nct'
print(x, n = Inf, alpha = 0.05, ...)
```

## Arguments

- x:

  A `psychnet_nct` object.

- n:

  Number of edges to list, largest difference first. Default `Inf`
  (all).

- alpha:

  Significance level for the `*` marks. Default 0.05.

- ...:

  Unused.

## Value

`x`, invisibly.
