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

## Examples

``` r
set.seed(1)
# group a has a stronger chain than group b
a <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
b <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.2^abs(outer(1:4, 1:4, "-")))
colnames(a) <- colnames(b) <- paste0("V", 1:4)
fit <- net_compare(a, b, iter = 10)  # iter >= 1000 for real use
print(fit, n = 3)   # the three largest edge differences
#> Network Comparison Test: a vs b
#> 10 permutations, independent groups, n = 100 and 100
#> 
#>   Test                  Statistic       p
#>   M  network structure      0.465   0.091 
#>   S  global strength        1.194   0.091 
#>        a 1.194 vs b 0.000
#>   E  edges differing       0 of 6
#> 
#>   Largest edge differences (3 of 6)
#>   Edge          a       b  |diff|       p
#>   V3 - V4   0.465   0.000   0.465   0.091 
#>   V1 - V2   0.391   0.000   0.391   0.091 
#>   V2 - V3   0.311   0.000   0.311   0.182 
#> 
#>   * p < 0.05
```
