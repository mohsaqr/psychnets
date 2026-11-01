# Print a network bootstrap

Print a network bootstrap

## Usage

``` r
# S3 method for class 'psychnet_bootstrap'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_bootstrap` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the method, the number of resamples, the interval
level, and the counts of edges, significant edges and nodes.

## Examples

``` r
set.seed(1)
print(net_boot(SRL_GPT, method = "pcor", n_boot = 50, cores = 1))
#> <psychnet_bootstrap> pcor, 50 resamples, 95% CI
#>   10 edges (9 significant), 5 nodes, measures: strength, expected_influence
```
