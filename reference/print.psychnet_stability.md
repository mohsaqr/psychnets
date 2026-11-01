# Print a centrality-stability result

Print a centrality-stability result

## Usage

``` r
# S3 method for class 'psychnet_stability'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_stability` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the method, the number of subsets per drop
proportion, and the CS-coefficient of each measure.

## Examples

``` r
set.seed(1)
cs <- net_stability(SRL_GPT, method = "pcor", drop_prop = c(0.3, 0.6),
                    iter = 10)
print(cs)
#> <psychnet_stability> pcor, 10 subsets/proportion
#>   CS-coefficient (cor >= 0.70 with 95% certainty):
#>     strength             0.30
#>     expected_influence   0.60
```
