# Print a split-half reliability result

Print a split-half reliability result

## Usage

``` r
# S3 method for class 'psychnet_reliability'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_reliability` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the table of between-halves edge metrics.

## Examples

``` r
print(net_split_reliability(SRL_Claude, method = "pcor", iter = 10))
#> # split-half reliability: pcor | 10 iterations (50/50 split)
#>           metric       mean         sd      lower      upper
#> 1   mean_abs_dev 0.07808955 0.02457258 0.03846060 0.10963050
#> 2 median_abs_dev 0.06588545 0.02220565 0.02991325 0.09716984
#> 3    correlation 0.97498013 0.01592382 0.94570362 0.99516211
#> 4    max_abs_dev 0.19047875 0.07735301 0.07527392 0.32270329
```
