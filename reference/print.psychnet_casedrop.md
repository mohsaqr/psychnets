# Print an edge-weight stability result

Print an edge-weight stability result

## Usage

``` r
# S3 method for class 'psychnet_casedrop'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_casedrop` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the table of edge-similarity metrics by drop
proportion and the edge-weight CS-coefficient.

## Examples

``` r
print(net_casedrop_reliability(SRL_Claude, method = "pcor", iter = 5,
                               drop_prop = c(0.25, 0.5)))
#> # edge-weight stability: pcor | CS = 0.50 (spearman cor >= 0.70 at 95%)
#>           metric drop_prop       mean          sd
#> 1   mean_abs_dev      0.25 0.02262283 0.003969336
#> 2   mean_abs_dev      0.50 0.04285114 0.018510713
#> 3 median_abs_dev      0.25 0.02265162 0.005387633
#> 4 median_abs_dev      0.50 0.03559012 0.014160221
#> 5    correlation      0.25 0.98787879 0.008570991
#> 6    correlation      0.50 0.96848485 0.016262313
#> 7    max_abs_dev      0.25 0.05128306 0.007629198
#> 8    max_abs_dev      0.50 0.10053851 0.033816187
```
