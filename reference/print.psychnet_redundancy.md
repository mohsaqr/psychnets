# Print a redundancy result

Print a redundancy result

## Usage

``` r
# S3 method for class 'psychnet_redundancy'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_redundancy` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the flagged pairs, or a note that none were
flagged.

## Examples

``` r
print(redundancy(SRL_Claude))
#> # redundant pairs (proportion < 0.25, r > 0.50): 0 found
#>   none
```
