# Tidy edge list for a psychnet network

Tidy edge list for a psychnet network

## Usage

``` r
# S3 method for class 'psychnet'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., include_zero = FALSE)
```

## Arguments

- x:

  A `psychnet` object.

- row.names, optional:

  Ignored (for S3 consistency).

- ...:

  Unused.

- include_zero:

  If TRUE, keep zero-weight (absent) edges. Default FALSE.

## Value

A one-row-per-edge `data.frame` with columns `from`, `to`, `weight`.

## Examples

``` r
fit <- pcor_network(SRL_GPT)
as.data.frame(fit)
#>    from to      weight
#> 1   CSU IV  0.41366984
#> 2   CSU SE  0.38511179
#> 3    IV SE  0.14991214
#> 4   CSU SR  0.36093733
#> 5    IV SR  0.32444290
#> 6    SE SR  0.20913531
#> 7   CSU TA  0.06697959
#> 8    IV TA  0.12908925
#> 9    SE TA  0.11123319
#> 10   SR TA -0.39003887
as.data.frame(fit, include_zero = TRUE)
#>    from to      weight
#> 1   CSU IV  0.41366984
#> 2   CSU SE  0.38511179
#> 3    IV SE  0.14991214
#> 4   CSU SR  0.36093733
#> 5    IV SR  0.32444290
#> 6    SE SR  0.20913531
#> 7   CSU TA  0.06697959
#> 8    IV TA  0.12908925
#> 9    SE TA  0.11123319
#> 10   SR TA -0.39003887
```
