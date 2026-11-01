# Tidy a network bootstrap

Tidy a network bootstrap

## Usage

``` r
# S3 method for class 'psychnet_bootstrap'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., significant = FALSE)
```

## Arguments

- x:

  A `psychnet_bootstrap` object.

- row.names, optional:

  Ignored (S3 consistency).

- ...:

  Unused.

- significant:

  If `TRUE`, return only the edges whose confidence interval excludes
  zero. Default `FALSE` (all edges).

## Value

A `data.frame` with one row per edge and columns `from`, `to`,
`observed`, `mean` (bootstrap mean), `lower`, `upper` (percentile
interval), `prop_nonzero` (share of resamples in which the edge is
nonzero) and `significant` (interval excludes zero).

## Examples

``` r
set.seed(1)
bs <- net_boot(SRL_GPT, method = "pcor", n_boot = 50, cores = 1)
as.data.frame(bs)
#>    from to    observed        mean        lower      upper prop_nonzero
#> 1   CSU IV  0.41366984  0.42821646  0.320910107  0.5155981            1
#> 2   CSU SE  0.38511179  0.37843279  0.290250707  0.4657695            1
#> 3    IV SE  0.14991214  0.13991828  0.012979431  0.2477617            1
#> 4   CSU SR  0.36093733  0.35485385  0.249544420  0.4314031            1
#> 5    IV SR  0.32444290  0.31510972  0.225618123  0.3890399            1
#> 6    SE SR  0.20913531  0.22106026  0.116065498  0.3120899            1
#> 7   CSU TA  0.06697959  0.06117521 -0.034105733  0.1912335            1
#> 8    IV TA  0.12908925  0.13204858  0.016655264  0.2346835            1
#> 9    SE TA  0.11123319  0.11487828  0.006240649  0.2333458            1
#> 10   SR TA -0.39003887 -0.39130789 -0.480946467 -0.2961468            1
#>    significant
#> 1         TRUE
#> 2         TRUE
#> 3         TRUE
#> 4         TRUE
#> 5         TRUE
#> 6         TRUE
#> 7        FALSE
#> 8         TRUE
#> 9         TRUE
#> 10        TRUE
as.data.frame(bs, significant = TRUE)
#>    from to   observed       mean        lower      upper prop_nonzero
#> 1   CSU IV  0.4136698  0.4282165  0.320910107  0.5155981            1
#> 2   CSU SE  0.3851118  0.3784328  0.290250707  0.4657695            1
#> 3    IV SE  0.1499121  0.1399183  0.012979431  0.2477617            1
#> 4   CSU SR  0.3609373  0.3548538  0.249544420  0.4314031            1
#> 5    IV SR  0.3244429  0.3151097  0.225618123  0.3890399            1
#> 6    SE SR  0.2091353  0.2210603  0.116065498  0.3120899            1
#> 8    IV TA  0.1290892  0.1320486  0.016655264  0.2346835            1
#> 9    SE TA  0.1112332  0.1148783  0.006240649  0.2333458            1
#> 10   SR TA -0.3900389 -0.3913079 -0.480946467 -0.2961468            1
#>    significant
#> 1         TRUE
#> 2         TRUE
#> 3         TRUE
#> 4         TRUE
#> 5         TRUE
#> 6         TRUE
#> 8         TRUE
#> 9         TRUE
#> 10        TRUE
```
