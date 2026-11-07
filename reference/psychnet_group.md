# Per-group psychometric networks

The object returned by
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
when `group` is supplied: a named list of
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
networks, one per level of the grouping variable. It plots as a grid
with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
and is consumed per level by the framework verbs
([`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md),
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)),
each returning a `*_group` result. Data-based resampling and comparison
require independent cross-sectional rows. Groups built from several
occasions per actor need resampling at the actor level and are refused
by these row-based diagnostics.

## Usage

``` r
# S3 method for class 'psychnet_group'
print(x, ...)

# S3 method for class 'psychnet_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'psychnet_group'
summary(object, ...)
```

## Arguments

- x:

  A `psychnet_group` object.

- ...:

  Ignored.

- row.names, optional:

  Ignored (S3 consistency).

- object:

  A `psychnet_group` object.

## Value

`x`, invisibly (for `print`), after printing one row per group with its
node and edge counts and sample size.

For `as.data.frame`, the per-group edge lists stacked with a leading
`group` column: one row per edge per group, columns `group`, `from`,
`to`, `weight`.

For `summary`, a data frame with node/edge counts and mean absolute edge
weight: one row per group, or per group and `level` when a group
contains within/between networks.

## Details

Networks are listed in the grouping variable's level order when it is a
factor (so its first level is group 1 of
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)),
otherwise in order of first appearance. Rows with a missing group value
are dropped with a `psychnet_rows_dropped` warning; an empty-string
level, or a cross-sectional level with fewer than 3 rows, is a
`psychnet_bad_input` error; a level with no more rows than variables
raises a `psychnet_small_sample` warning.

## Examples

``` r
d <- rbind(data.frame(SRL_GPT, source = "GPT"),
           data.frame(SRL_Claude, source = "Claude"))
nets <- psychnet(d, method = "glasso", group = "source")
nets
#> <psychnet_group> 2 networks by source (method: glasso)
#>   group nodes edges   n
#>     GPT     5    10 300
#>  Claude     5    10 300
summary(nets)
#>    group nodes edges mean_abs_weight
#> 1    GPT     5    10       0.2444301
#> 2 Claude     5    10       0.3536997
as.data.frame(nets)
#>     group from to       weight
#> 1     GPT  CSU IV  0.411668665
#> 2     GPT  CSU SE  0.382902233
#> 3     GPT   IV SE  0.159926047
#> 4     GPT  CSU SR  0.354813585
#> 5     GPT   IV SR  0.319094429
#> 6     GPT   SE SR  0.207672075
#> 7     GPT  CSU TA  0.049066683
#> 8     GPT   IV TA  0.110587351
#> 9     GPT   SE TA  0.098718183
#> 10    GPT   SR TA -0.349851570
#> 11 Claude  CSU IV  0.001744267
#> 12 Claude  CSU SE  0.223973564
#> 13 Claude   IV SE -0.402239002
#> 14 Claude  CSU SR  0.728177591
#> 15 Claude   IV SR  0.574170376
#> 16 Claude   SE SR  0.221496430
#> 17 Claude  CSU TA  0.525584044
#> 18 Claude   IV TA -0.013854928
#> 19 Claude   SE TA -0.417575459
#> 20 Claude   SR TA -0.428181457
```
