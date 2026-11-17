# Networks estimated separately for each group

A `psychnet_group` holds networks estimated separately for each group:
one fitted network per level of the grouping variable, estimated with
the same method and settings. When the same variables are measured in
several groups, such as clinical and control samples, cohorts or
occasions, these networks show how the association structure differs
between them.

It is returned by
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `group`. [`print()`](https://rdrr.io/r/base/print.html) lists the
groups with their numbers of nodes and edges and their sample sizes,
[`summary()`](https://rdrr.io/r/base/summary.html) returns these counts
as a table, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
the edges of all groups;
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the networks in a grid. The analysis functions
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md),
[`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
and
[`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
accept the object and return one result per group (a
[psychnet_result_group](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md));
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
tests the difference between two of its levels.

Resampling and comparison assume independent rows. Groups formed from
several occasions per actor of an event log are refused by these
functions, because they would need resampling at the actor level.

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

  Ignored.

- object:

  A `psychnet_group` object.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly,
after printing one row per group with its numbers of nodes and edges and
its sample size.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with one row per edge per group and columns `group`, `from`,
`to` and `weight`.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per group (or per group and `level` when a group holds
within-actor and between-actor networks) and columns `group`, `nodes`,
`edges` and `mean_abs_weight`.

## Details

Groups are ordered by the levels of the grouping variable when it is a
factor, so its first level is group 1 in
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md),
and otherwise by first appearance. Each network keeps the rows it was
estimated from, so two members can be passed to
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
directly.

## Conditions

`psychnet_rows_dropped` (warning) when rows with a missing group value
are dropped; `psychnet_bad_input` when a group column is missing, a
level is the empty string, there are fewer than 2 levels, or a level of
a numeric table has fewer than 3 rows; `psychnet_small_sample` (warning)
when a level has no more rows than variables.

## See also

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
to estimate the networks;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for a single network;
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
to compare two groups.

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
