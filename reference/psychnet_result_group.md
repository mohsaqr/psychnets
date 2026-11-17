# Results of an analysis applied to networks per group

A `psychnet_result_group` holds the results of an analysis applied to
networks estimated per group, such as their centralities or bootstrapped
edges: one result per group, keyed by group level, so that they can be
printed together and combined into one table.

It is returned when a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
is passed to
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md),
[`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
or
[`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md).
[`print()`](https://rdrr.io/r/base/print.html) shows each group's result
under a heading, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
the groups' tables with a leading `group` column.

## Usage

``` r
# S3 method for class 'psychnet_result_group'
print(x, ...)

# S3 method for class 'psychnet_result_group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `psychnet_result_group` object.

- ...:

  Ignored by [`print()`](https://rdrr.io/r/base/print.html); passed by
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
  of each group's result (for example `what`).

- row.names, optional:

  Ignored.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly,
after printing each group's result under a group heading.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the groups' tables, each converted by its own
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method,
stacked with a leading `group` column. An error in converting any
group's result is raised, not skipped.

## See also

[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
for the networks; the analysis functions listed above for the columns of
each table.

## Examples

``` r
d <- rbind(data.frame(SRL_GPT, source = "GPT"),
           data.frame(SRL_Claude, source = "Claude"))
nets <- psychnet(d, method = "glasso", group = "source")
cent <- net_centralities(nets)
cent
#> <psychnet_centrality_group> 2 groups: GPT, Claude
#> 
#> --- GPT ---
#>   node  strength expected_influence
#> 1  CSU 1.1984512         1.19845117
#> 2   IV 1.0012765         1.00127649
#> 3   SE 0.8492185         0.84921854
#> 4   SR 1.2314317         0.53172852
#> 5   TA 0.6082238        -0.09147935
#> 
#> --- Claude ---
#>   node  strength expected_influence
#> 1  CSU 1.4794795          1.4794795
#> 2   IV 0.9920086          0.1598207
#> 3   SE 1.2652845         -0.3743445
#> 4   SR 1.9520259          1.0956629
#> 5   TA 1.3851959         -0.3340278
as.data.frame(cent)
#>     group node  strength expected_influence
#> 1     GPT  CSU 1.1984512         1.19845117
#> 2     GPT   IV 1.0012765         1.00127649
#> 3     GPT   SE 0.8492185         0.84921854
#> 4     GPT   SR 1.2314317         0.53172852
#> 5     GPT   TA 0.6082238        -0.09147935
#> 6  Claude  CSU 1.4794795         1.47947947
#> 7  Claude   IV 0.9920086         0.15982071
#> 8  Claude   SE 1.2652845        -0.37434447
#> 9  Claude   SR 1.9520259         1.09566294
#> 10 Claude   TA 1.3851959        -0.33402780
```
