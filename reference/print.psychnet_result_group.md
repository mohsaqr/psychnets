# Per-group framework results

The list returned by a framework verb
([`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md))
applied to a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/print.psychnet_group.md):
one result per group level, keyed by level.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) stacks
the per-group tables with a leading `group` column.

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

  Ignored.

- row.names, optional:

  Ignored (S3 consistency).

## Value

`x`, invisibly (for `print`), after printing each group's result under a
group heading.

For `as.data.frame`, the per-group tables stacked with a `group` column.

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
#>           group node  strength expected_influence
#> GPT.1       GPT  CSU 1.1984512         1.19845117
#> GPT.2       GPT   IV 1.0012765         1.00127649
#> GPT.3       GPT   SE 0.8492185         0.84921854
#> GPT.4       GPT   SR 1.2314317         0.53172852
#> GPT.5       GPT   TA 0.6082238        -0.09147935
#> Claude.1 Claude  CSU 1.4794795         1.47947947
#> Claude.2 Claude   IV 0.9920086         0.15982071
#> Claude.3 Claude   SE 1.2652845        -0.37434447
#> Claude.4 Claude   SR 1.9520259         1.09566294
#> Claude.5 Claude   TA 1.3851959        -0.33402780
```
