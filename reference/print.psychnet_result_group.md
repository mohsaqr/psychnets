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

`x`, invisibly (for `print`).

For `as.data.frame`, the per-group tables stacked with a `group` column.
