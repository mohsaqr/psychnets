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

`x`, invisibly (for `print`).

For `as.data.frame`, the per-group edge lists stacked with a `group`
column.

For `summary`, a data frame with node/edge counts and mean absolute edge
weight: one row per group, or per group and `level` when a group
contains within/between networks.
