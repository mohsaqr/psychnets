# Plot node centralities

Draws the centrality table returned by
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).
`type = "bar"` gives one sorted horizontal lollipop panel per measure;
`type = "line"` gives the qgraph/bootnet centrality plot — one faceted
panel per measure, nodes on a shared vertical axis, a line-and-marker
series within each panel.

## Usage

``` r
# S3 method for class 'psychnet_centrality'
plot(
  x,
  type = c("bar", "line"),
  scale = c("raw", "z", "relative"),
  measures = NULL,
  ...
)
```

## Arguments

- x:

  A `psychnet_centrality` data frame from
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).

- type:

  `"bar"` (default) for one sorted lollipop panel per measure, or
  `"line"` for the faceted qgraph-style centrality plot.

- scale:

  For `type = "line"`, the per-measure transform: `"raw"` (default —
  each faceted panel keeps its own axis, so no rescaling is needed),
  `"z"` (z-score per measure, centred at 0), or `"relative"` (each
  measure min-max scaled to \[0, 1\]). `type = "bar"` always shows raw
  values.

- measures:

  Which measure columns to draw. Default: all of them.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 300)
plot(net_centralities(fit))

plot(net_centralities(fit), type = "line")
```
