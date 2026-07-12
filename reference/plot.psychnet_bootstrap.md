# Plot a network bootstrap

Visualises a
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
result. `type = "edges"` (default) draws the bootstrapped edge-weight
confidence intervals sorted by the observed weight (bootnet's
edge-accuracy plot); `type = "centrality"` draws the bootstrapped
centrality intervals, one sorted panel per measure; `type = "edge_diff"`
and `type = "centrality_diff"` draw the bootstrapped difference
"significance box" matrix for edges or for one centrality;
`type = "predictability"` draws the node predictability intervals (only
when
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
was run with `predictability = TRUE`).

## Usage

``` r
# S3 method for class 'psychnet_bootstrap'
plot(
  x,
  type = c("edges", "centrality", "edge_diff", "centrality_diff", "predictability"),
  measure = NULL,
  ...
)
```

## Arguments

- x:

  A `psychnet_bootstrap` object from
  [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md).

- type:

  One of `"edges"`, `"centrality"`, `"edge_diff"`, `"centrality_diff"`,
  `"predictability"`.

- measure:

  For `"centrality"`/`"centrality_diff"`, which measure(s) to draw.
  Default: all bootstrapped measures (`"centrality"`) or the first
  (`"centrality_diff"`).

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plot it draws.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
bs <- net_boot(x, n_boot = 50, cores = 1)   # n_boot >= 1000 for real use
plot(bs)                       # edge-weight CIs

plot(bs, type = "centrality")  # centrality CIs
```
