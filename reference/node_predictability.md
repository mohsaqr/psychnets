# Node predictability as a plotting vector

A thin companion to
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
that returns predictability as a plain numeric vector in node order,
clamped to `[0, 1]` – the form
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
expects for `pie_values` (the predictability ring drawn around each
node). Use
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
when you want the full tidy table.

## Usage

``` r
node_predictability(x, data = NULL)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object.

- data:

  The data the network was estimated from; required for the nodewise
  models (ising / ising_sampler / mgm), ignored for the GGMs.

## Value

A named numeric vector, one value per node (node order), each in
`[0, 1]`.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
node_predictability(ebic_glasso(cor_matrix = S, n = 250))
#>        V1        V2        V3        V4        V5        V6 
#> 0.1568160 0.2711166 0.2711166 0.2711166 0.2711166 0.1568160 
```
