# Return node predictability as a vector for plotting

`node_predictability()` returns node predictability as a numeric vector
in node order, limited to the range 0 to 1, as
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
expects for `pie_values`. Node predictability is the proportion of a
node's variance explained by its neighbours in the network (see
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)),
and is often drawn as a ring around each node of a network plot.

For binary nodes the value is the normalized accuracy (nCC). Use
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
for the full table, including the type of each node and the raw
accuracy.

## Usage

``` r
node_predictability(x, data = NULL)
```

## Arguments

- x:

  A fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

- data:

  The data the network was estimated from; required for the nodewise
  models (`ising`, `ising_sampler`, `mgm`), ignored otherwise.

## Value

A numeric vector with one value per node, in node order, named by node
label, each value between 0 and 1.

## Conditions

The conditions of
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md).

## See also

[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
for the table.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
node_predictability(ebic_glasso(cor_matrix = S, n = 250))
#>        V1        V2        V3        V4        V5        V6 
#> 0.1568160 0.2711166 0.2711166 0.2711166 0.2711166 0.1568160 
```
