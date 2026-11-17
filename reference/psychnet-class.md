# Fitted psychometric network

A `psychnet` object holds a fitted psychometric network; every estimator
in the package returns its network in this form. The network consists of
nodes, which are the variables, and weighted edges, which are the
estimated associations: partial correlations for the Gaussian graphical
models, correlations for
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md),
symmetrized regression weights for the Ising and mixed graphical models,
and directed importance shares for
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md).

[`print()`](https://rdrr.io/r/base/print.html) reports the estimator,
the numbers of nodes and edges, the selected penalty and the
certificate; [`summary()`](https://rdrr.io/r/base/summary.html) adds the
range of the edge weights;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the edges as a table; and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
network with the `cograph` package. The certificate is read with
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md).

Analyses of a fitted network, such as
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
and
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md),
return their own result classes, documented on their pages. Networks
estimated per group are collected in a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md).

## Usage

``` r
# S3 method for class 'psychnet'
plot(x, ...)

# S3 method for class 'psychnet'
x$name

# S3 method for class 'psychnet'
print(x, ...)

# S3 method for class 'psychnet'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., include_zero = FALSE)

# S3 method for class 'psychnet'
summary(object, ...)
```

## Arguments

- x, object:

  A `psychnet` object.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), passed
  to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html);
  ignored by the other methods.

- name:

  For `$`: a field name, or one of the earlier names listed in Details.

- row.names, optional:

  Ignored.

- include_zero:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
  keep absent (zero-weight) edges. Default `FALSE`.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with one row per edge (each node pair once for an undirected
network, each ordered pair for a directed one) and columns `from`, `to`
(node labels) and `weight`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) prints the network
with its edge-weight range and returns the edge table invisibly.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
value of
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
invisibly. `$` returns the requested field, or `NULL` when the name is
neither a field nor an earlier name.

## Details

The object is a list that follows the network object of the `cograph`
and `Nestimate` packages, so it can be passed to
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
directly. It stores the weighted adjacency matrix with a zero diagonal
(`weights`), the node table (`nodes`, with `id`, `label` and `name`),
the edge table in the cograph format with integer node indices
(`edges`), whether the network is directed (`directed`), the estimator
(`method`) and the sample size (`n`), together with fields specific to
each estimator, such as the precision matrix (`precision`), the selected
penalty (`lambda`), the EBIC hyperparameter (`gamma`), the correlation
matrix (`cor_matrix`), the certificate (`kkt`), the Ising thresholds
(`thresholds`) and the node types of a mixed model (`types`). Field
access with `$` also recognises the names used by earlier versions:
`graph` (the weights), `labels`, `n_nodes`, `n_edges` and `n_obs`.

## Functions

- `plot(psychnet)`: Draws the network with
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html),
  passing `...` on; needs the suggested package `cograph`.

- `$`: Field access, including the legacy aliases.

- `print(psychnet)`: The estimator, the node and edge counts, the
  selected penalty and the certificate (see
  [`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)).

- `as.data.frame(psychnet)`: The tidy edge list, one row per edge.

- `summary(psychnet)`: The printout plus the edge-weight range.

## Conditions

`psychnet_bad_input` from
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) when the
`cograph` package is not installed.

## See also

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
and the estimators it routes to;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md);
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
for networks per group.

## Examples

``` r
fit <- pcor_network(SRL_GPT)
fit
#> <psychnet> pcor network
#>   nodes: 5   edges: 10   (undirected)
summary(fit)
#> <psychnet> pcor network
#>   nodes: 5   edges: 10   (undirected)
#>   edge weight: range [-0.390, 0.414], mean 0.176
as.data.frame(fit)
#>    from to      weight
#> 1   CSU IV  0.41366984
#> 2   CSU SE  0.38511179
#> 3    IV SE  0.14991214
#> 4   CSU SR  0.36093733
#> 5    IV SR  0.32444290
#> 6    SE SR  0.20913531
#> 7   CSU TA  0.06697959
#> 8    IV TA  0.12908925
#> 9    SE TA  0.11123319
#> 10   SR TA -0.39003887
as.data.frame(fit, include_zero = TRUE)
#>    from to      weight
#> 1   CSU IV  0.41366984
#> 2   CSU SE  0.38511179
#> 3    IV SE  0.14991214
#> 4   CSU SR  0.36093733
#> 5    IV SR  0.32444290
#> 6    SE SR  0.20913531
#> 7   CSU TA  0.06697959
#> 8    IV TA  0.12908925
#> 9    SE TA  0.11123319
#> 10   SR TA -0.39003887
if (requireNamespace("cograph", quietly = TRUE)) plot(fit)
```
