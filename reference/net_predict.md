# Compute how well each node is predicted by the network

`net_predict()` computes the predictability of each node in a network.
The predictability of a node is the proportion of its variance explained
by the nodes it is connected to (Haslbeck & Waldorp, 2018): where
centrality describes how a node is connected, predictability describes
how much of the node's variation its connections account for. It places
the edges on an absolute scale: a node with low predictability depends
largely on variables outside the network, which limits how far it can be
changed through its neighbours.

For Gaussian graphical models, predictability is the explained variance
(R^2), which follows from the estimated precision matrix without the
data. For a relative-importance network
([`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md))
it is the R^2 stored with the fit. For the nodewise models
([`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md))
it is computed from the data by predicting each node from the others:
R^2 for continuous nodes and, for binary nodes, the proportion of
correct classifications (CC) and the normalized accuracy (nCC), the gain
over always predicting the more frequent category as a share of the gain
that was possible.

It returns one row per node, with the type of the node, the metric used
and its value, and the raw accuracy of binary nodes.

[`node_predictability()`](https://pak.dynasite.org/psychnets/reference/node_predictability.md)
returns the same values as a vector for plotting;
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
with `predictability = TRUE` gives bootstrap intervals for Gaussian
graphical models.

## Usage

``` r
net_predict(x, data = NULL, ...)
```

## Arguments

- x:

  A fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)),
  or a
  [psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  (one table per group).

- data:

  The data the network was estimated from; required for the nodewise
  models (`ising`, `ising_sampler`, `mgm`), ignored otherwise.

- ...:

  Unused.

## Value

A data frame with one row per node and columns `node`, `type`
(`"gaussian"` or `"binary"`), `metric` (`"R2"` or `"nCC"`),
`predictability` (the value of `metric`) and `accuracy` (CC for binary
nodes, `NA` for Gaussian ones). For a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md),
one table per group (see
[psychnet_result_group](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md)).

## Details

For a Gaussian graphical model with precision matrix \\\Theta\\ and
correlation matrix \\S\\, the predictability of node \\j\\ is \\R^2_j =
1 - 1 / (\Theta\_{jj} S\_{jj})\\.

For the nodewise models, the columns of `data` are matched to the nodes
by name (the original column names, then the node labels) when the names
allow it, and otherwise by position; binary categories are coded as when
the network was fitted, and a category not seen in fitting is refused.
Incomplete rows are left out.

## Conditions

`psychnet_bad_input` for a network without predictability (for example
from
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md)
or
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md));
for a nodewise model when `data` is missing, lacks columns for some
nodes, contains a category not seen in fitting, or has no complete row;
for a mixed graphical model with a categorical node of more than two
levels; and for a group of within- and between-person networks.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853-861.
[doi:10.3758/s13428-017-0910-x](https://doi.org/10.3758/s13428-017-0910-x)

## See also

[`node_predictability()`](https://pak.dynasite.org/psychnets/reference/node_predictability.md);
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
for centrality;
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
for intervals.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_predict(ebic_glasso(cor_matrix = S, n = 250))
#>   node     type metric predictability accuracy
#> 1   V1 gaussian     R2      0.1568160       NA
#> 2   V2 gaussian     R2      0.2711166       NA
#> 3   V3 gaussian     R2      0.2711166       NA
#> 4   V4 gaussian     R2      0.2711166       NA
#> 5   V5 gaussian     R2      0.2711166       NA
#> 6   V6 gaussian     R2      0.1568160       NA
# Nodewise models need the data they were fitted to.
b <- dichotomize(SRL_GPT)
net_predict(ising_fit(b), data = b)
#>   node   type metric predictability  accuracy
#> 1  CSU binary    nCC      0.7744361 0.9000000
#> 2   IV binary    nCC      0.7769784 0.8966667
#> 3   SE binary    nCC      0.6641791 0.8500000
#> 4   SR binary    nCC      0.7101449 0.8666667
#> 5   TA binary    nCC      0.0000000 0.5000000
```
