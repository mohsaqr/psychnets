# Filter a correlation network with the triangulated maximally filtered graph

`tmfg_network()` filters a correlation matrix to a sparse backbone with
the triangulated maximally filtered graph (TMFG). The TMFG retains, for
\\p\\ variables, exactly \\3(p - 2)\\ edges, chosen greedily to carry as
much of the absolute correlation as possible while the graph stays
planar and built from tetrahedra (four-variable cliques) (Massara et
al., 2016).

Like other information-filtering networks, the TMFG keeps the strongest
associations while imposing a fixed, interpretable structure. It shows
which associations form the core of the correlation structure, without a
sample-size-dependent selection rule. `tmfg_network()` builds the TMFG
from raw data or a correlation matrix.

It returns the graph with the signed correlation of each retained pair
as the edge weight; the edges are marginal correlations, not partial
correlations. The structure that defines the graph is checked by
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md).

For a Gaussian graphical model with partial correlations on the TMFG
graph use
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md);
for an unfiltered correlation network use
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md).
The graph is the one built by `NetworkToolbox::TMFG()`.

## Usage

``` r
tmfg_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given.

- cor_matrix:

  A correlation matrix to use instead of `data`; a covariance matrix is
  rescaled to a correlation matrix.

- n:

  Accepted and ignored: the TMFG needs no sample size. It exists so that
  a call written for the other estimators, `(cor_matrix = , n = )`, also
  works here.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) or
  `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per node of the network.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
whose edge weights are the signed correlations on the retained edges.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the structural check. The object also stores the TMFG graph as a
logical matrix, including any retained edge of weight 0 (`adjacency`),
its four-variable cliques and three-variable separators (`cliques`,
`separators`; used by
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md))
and the correlation matrix (`cor_matrix`). The sample size is recorded
when the network is built from `data` and is `NA` with `cor_matrix`.

## Details

The construction starts from the four variables with the largest sums of
above-average absolute correlations, joined as a tetrahedron. It then
adds the remaining variables one at a time, each time choosing the
variable and the triangular face of the current graph for which the
variable's summed absolute correlation with the face's three variables
is largest, and joining the variable to those three. The result is a
chordal planar graph with \\3(p - 2)\\ edges, which
[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md)
verifies.

The graph always has \\3(p - 2)\\ edges, but a retained pair whose
correlation is exactly zero has edge weight 0 and is therefore absent
from [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
from the printed edge count. With continuous data this happens only for
a correlation matrix with exact zeros, such as a constructed or
block-diagonal `cor_matrix`.

## Conditions

`psychnet_bad_input` when there are fewer than four variables,
`cor_matrix` is not a finite, symmetric, positive semi-definite matrix
with a positive diagonal, an argument is invalid, or `labels` do not
match the nodes; `psychnet_constant_variable` (warning) when a constant
or all-missing column is dropped.

## References

Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network filtering for
big data: Triangulated maximally filtered graph. *Journal of Complex
Networks*, 5(2), 161-178.
[doi:10.1093/comnet/cnw015](https://doi.org/10.1093/comnet/cnw015)

## See also

[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md)
and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
for the check;
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)
for the Gaussian graphical model on this graph;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "tmfg"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 6), 200, 6)
tmfg_network(x)
#> <psychnet> tmfg network
#>   nodes: 6   edges: 12   (undirected)
#>   structural certificate (violated TMFG invariants): 0
```
