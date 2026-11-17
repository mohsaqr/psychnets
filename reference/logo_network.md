# Estimate a sparse Gaussian graphical model on a filtered graph (LoGo)

`logo_network()` estimates a sparse Gaussian graphical model on a
filtered graph with the Local-Global (LoGo) method. LoGo works in two
steps: it first chooses the network structure by information filtering
with the triangulated maximally filtered graph (TMFG), and then
estimates the partial correlations as the maximum-likelihood Gaussian
graphical model on that graph (Barfuss et al., 2016).

The TMFG retains the strongest correlations under the constraint that
the graph stays planar and is built from tetrahedra (four-variable
cliques), which fixes the number of edges at \\3(p - 2)\\ for \\p\\
variables. Because the graph is decomposable, the maximum-likelihood
estimate has a closed form, assembled from the inverses of small local
correlation matrices, so no penalty or iterative fitting is needed. The
model is estimated from raw data or from a correlation matrix with its
sample size.

`logo_network()` returns the partial-correlation network of this model,
the maximum-likelihood network on the TMFG graph, which
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
checks.

The graph alone, with correlations as edge weights, is returned by
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md).
For a network whose sparsity is chosen by the data use
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
or
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md).
The model is the one fitted by `NetworkToolbox::LoGo()`.

## Usage

``` r
logo_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  threshold = 0,
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
  rescaled to a correlation matrix. `n` is then required.

- n:

  Sample size, required with `cor_matrix` and recorded on the result: a
  single finite number \>= 1. Ignored when `data` is given.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- threshold:

  Absolute partial correlation below which an edge is set to zero after
  estimation: a single number \>= 0. Default `0`. A logical value is
  refused. A positive value leaves the certificate missing.

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) or
  `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per node of the network.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
whose edge weights are partial correlations.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the check. The object also stores the precision matrix
(`precision`), the TMFG graph as a logical matrix (`support`) and the
correlation matrix (`cor_matrix`).

## Details

With \\C\\ the four-variable cliques and \\T\\ the three-variable
separators of the TMFG, the precision matrix is \$\$\Theta = \sum\_{c
\in C} \[S\_{cc}^{-1}\]^0 - \sum\_{t \in T} \[S\_{tt}^{-1}\]^0,\$\$
where \\\[\cdot\]^0\\ pads a block with zeros to the full size. The
certificate is computed by
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
with the TMFG as the graph. After a positive `threshold` it is missing;
[`print()`](https://rdrr.io/r/base/print.html) still shows the value for
the network before thresholding.

The closed form needs every clique and separator block to be invertible.
A block that is singular or nearly so (reciprocal condition number below
`1e-10`, as with duplicate or collinear variables) is refused rather
than approximated, as is a precision matrix whose certificate exceeds
`1e-6`.

## Conditions

`psychnet_singular` (error, also of class `psychnet_bad_input`) when a
clique or separator block is singular or nearly so, or the precision
matrix fails its certificate; `psychnet_bad_input` when there are fewer
than four variables, an argument is out of range, `threshold` is
logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
matrix with a positive diagonal, `n` is missing with `cor_matrix`, or
`labels` do not match the nodes; `psychnet_constant_variable` (warning)
when a constant or all-missing column is dropped.

## References

Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
Parsimonious modeling with information filtering networks. *Physical
Review E*, 94(6), 062306.
[doi:10.1103/PhysRevE.94.062306](https://doi.org/10.1103/PhysRevE.94.062306)

Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network filtering for
big data: Triangulated maximally filtered graph. *Journal of Complex
Networks*, 5(2), 161-178.
[doi:10.1093/comnet/cnw015](https://doi.org/10.1093/comnet/cnw015)

## See also

[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md)
for the filtered graph;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
for the check;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "logo"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(300 * 6), 300, 6)
logo_network(x)
#> <psychnet> logo network
#>   nodes: 6   edges: 12   (undirected)
#>   optimality (KKT residual): 4.44e-16
```
