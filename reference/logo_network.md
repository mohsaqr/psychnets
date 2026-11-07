# Local-Global sparse inverse covariance (LoGo)

Estimates a sparse Gaussian graphical model whose
conditional-independence structure is the chordal TMFG: the precision is
the closed-form Gaussian Markov random field on that graph (Barfuss et
al. 2016). Equivalent in purpose to `NetworkToolbox::LoGo()` /
`bootnet`'s `"LoGo"` default, pure base R and self-certified via
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
(the precision reproduces `S` exactly on the TMFG support).

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

  Numeric data frame or matrix (rows = observations). Optional if
  `cor_matrix` is supplied.

- cor_matrix:

  Optional correlation matrix; if given, `n` is required.

- n:

  Sample size (recorded on the result; required with `cor_matrix`): a
  single finite number \>= 1.

- cor_method:

  Correlation when `data` is supplied: `"pearson"` (default),
  `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- threshold:

  Partial correlations with absolute value below this are zeroed: a
  single number \>= 0, default 0. With a positive value `$kkt` is `NA`
  for the returned graph and the pre-threshold residual is stored as
  `$fit_kkt`.

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default)
  or `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per network node.

## Value

A `psychnet` object whose `$weights` is the partial-correlation matrix,
with `$precision`, `$support` (the TMFG graph), `$cor_matrix`, and
`$kkt`. A clique or separator block that is singular or nearly so
(reciprocal condition number below 1e-10, e.g. a duplicate column), or a
precision whose certificate exceeds 1e-6, is refused with an error of
class `psychnet_singular` rather than returned as an uncertified
network.

## References

Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
Parsimonious modeling with information filtering networks. *Physical
Review E*, 94(6), 062306.
[doi:10.1103/PhysRevE.94.062306](https://doi.org/10.1103/PhysRevE.94.062306)

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(300 * 6), 300, 6)
logo_network(x)
#> <psychnet> logo network
#>   nodes: 6   edges: 12   (undirected)
#>   optimality (KKT residual): 4.44e-16
```
