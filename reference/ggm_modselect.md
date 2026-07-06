# Stepwise Gaussian graphical model selection (ggmModSelect)

Selects a GGM by extended-BIC model search over edge sets generated from
the glasso path, refitting the *unregularized* maximum-likelihood
precision on each candidate graph, with an optional stepwise add/drop
search. Unlike the graphical lasso, retained edges are not shrunk.
Equivalent in purpose to
[`qgraph::ggmModSelect()`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html),
but pure base R and self-certified via
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md).

## Usage

``` r
ggm_modselect(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  gamma = 0,
  stepwise = TRUE,
  nlambda = 100L,
  lambda_min_ratio = 0.01,
  threshold = 0,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  na_method = c("pairwise", "listwise"),
  native = TRUE,
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

  Sample size (required when `cor_matrix` is supplied).

- gamma:

  EBIC hyperparameter. Default `0`, matching
  [`qgraph::ggmModSelect()`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html):
  because the selected graph is refit with the *unregularized* MLE (no
  edge shrinkage), the extra EBIC penalty is not needed, so gamma 0
  (plain BIC) is the method's intended setting. (The regularized GGMs
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
  and
  [`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)
  keep gamma 0.5.)

- stepwise:

  If `TRUE` (default), refine the best glasso-path graph by a greedy
  single-edge add/drop search.

- nlambda:

  Number of glasso penalties scanned for candidate graphs.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest.

- threshold:

  Partial correlations with absolute value below this are zeroed.
  Default 0.

- cor_method:

  Correlation used when `data` is supplied: `"pearson"` (default),
  `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial, the
  `qgraph`/`bootnet` default for ordinal items). See
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md).

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default)
  or `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- native:

  Solver switch for generating candidate supports: `TRUE` (default) uses
  the pure-R solver; `FALSE` delegates to the `glasso` Fortran package
  (in `Suggests`). The reported precision is the unregularized refit
  either way. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the partial-correlation matrix,
with `$precision`, `$support` (the selected graph), `$gamma`, `$ebic`,
`$cor_matrix`, and `$kkt`.

## Examples

``` r
S <- 0.5^abs(outer(1:6, 1:6, "-"))
ggm_modselect(cor_matrix = S, n = 250)
#> <psychnet> ggm network
#>   nodes: 6   edges: 5   (undirected)
#>   optimality (KKT residual): 2.22e-16
```
