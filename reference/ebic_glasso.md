# EBIC-regularized Gaussian graphical model (graphical lasso)

Selects an L1 penalty by the extended BIC (Foygel & Drton 2010) over a
log-spaced path, then refits the chosen penalty to machine precision so
the returned network is the certified global optimum of the convex
objective. Equivalent in purpose to
[`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html)
/ `bootnet`'s `"EBICglasso"` default, but pure base R and self-certified
(see
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)).

## Usage

``` r
ebic_glasso(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  gamma = 0.5,
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

  Optional correlation matrix; if given, `n` is required and `data` is
  ignored.

- n:

  Sample size (required when `cor_matrix` is supplied).

- gamma:

  EBIC hyperparameter. Default 0.5.

- nlambda:

  Number of penalties on the path. Default 100.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest. Default 0.01.

- threshold:

  Partial correlations with absolute value below this are set to zero.
  Default 0.

- cor_method:

  Correlation used when `data` is supplied: `"pearson"` (default),
  `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial for
  ordinal items, the
  [`qgraph::cor_auto`](https://rdrr.io/pkg/qgraph/man/cor_auto.html) /
  `bootnet` default). See
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md).

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default,
  pairwise-complete correlations + nearest-PD projection) or
  `"listwise"` (drop incomplete rows). Identical for complete data.

- native:

  Solver switch. `TRUE` (default) uses psychnet's own pure-R,
  dependency-free, self-certified solver. `FALSE` delegates each
  fixed-penalty solve to the established `glasso` Fortran package (in
  `Suggests`) for speed and byte-identical `glasso`/`qgraph` output, at
  its looser convergence (the reported `$kkt` then shows glasso's
  tolerance rather than ~1e-11).

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the partial-correlation matrix,
with `$precision`, `$lambda`, `$gamma`, `$cor_matrix`, `$ebic`,
`$native`, and `$kkt` (the stationarity residual of the returned
network).

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
fit <- ebic_glasso(cor_matrix = S, n = 250)
fit
#> <psychnet> glasso network
#>   nodes: 6   edges: 5   (undirected)
#>   lambda: 0.004   gamma: 0.5
#>   optimality (KKT residual): 1.10e-10
as.data.frame(fit)
#>   from to    weight
#> 1   V1 V2 0.3681824
#> 2   V2 V3 0.3423189
#> 3   V3 V4 0.3423189
#> 4   V4 V5 0.3423189
#> 5   V5 V6 0.3681824
```
