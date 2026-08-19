# EBIC-regularized Gaussian graphical model (graphical lasso)

Selects an L1 penalty by the extended BIC (Foygel & Drton 2010) over a
log-spaced path, then refits the chosen penalty to machine precision so
the fitted precision matrix is the certified global optimum of the
convex objective. Equivalent in purpose to
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
  lambda_path = NULL,
  penalize_diagonal = FALSE,
  refit = TRUE,
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
  Default 0. A positive value is post-estimation processing, so the
  returned graph has `$kkt = NA`; the fitted precision's residual
  remains in `$fit_kkt`.

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

- lambda_path:

  Optional numeric vector of penalties to scan, in the order given. When
  supplied it overrides `nlambda` / `lambda_min_ratio` and is echoed
  back unchanged. Intended for resampling callers: computing the path
  once from the original correlation matrix and reusing it for every
  resample keeps the selected penalties comparable across draws, whereas
  recomputing the path per resample lets the grid drift with each
  sample's largest absolute correlation.

- penalize_diagonal:

  Apply the L1 penalty to the diagonal of the precision matrix as well.
  Default `FALSE` (the working covariance diagonal is held at
  `diag(S)`); `TRUE` holds it at `diag(S) + lambda`, and `$kkt` is
  graded against the matching optimality condition.

- refit:

  What to return at the selected penalty. `TRUE` (default) re-solves it
  to machine precision, so the returned network is the certified
  optimum. `FALSE` returns the path fit exactly as scanned (at the loose
  scan tolerance). This is a bit-compatibility switch, not a speed one:
  the refit is the strictly better answer and costs little, but a
  consumer holding frozen baselines against an unrefitted path scan
  needs the scanned matrix to reproduce them. `"unregularized"` keeps
  the selected support and refits the *unpenalized* graph-restricted MLE
  on it, which is no longer a penalized optimum and is therefore
  certified by
  [`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
  instead.

- labels:

  Optional node labels.

## Value

A `psychnet` object whose `$weights` is the partial-correlation matrix,
with `$precision`, `$lambda`, `$gamma`, `$cor_matrix`, `$ebic`,
`$native`, and `$kkt` (the stationarity residual, or `NA` after
thresholding). Node order in `$weights`, `$nodes`, and `$precision` is
always the column order of the input `data` / `cor_matrix`, never
sorted. With `refit = "unregularized"` the object also carries
`$support`, and `$kkt` is a
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
residual; with `refit = FALSE` the reported `$kkt` reflects the loose
scan tolerance (~1e-4) rather than the refit's ~1e-11.

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
