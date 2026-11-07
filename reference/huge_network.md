# Nonparanormal graphical model (huge)

Estimates a Gaussian graphical model after a rank-based nonparanormal
transform that relaxes the multivariate-normal assumption, then selects
the L1 penalty by EBIC and refits to the certified optimum. Equivalent
in purpose to `huge::huge()` (nonparanormal) / `bootnet`'s `"huge"`
default, but pure base R and self-certified via
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
on the transformed correlation.

## Usage

``` r
huge_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  npn = c("shrinkage", "truncation", "skeptic"),
  gamma = 0.5,
  nlambda = 100L,
  lambda_min_ratio = 0.01,
  threshold = 0,
  na_method = c("pairwise", "listwise"),
  native = TRUE,
  labels = NULL
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations). Optional if
  `cor_matrix` is supplied (then the transform is skipped).

- cor_matrix:

  Optional pre-transformed correlation matrix; if given, `n` is
  required, `data` and `npn` are ignored.

- n:

  Sample size (required when `cor_matrix` is supplied): a single finite
  number \>= 1.

- npn:

  Nonparanormal transform: `"shrinkage"` (default), `"truncation"`, or
  `"skeptic"` (Spearman).

- gamma:

  EBIC hyperparameter. Default 0.5.

- nlambda:

  Number of penalties on the path. Default 100.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest.

- threshold:

  Partial correlations with absolute value below this are zeroed: a
  single number \>= 0, default 0. With a positive value `$kkt` is `NA`
  for the returned graph and the pre-threshold residual is stored as
  `$fit_kkt`.

- na_method:

  Missing-data handling when `data` is supplied: `"pairwise"` (default,
  with the nonparanormal transform applied per column over observed
  values) or `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- native:

  Solver switch for the glasso path: `TRUE` (default) uses the pure-R
  solver; `FALSE` delegates to the `glasso` Fortran package (in
  `Suggests`). See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per network node.

## Value

A `psychnet` object whose `$weights` is the partial-correlation matrix,
with `$precision`, `$lambda`, `$gamma`, `$cor_matrix` (the transformed
correlation), `$npn`, `$ebic`, and `$kkt`.

## References

Liu, H., Lafferty, J., & Wasserman, L. (2009). The nonparanormal:
Semiparametric estimation of high dimensional undirected graphs.
*Journal of Machine Learning Research*, 10, 2295-2328.

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(300 * 5), 300, 5)
x <- exp(x %*% chol(0.4^abs(outer(1:5, 1:5, "-"))))   # break normality
huge_network(x)
#> <psychnet> huge network
#>   nodes: 5   edges: 5   (undirected)
#>   lambda: 0.06688   gamma: 0.5
#>   optimality (KKT residual): 5.30e-13
```
