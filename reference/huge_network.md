# Estimate a nonparanormal Gaussian graphical model

`huge_network()` estimates a nonparanormal Gaussian graphical model: a
regularized partial-correlation network for variables that are normal
only after an unknown monotone transformation. The nonparanormal model
(Liu et al., 2009) relaxes the assumption of the Gaussian graphical
model that the variables are jointly normal, an assumption that skewed
or heavy-tailed questionnaire and symptom data often violate, so the
network of partial correlations is defined on the transformed scale and
is unaffected by the shape of each variable's distribution.

The transformation is estimated from the ranks of the observations.
`huge_network()` transforms the data, computes the correlation matrix of
the transformed variables and estimates a regularized
partial-correlation network on it with the EBIC graphical lasso, as in
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).
`npn` chooses the transformation: normal scores of the ranks with
shrinkage (default) or truncation, or the "skeptic", a correlation
derived from Spearman's rank correlation.

The returned network is the optimum of the penalized likelihood at the
selected penalty, which
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
checks.

For approximately normal data,
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
fits the same model without the transformation; its
`cor_method = "auto"` treats ordinal items through polychoric
correlations instead. The model is the one fitted by `huge::huge()` with
the nonparanormal transformation.

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

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given.

- cor_matrix:

  A correlation matrix of already transformed data, used instead of
  `data`; no transformation is applied and `n` is required.

- n:

  Sample size, required with `cor_matrix`: a single finite number= 1.
  Ignored when data is given.

- npn:

  The transformation: `"shrinkage"` (default), `"truncation"` or
  `"skeptic"`. Ignored with `cor_matrix`.

- gamma:

  EBIC hyperparameter, a number \>= 0; larger values favour sparser
  networks. Default `0.5`.

- nlambda:

  Number of penalties scanned, a whole number \>= 2. Default `100`.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest, in (0, 1). Default
  `0.01`.

- threshold:

  Absolute partial correlation below which an edge is set to zero after
  estimation: a single number \>= 0. Default `0`. A logical value is
  refused. A positive value leaves the certificate missing.

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default)
  transforms each variable over its observed values and uses
  pairwise-complete correlations; `"listwise"` drops incomplete rows
  first.

- native:

  `TRUE` (default) uses the package's own graphical lasso solver;
  `FALSE` uses the `glasso` package, which must be installed. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per node of the network.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
whose edge weights are partial correlations on the transformed scale.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the check. The object also stores the precision matrix
(`precision`), the selected penalty (`lambda`), `gamma`, the transformed
correlation matrix (`cor_matrix`), the transformation used (`npn`,
`"precomputed"` with `cor_matrix`) and the EBIC of the selected network
(`ebic`).

## Details

The shrinkage transform replaces each observation by \\\Phi^{-1}(r /
(n + 1))\\, with \\r\\ its rank among the \\n\\ observed values of the
variable. The truncation transform uses \\\Phi^{-1}(r / n)\\ truncated
to \\\[\delta, 1 - \delta\]\\, with \\\delta = 1 / (4 n^{1/4} \sqrt{\pi
\log n})\\. The skeptic uses \\2 \sin(\pi \rho_S / 6)\\, with \\\rho_S\\
Spearman's correlation. A correlation matrix that is not positive
definite (possible for the skeptic and for pairwise data) is adjusted to
the nearest positive-definite matrix before estimation.

The penalty path, the EBIC and the certificate are those of
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
with its default settings; the certificate is computed by
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
on the transformed correlation matrix. After a positive `threshold` the
certificate is missing; [`print()`](https://rdrr.io/r/base/print.html)
still shows the value for the network before thresholding.

## Conditions

`psychnet_bad_input` when an argument is out of range, `threshold` is
logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
matrix with a positive diagonal, `n` is missing with `cor_matrix`,
`data` has fewer than two usable variables, a pair of variables is never
observed together, `labels` do not match the nodes, or `native = FALSE`
without the `glasso` package; `psychnet_constant_variable` (warning)
when a constant or all-missing column is dropped.

## References

Liu, H., Lafferty, J., & Wasserman, L. (2009). The nonparanormal:
Semiparametric estimation of high dimensional undirected graphs.
*Journal of Machine Learning Research*, 10, 2295-2328.

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

## See also

[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
for the untransformed model;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
for the check;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "huge"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

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
