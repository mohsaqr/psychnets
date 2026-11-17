# Estimate a Gaussian graphical model with the EBIC graphical lasso

`ebic_glasso()` estimates a Gaussian graphical model: a network in which
variables are nodes joined by partial correlations, and a missing edge
means two variables are independent given all the others. It fits the
graphical lasso, which shrinks weak partial correlations to zero, and
chooses the penalty by the extended Bayesian information criterion
(Epskamp & Fried, 2018).

An edge is the association between two variables after controlling for
all other variables in the network. The penalty sets weak partial
correlations to exactly zero, so that only the associations supported by
the data are retained, and the EBIC balances the fit of the model
against the number of edges. `ebic_glasso()` fits the graphical lasso
over a sequence of penalties, selects the penalty with the lowest EBIC
and returns the partial-correlation network at that penalty. It takes
raw data, or a correlation matrix together with its sample size. `gamma`
sets how strongly the EBIC favours sparse networks,
`refit = "unregularized"` keeps the selected edges but estimates them
without shrinkage, and `threshold` removes weak edges after estimation.

The returned network is the optimum of the penalized likelihood at the
selected penalty, which
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
checks.

This regularized partial-correlation network is the standard estimator
for continuous and ordinal data in psychological network analysis. For
data far from normality,
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)
fits the same model after a rank-based transformation.
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
selects an unregularized network by stepwise search,
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)
keeps every partial correlation, and
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
estimate networks of binary and mixed variables. The model is the one
fitted by
[`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html).

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

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given.

- cor_matrix:

  A correlation matrix to use instead of `data`; a covariance matrix is
  rescaled to a correlation matrix. It must be symmetric and positive
  semi-definite. `n` is then required.

- n:

  Sample size, required with `cor_matrix`: a single finite number= 1.
  Ignored when data is given.

- gamma:

  EBIC hyperparameter, a number \>= 0; larger values favour sparser
  networks and 0 gives the ordinary BIC. Default `0.5`.

- nlambda:

  Number of penalties scanned, a whole number \>= 2. Default `100`.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest, in (0, 1). Default
  `0.01`.

- threshold:

  Absolute partial correlation below which an edge is set to zero after
  estimation: a single number \>= 0. Default `0`. A logical value (the
  `threshold = TRUE` rule of `qgraph`) is refused. A positive value
  leaves the certificate missing.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) uses
  pairwise-complete correlations adjusted to the nearest
  positive-definite matrix; `"listwise"` drops incomplete rows. The two
  agree for complete data.

- native:

  `TRUE` (default) uses the package's own solver. `FALSE` solves each
  penalty with the `glasso` package, which must be installed, to
  reproduce results of the `glasso` and `qgraph` packages.

- lambda_path:

  Optional numeric vector of positive penalties to scan, in the given
  order, replacing `nlambda` and `lambda_min_ratio`. Useful when
  resampling: the path of the full data, reused for every resample,
  keeps the selected penalties comparable. Default `NULL`.

- penalize_diagonal:

  Whether the penalty also applies to the diagonal of the precision
  matrix. Default `FALSE`.

- refit:

  What to return at the selected penalty: `TRUE` (default) the optimum
  at that penalty; `FALSE` the fit obtained while scanning the
  penalties, less precise and kept to reproduce earlier results;
  `"unregularized"` the unpenalized maximum-likelihood network
  restricted to the selected edges, without shrinkage.

- labels:

  Optional node labels, one unique label per node of the network (the
  columns left after non-numeric and constant columns are dropped).

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
whose edge weights are partial correlations.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from`, `to` and `weight`;
[`print()`](https://rdrr.io/r/base/print.html) shows the selected
penalty and the certificate;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the check. The object also stores the precision matrix
(`precision`), the selected penalty (`lambda`), `gamma`, the correlation
matrix used (`cor_matrix`), the EBIC of the selected network (`ebic`),
the penalties scanned and their EBIC (`lambda_path`, `ebic_path`) and,
with `refit = "unregularized"`, the selected graph (`support`). Nodes
keep the column order of the input.

## Details

For a precision matrix \\\Theta\\ with \\E\\ nonzero off-diagonal
entries, fitted to a correlation matrix \\S\\ of \\p\\ variables and
\\n\\ observations, the criterion is \$\$\mathrm{EBIC} = -2L + E \log
n + 4 \gamma E \log p, \quad L = \frac{n}{2}\left(\log\det\Theta -
\mathrm{tr}(S\Theta)\right).\$\$ The penalties form a log-spaced
sequence of `nlambda` values from the largest absolute off-diagonal
correlation, at which the network is empty, down to `lambda_min_ratio`
times that value. When two penalties give the same EBIC, the one earlier
in the path is selected, which on the default path is the larger penalty
and the sparser network. Unless `lambda_path` is given, a correlation
matrix with no off-diagonal association returns the empty network with
penalty 0.

With `na_method = "pairwise"` and incomplete data, the sample size
entering the EBIC is the median number of observations available for a
pair of variables.

The certificate is computed by
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md).
With `refit = "unregularized"` the network is no longer a penalized
optimum and is checked by
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
instead. With `refit = FALSE` or `native = FALSE` the network is solved
to a looser tolerance, and its certificate can exceed the default
tolerance of
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md).
After a positive `threshold` the certificate is missing;
[`print()`](https://rdrr.io/r/base/print.html) still shows the value for
the network before thresholding.

## Conditions

`psychnet_bad_input` when an argument is out of range, `threshold` is
logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
matrix with a positive diagonal, `n` is missing with `cor_matrix`,
`data` has fewer than two usable variables, a pair of variables is never
observed together, `labels` do not match the nodes, `native = FALSE`
without the `glasso` package, or every penalty fails to fit;
`psychnet_constant_variable` (warning) when a constant or all-missing
column is dropped. A `psychnet_dropped_columns` message names
non-numeric columns that were dropped.

## References

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
correlation networks. *Psychological Methods*, 23(4), 617-634.
[doi:10.1037/met0000167](https://doi.org/10.1037/met0000167)

## See also

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
for the check;
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md),
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
and
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)
for related Gaussian networks;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "glasso"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

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
certificate(fit)
#>   method certificate kind certified
#> 1 glasso 1.10463e-10  kkt      TRUE
```
