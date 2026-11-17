# Select an unregularized Gaussian graphical model by stepwise search

`ggm_modselect()` selects an unregularized Gaussian graphical model: it
chooses which edges are present and then estimates them by maximum
likelihood, so the partial correlations it reports are not shrunk. In a
Gaussian graphical model variables are nodes joined by partial
correlations, with a missing edge meaning conditional independence; a
regularized estimate such as the graphical lasso instead shrinks every
retained edge towards zero, and Williams et al. (2019) discuss when the
unregularized approach is preferable to regularization.

`ggm_modselect()` generates candidate graphs from the graphical lasso
over a sequence of penalties, fits the unregularized maximum-likelihood
model on each candidate and keeps the one with the lowest extended
Bayesian information criterion (EBIC). With `stepwise = TRUE` (default)
it then adds or removes single edges as long as the EBIC improves.

The returned network is the maximum-likelihood network on the selected
graph, which
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
checks.

For a regularized network use
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md);
for all partial correlations without selection use
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md).
The procedure is the one of
[`qgraph::ggmModSelect()`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html).

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

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given.

- cor_matrix:

  A correlation matrix to use instead of `data`; a covariance matrix is
  rescaled to a correlation matrix. `n` is then required.

- n:

  Sample size, required with `cor_matrix`: a single finite number= 1.
  Ignored when data is given.

- gamma:

  EBIC hyperparameter, a number \>= 0. Default `0` (the ordinary BIC).

- stepwise:

  Whether to refine the best candidate graph by adding and removing
  single edges. Default `TRUE`.

- nlambda:

  Number of graphical lasso penalties used to generate candidate graphs,
  a whole number \>= 2. Default `100`.

- lambda_min_ratio:

  Smallest penalty as a fraction of the largest, in (0, 1). Default
  `0.01`.

- threshold:

  Absolute partial correlation below which an edge is set to zero after
  estimation: a single number \>= 0. Default `0`. A logical value is
  refused. A positive value leaves the certificate missing.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) or
  `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- native:

  Solver for the candidate graphs: `TRUE` (default) the package's own
  graphical lasso, `FALSE` the `glasso` package, which must be
  installed. The returned network is the unregularized fit either way.

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
(`precision`), the selected graph as a logical matrix (`support`),
`gamma`, the correlation matrix (`cor_matrix`) and the EBIC of the
selected graph (`ebic`).

## Details

The EBIC of a graph with \\E\\ edges is \\-2L + E \log n + 4 \gamma E
\log p\\, where \\L\\ is the log-likelihood of the maximum-likelihood
fit on that graph. Because the retained edges are not shrunk, the
default `gamma = 0` (the ordinary BIC) is the intended setting, as in
[`qgraph::ggmModSelect()`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html);
the regularized
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)
default to `0.5`.

The stepwise search visits every pair of variables in turn, adds or
removes the edge, and keeps the change whenever the EBIC decreases; it
stops when a full pass brings no improvement (at most 200 passes). A
correlation matrix with no off-diagonal association returns the empty
network without a search.

The certificate is computed by
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md).
After a positive `threshold` it is missing;
[`print()`](https://rdrr.io/r/base/print.html) still shows the value for
the network before thresholding.

## Conditions

`psychnet_bad_input` when an argument is out of range, `threshold` is
logical, `cor_matrix` is not a finite, symmetric, positive semi-definite
matrix with a positive diagonal, `n` is missing with `cor_matrix`,
`data` has fewer than two usable variables, `labels` do not match the
nodes, or `native = FALSE` without the `glasso` package;
`psychnet_no_converge` (error) when no candidate graph can be fitted;
`psychnet_fit_failed` (warning) when some candidate fits failed and were
skipped; `psychnet_no_converge` (warning) when the selected fit's
certificate exceeds `1e-6`; `psychnet_constant_variable` (warning) when
a constant or all-missing column is dropped.

## References

Williams, D. R., Rhemtulla, M., Wysocki, A. C., & Rast, P. (2019). On
nonregularized estimation of psychological networks. *Multivariate
Behavioral Research*, 54(5), 719-750.
[doi:10.1080/00273171.2019.1575716](https://doi.org/10.1080/00273171.2019.1575716)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

## See also

[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
for the regularized network;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
for the check;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "ggm"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

## Examples

``` r
S <- 0.5^abs(outer(1:6, 1:6, "-"))
ggm_modselect(cor_matrix = S, n = 250)
#> <psychnet> ggm network
#>   nodes: 6   edges: 5   (undirected)
#>   optimality (KKT residual): 2.22e-16
```
