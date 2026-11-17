# Estimate a directed relative-importance network

`relimp_network()` estimates a directed relative-importance network: for
each variable, how much of its variance the other variables explain
together and how that explained variance is shared among them. Each
variable is regressed on all the others, and its \\R^2\\ is divided
among the predictors by the LMG method (Lindeman et al., 1980;
Groemping, 2006): a predictor's share is its increase in \\R^2\\ when it
enters the regression, averaged over all orders in which the predictors
can enter, which is the Shapley value of the predictor.

The shares are non-negative and sum to the \\R^2\\. The network is
directed: the edge from a predictor to an outcome carries the
predictor's share of the outcome's explained variance, so the two edges
between a pair of variables can differ. `relimp_network()` computes
these shares for every variable from raw data or a correlation matrix.
With `normalized = TRUE` each outcome's shares are rescaled to sum to 1,
giving proportions of explained variance rather than proportions of
total variance. The computation considers every subset of predictors, so
its cost doubles with each added variable; `max_nodes` limits the
network size.

It returns the network of shares; the identity that defines the
decomposition is checked by
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md).

For undirected conditional associations use
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
or
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md).
The shares are those of `relaimpo::calc.relimp(type = "lmg")` applied to
each variable in turn.

## Usage

``` r
relimp_network(
  data = NULL,
  cor_matrix = NULL,
  n = NULL,
  cor_method = c("pearson", "spearman", "kendall", "auto"),
  max_nodes = 21L,
  normalized = FALSE,
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  A numeric data frame or matrix with one row per observation. Not
  needed when `cor_matrix` is given.

- cor_matrix:

  A correlation matrix to use instead of `data`; it must be symmetric
  and positive semi-definite, and a covariance matrix is rescaled to a
  correlation matrix.

- n:

  Accepted and ignored: the decomposition needs only the correlation
  matrix. It exists so that a call written for the other estimators,
  `(cor_matrix = , n = )`, also works here.

- cor_method:

  Correlation computed from `data`: `"pearson"` (default), `"spearman"`,
  `"kendall"`, or `"auto"` (polychoric and polyserial correlations for
  ordinal items; see
  [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)).

- max_nodes:

  Largest number of variables accepted, a whole number \>= 2. Default
  `21`; larger networks are refused.

- normalized:

  Whether to rescale each outcome's incoming shares to sum to 1, as
  `bootnet` and `relaimpo` report them. Default `FALSE`, which keeps
  shares that sum to the outcome's \\R^2\\.

- na_method:

  Handling of missing values in `data`: `"pairwise"` (default) or
  `"listwise"`. See
  [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md).

- labels:

  Optional node labels, one unique label per node of the network.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)),
directed, whose edge weights are importance shares.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per edge with columns `from` (the predictor), `to` (the outcome)
and `weight`;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the structural check. The object also stores each node's
full-model \\R^2\\ (`r2`), the correlation matrix (`cor_matrix`),
whether the shares were normalized (`normalized`) and, with
`normalized = TRUE`, the unscaled shares (`raw_importance`), from which
the certificate is computed. The sample size is recorded when the
network is built from `data` and is `NA` with `cor_matrix`.

## Details

The \\R^2\\ of each subset regression is computed from the correlation
matrix. The certificate, computed by
[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md),
checks that each outcome's shares sum to its full-model \\R^2\\ and that
this \\R^2\\ agrees with its closed form from the inverse correlation
matrix.

## Conditions

`psychnet_bad_input` when there are more variables than `max_nodes`,
`max_nodes` or `normalized` is invalid, `cor_matrix` is not a finite,
symmetric, positive semi-definite matrix with a positive diagonal, or
`labels` do not match the nodes; `psychnet_constant_variable` (warning)
when a constant or all-missing column is dropped.

## References

Lindeman, R. H., Merenda, P. F., & Gold, R. Z. (1980). *Introduction to
bivariate and multivariate analysis*. Scott, Foresman.

Groemping, U. (2006). Relative importance for linear regression in R:
The package relaimpo. *Journal of Statistical Software*, 17(1), 1-27.
[doi:10.18637/jss.v017.i01](https://doi.org/10.18637/jss.v017.i01)

## See also

[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md)
and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
for the check;
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
and
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)
for undirected networks;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md),
which routes `method = "relimp"` here;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the result.

## Examples

``` r
S <- 0.4^abs(outer(1:5, 1:5, "-"))
relimp_network(cor_matrix = S)
#> <psychnet> relimp network
#>   nodes: 5   edges: 20   (directed)
#>   structural certificate (decomposition residual): 5.55e-17
```
