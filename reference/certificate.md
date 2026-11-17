# Check that a fitted network solves its estimator

`certificate()` measures how far a fitted network is from the solution
of its estimator. For a regularized or constrained estimator such as the
graphical lasso, this is the largest violation of the optimality
(Karush–Kuhn–Tucker) conditions; a value near zero means the network is
the optimum, up to numerical precision.

Such an estimator defines the network as the solution of an optimization
problem; the graphical lasso, for example, maximizes a penalized
likelihood. Any solution must satisfy the problem's optimality
conditions, so whether a returned network is that solution can be
checked without a second solver; the amount by which the network
violates the conditions is its certificate. The information-filtering
(TMFG) and relative-importance networks are not optima of an objective;
for them the certificate checks the identity that defines the graph. A
network changed after fitting, by a threshold or by pruning, is no
longer the solution that was checked, and its certificate is reported as
missing.

`certificate()` returns the value together with the kind of certificate
and whether it falls below a tolerance.

Correlation and partial-correlation networks involve no optimization and
have no certificate. To check a precision matrix or regression
coefficients that are not stored in a fitted network, call
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md),
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
or
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
directly.

## Usage

``` r
certificate(x, tol = 1e-06)
```

## Arguments

- x:

  A fitted network, an object of class `psychnet`
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

- tol:

  Tolerance at or below which the network is reported as certified: a
  single finite number \>= 0. Default `1e-6`.

## Value

A data frame with one row and columns `method` (the estimator),
`certificate` (the violation; smaller is better), `kind` (`"kkt"` for an
optimality certificate, `"structural"` for TMFG and relative importance,
`"none"` for correlation and partial-correlation networks) and
`certified` (`TRUE` when `certificate` is at or below `tol`). For a
network changed after fitting, or a fit that did not converge,
`certificate` and `certified` are `NA`. For correlation and
partial-correlation networks `certificate` is `NA` and `certified` is
`TRUE`, because no optimization claim is made.

## Details

The optimality certificates are computed by
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
(graphical lasso networks),
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
(unregularized networks restricted to a graph) and
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
(nodewise Ising and mixed graphical models). The structural certificates
are
[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md),
the number of violated TMFG properties (0 to 4), and
[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md),
the largest difference between the summed importance shares of a node's
predictors and that node's explained variance.

## Conditions

`psychnet_bad_input` when `x` is not a `psychnet` network or `tol` is
not a single finite number \>= 0.

## See also

[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md),
[`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md),
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md),
[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md)
and
[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md)
for the individual checks;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
for the estimators.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
certificate(ebic_glasso(cor_matrix = S, n = 250))
#>   method certificate kind certified
#> 1 glasso 1.10463e-10  kkt      TRUE
certificate(tmfg_network(cor_matrix = S))
#>   method certificate       kind certified
#> 1   tmfg           0 structural      TRUE
```
