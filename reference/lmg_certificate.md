# Check the decomposition of a relative-importance network

`lmg_certificate()` checks the decomposition stored in a
relative-importance network from
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md).
Such a network divides the \\R^2\\ of each variable, regressed on all
the others, among its predictors; it is defined by a decomposition
rather than as the optimum of an objective, and two identities must hold
for that decomposition to be correct.

The shares a node receives must sum exactly to its \\R^2\\, which is the
efficiency property of the Shapley value. And that \\R^2\\ must be the
one implied by the correlation matrix, which has a closed form. Together
they show that the stored shares add up to the variance each node's
predictors actually explain.

`lmg_certificate()` returns the larger of the two deviations: a value
near zero confirms the decomposition, up to numerical precision.

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the same value for a fitted network. The check does not verify
the \\R^2\\ of the individual subset regressions between the empty and
the full model, and the sum identity alone would hold for any set of
subset values; the closed form is what ties the full model to the data.

## Usage

``` r
lmg_certificate(x)
```

## Arguments

- x:

  A relative-importance network from
  [`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md)
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

## Value

A single number, the larger of the two deviations; 0 for an exact
decomposition.

## Details

For outcome \\j\\ with shares \\w\_{kj}\\ and full-model \\R^2_j\\, the
two deviations are \\\max_j \|\sum_k w\_{kj} - R^2_j\|\\ and \\\max_j
\|R^2_j - (1 - 1 / (S^{-1})\_{jj})\|\\. When the correlation matrix
\\S\\ is singular (reciprocal condition number below `1e-12`) the closed
form does not exist and only the first deviation is returned. For a
network with `normalized = TRUE` the unscaled shares are used.

## Conditions

`psychnet_bad_input` when `x` is not a network from
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md).

## See also

[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md);
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md);
[`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md),
the other structural check.

## Examples

``` r
S <- 0.4^abs(outer(1:5, 1:5, "-"))
lmg_certificate(relimp_network(cor_matrix = S))
#> [1] 5.551115e-17
```
