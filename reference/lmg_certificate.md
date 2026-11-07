# Relative-importance (LMG / Shapley) certificate

Checks two identities. By Shapley efficiency, the importance shares a
node receives from its predictors sum exactly to that node's full-model
R-squared; and that full-model R-squared must equal the closed form
\\1 - 1 / (S^{-1})\_{jj}\\. The first alone cannot detect a wrong subset
solver (efficiency holds for any set function), so the second anchors
the full model independently. When `S` is singular (reciprocal condition
number below 1e-12) the closed form does not exist and only efficiency
is checked. The certificate does not verify the individual subset
R-squared values between the empty and the full model.

## Usage

``` r
lmg_certificate(x)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object produced by
  [`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md).

## Value

The larger of the two deviations (scalar); 0 = exact decomposition. An
`x` that is not a relimp network is an error of class
`psychnet_bad_input`.

## Examples

``` r
S <- 0.4^abs(outer(1:5, 1:5, "-"))
lmg_certificate(relimp_network(cor_matrix = S))
#> [1] 5.551115e-17
```
