# Correctness certificate of a fitted network

Every regularized or constrained estimator in `psychnet` self-certifies:
it reports how far the returned network sits from the unique optimum of
its own convex objective (a KKT / stationarity residual), or – for the
structural methods – whether the graph satisfies the identity that
defines it. This verb returns that certificate as a tidy one-row
`data.frame`, so correctness is read the same way for every method.

## Usage

``` r
certificate(x, tol = 1e-06)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object.

- tol:

  Tolerance below which the fit is flagged `certified = TRUE`. Default
  `1e-6`.

## Value

A one-row `data.frame` with columns `method`, `certificate` (the
residual; smaller is better), `kind` (`"kkt"` for the optimization
certificates, `"structural"` for TMFG/relimp, `"none"` for cor/pcor),
and `certified` (logical: residual at or below `tol`).

## Details

The residual is near machine zero for a correctly solved problem. `cor`
and `pcor` have no optimization to certify and report `NA`.

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
