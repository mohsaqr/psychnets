# Correctness certificate of a fitted network

Regularized or constrained estimators in `psychnets` self-certify: each
fitted network carries a residual that reports how far it sits from the
optimum of its own objective (a KKT / stationarity residual), or – for
the structural methods (TMFG, relative importance) – whether the graph
satisfies the identity that defines it. This verb returns that
certificate as a tidy one-row `data.frame`, so correctness is read the
same way for every method.

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
and `certified` (logical: residual at or below `tol`; `NA` when no
applicable certificate is available). Correlation networks use kind
`"none"` and are reported as `TRUE` because no optimization claim is
made.

## Details

The residual is near machine zero for a correctly solved problem. A
graph altered after fitting (post-estimation `threshold`, Wald pruning
in
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md))
or an external-engine fit without an available residual reports `NA`
rather than being presented as certified. `cor` and `pcor` have no
optimization to certify and report `NA` with kind `"none"`.

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
