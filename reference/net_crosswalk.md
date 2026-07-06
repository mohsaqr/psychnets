# Argument crosswalk: psychnet as a substitute for qgraph / IsingFit / mgm

A tidy, one-row-per-argument map from each reference package's estimator
to its `psychnet` equivalent, so users migrating from
[`qgraph::EBICglasso`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html),
[`qgraph::cor_auto`](https://rdrr.io/pkg/qgraph/man/cor_auto.html),
[`qgraph::ggmModSelect`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html),
[`IsingFit::IsingFit`](https://rdrr.io/pkg/IsingFit/man/isingfit.html),
or [`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html) can find the
matching argument and see what `psychnet` changes by default.
Cross-sectional estimators only (temporal models are out of scope).

## Usage

``` r
net_crosswalk(
  reference = c("all", "EBICglasso", "cor_auto", "ggmModSelect", "IsingFit", "mgm")
)
```

## Arguments

- reference:

  Which reference function to show: `"all"` (default), `"EBICglasso"`,
  `"cor_auto"`, `"ggmModSelect"`, `"IsingFit"`, or `"mgm"`.

## Value

A tidy `data.frame`, one row per argument, with columns `reference` (the
`pkg::fn` being substituted), `psychnet` (the psychnet verb), `ref_arg`,
`psychnet_arg` (`"-"` when there is no counterpart), `status`
(`identical` / `renamed` / `default differs` / `semantics differ` /
`reference only` / `psychnet only`), and a short `note`.

## Examples

``` r
net_crosswalk("EBICglasso")
#>             reference    psychnet           ref_arg     psychnet_arg
#> 1  qgraph::EBICglasso ebic_glasso                 S       cor_matrix
#> 2  qgraph::EBICglasso ebic_glasso                 n                n
#> 3  qgraph::EBICglasso ebic_glasso             gamma            gamma
#> 4  qgraph::EBICglasso ebic_glasso           nlambda          nlambda
#> 5  qgraph::EBICglasso ebic_glasso  lambda.min.ratio lambda_min_ratio
#> 6  qgraph::EBICglasso ebic_glasso         threshold        threshold
#> 7  qgraph::EBICglasso ebic_glasso penalize.diagonal                -
#> 8  qgraph::EBICglasso ebic_glasso             refit                -
#> 9  qgraph::EBICglasso ebic_glasso           checkPD                -
#> 10 qgraph::EBICglasso ebic_glasso    penalizeMatrix                -
#> 11 qgraph::EBICglasso ebic_glasso     countDiagonal                -
#> 12 qgraph::EBICglasso ebic_glasso  returnAllResults                -
#> 13 qgraph::EBICglasso ebic_glasso           verbose                -
#> 14 qgraph::EBICglasso ebic_glasso                 -             data
#> 15 qgraph::EBICglasso ebic_glasso                 -       cor_method
#> 16 qgraph::EBICglasso ebic_glasso                 -        na_method
#> 17 qgraph::EBICglasso ebic_glasso                 -           native
#> 18 qgraph::EBICglasso ebic_glasso                 -           labels
#>              status                                                      note
#> 1           renamed     correlation matrix in; psychnet also accepts raw data
#> 2         identical                                                          
#> 3         identical                                  EBIC tuning, default 0.5
#> 4         identical                                               default 100
#> 5           renamed                                         same default 0.01
#> 6  semantics differ qgraph: logical sig-rule; psychnet: numeric weight cutoff
#> 7    reference only                       psychnet always FALSE (W_ii = S_ii)
#> 8    reference only     psychnet always refits the selected lambda (two-tier)
#> 9    reference only            psychnet always validates/PD-checks cor_matrix
#> 10   reference only         per-edge penalty matrix; psychnet uses scalar rho
#> 11   reference only        EBIC df toggle; psychnet never counts the diagonal
#> 12   reference only                 psychnet always stores $kkt, $lambda, ...
#> 13   reference only                                        psychnet is silent
#> 14    psychnet only                           estimate straight from raw data
#> 15    psychnet only                             pearson/spearman/kendall/auto
#> 16    psychnet only                            pairwise/listwise missing data
#> 17    psychnet only             native solver (TRUE) / glasso Fortran (FALSE)
#> 18    psychnet only                                               node labels
net_crosswalk("IsingFit")
#>             reference  psychnet           ref_arg     psychnet_arg
#> 1  IsingFit::IsingFit ising_fit                 x             data
#> 2  IsingFit::IsingFit ising_fit             gamma            gamma
#> 3  IsingFit::IsingFit ising_fit               AND             rule
#> 4  IsingFit::IsingFit ising_fit           min_sum          min_sum
#> 5  IsingFit::IsingFit ising_fit            family                -
#> 6  IsingFit::IsingFit ising_fit              plot                -
#> 7  IsingFit::IsingFit ising_fit       progressbar                -
#> 8  IsingFit::IsingFit ising_fit lowerbound.lambda                -
#> 9  IsingFit::IsingFit ising_fit                 -          nlambda
#> 10 IsingFit::IsingFit ising_fit                 - lambda_min_ratio
#> 11 IsingFit::IsingFit ising_fit                 -          weights
#> 12 IsingFit::IsingFit ising_fit                 -        na_method
#> 13 IsingFit::IsingFit ising_fit                 -           labels
#>             status                                   note
#> 1          renamed                        binary 0/1 data
#> 2        identical                           default 0.25
#> 3          renamed      AND=TRUE/FALSE -> rule='AND'/'OR'
#> 4  default differs      -Inf vs NULL (both keep all rows)
#> 5   reference only              psychnet is binomial-only
#> 6   reference only                   psychnet never plots
#> 7   reference only           psychnet has no progress bar
#> 8   reference only psychnet uses its own full lambda path
#> 9    psychnet only                     lambda-path length
#> 10   psychnet only              smallest penalty fraction
#> 11   psychnet only                    observation weights
#> 12   psychnet only                      pairwise/listwise
#> 13   psychnet only                            node labels
```
