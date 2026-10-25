# Argument crosswalk: psychnet as a substitute for qgraph / IsingFit / mgm / psychonetrics

A tidy, one-row-per-argument map from each reference package's estimator
to its `psychnet` equivalent, so users migrating from
[`qgraph::EBICglasso`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html),
[`qgraph::cor_auto`](https://rdrr.io/pkg/qgraph/man/cor_auto.html),
[`qgraph::ggmModSelect`](https://rdrr.io/pkg/qgraph/man/ggmModSelect.html),
[`IsingFit::IsingFit`](https://rdrr.io/pkg/IsingFit/man/isingfit.html),
[`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html), or
[`psychonetrics::panelgvar`](https://rdrr.io/pkg/psychonetrics/man/panelvar_family.html)
can find the matching argument and see what `psychnet` changes by
default. Intensive time-series models (`graphicalVAR`, `mlVAR`) are out
of scope; the few-wave panel model is covered by
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md).

## Usage

``` r
net_crosswalk(
  reference = c("all", "EBICglasso", "cor_auto", "ggmModSelect", "IsingFit", "mgm",
    "panelgvar")
)
```

## Arguments

- reference:

  Which reference function to show: `"all"` (default), `"EBICglasso"`,
  `"cor_auto"`, `"ggmModSelect"`, `"IsingFit"`, `"mgm"`, or
  `"panelgvar"`.

## Value

A tidy `data.frame`, one row per argument, with columns `reference` (the
`pkg::fn` being substituted), `psychnet` (the psychnet verb), `ref_arg`,
`psychnet_arg` (`"-"` when there is no counterpart), `status`
(`identical` / `renamed` / `default differs` / `semantics differ` /
`reference only` / `psychnet only`), and a short `note`.

## Examples

``` r
net_crosswalk("EBICglasso")
#>             reference    psychnet           ref_arg      psychnet_arg
#> 1  qgraph::EBICglasso ebic_glasso                 S        cor_matrix
#> 2  qgraph::EBICglasso ebic_glasso                 n                 n
#> 3  qgraph::EBICglasso ebic_glasso             gamma             gamma
#> 4  qgraph::EBICglasso ebic_glasso           nlambda           nlambda
#> 5  qgraph::EBICglasso ebic_glasso  lambda.min.ratio  lambda_min_ratio
#> 6  qgraph::EBICglasso ebic_glasso         threshold         threshold
#> 7  qgraph::EBICglasso ebic_glasso penalize.diagonal penalize_diagonal
#> 8  qgraph::EBICglasso ebic_glasso             refit             refit
#> 9  qgraph::EBICglasso ebic_glasso           checkPD                 -
#> 10 qgraph::EBICglasso ebic_glasso    penalizeMatrix                 -
#> 11 qgraph::EBICglasso ebic_glasso     countDiagonal                 -
#> 12 qgraph::EBICglasso ebic_glasso  returnAllResults                 -
#> 13 qgraph::EBICglasso ebic_glasso           verbose                 -
#> 14 qgraph::EBICglasso ebic_glasso                 -              data
#> 15 qgraph::EBICglasso ebic_glasso                 -        cor_method
#> 16 qgraph::EBICglasso ebic_glasso                 -         na_method
#> 17 qgraph::EBICglasso ebic_glasso                 -            native
#> 18 qgraph::EBICglasso ebic_glasso                 -            labels
#> 19 qgraph::EBICglasso ebic_glasso                 -       lambda_path
#>              status
#> 1           renamed
#> 2         identical
#> 3         identical
#> 4         identical
#> 5           renamed
#> 6  semantics differ
#> 7           renamed
#> 8  semantics differ
#> 9    reference only
#> 10   reference only
#> 11   reference only
#> 12   reference only
#> 13   reference only
#> 14    psychnet only
#> 15    psychnet only
#> 16    psychnet only
#> 17    psychnet only
#> 18    psychnet only
#> 19    psychnet only
#>                                                                  note
#> 1               correlation matrix in; psychnet also accepts raw data
#> 2                                                                    
#> 3                                            EBIC tuning, default 0.5
#> 4                                                         default 100
#> 5                                                   same default 0.01
#> 6           qgraph: logical sig-rule; psychnet: numeric weight cutoff
#> 7                                    same default FALSE (W_ii = S_ii)
#> 8  qgraph: logical; psychnet also takes "unregularized" (support MLE)
#> 9                      psychnet always validates/PD-checks cor_matrix
#> 10                  per-edge penalty matrix; psychnet uses scalar rho
#> 11                 EBIC df toggle; psychnet never counts the diagonal
#> 12                          psychnet always stores $kkt, $lambda, ...
#> 13                                                 psychnet is silent
#> 14                                    estimate straight from raw data
#> 15                                      pearson/spearman/kendall/auto
#> 16                                     pairwise/listwise missing data
#> 17                      native solver (TRUE) / glasso Fortran (FALSE)
#> 18                                                        node labels
#> 19                         fixed penalty grid, for resampling callers
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
net_crosswalk("panelgvar")
#>                   reference   psychnet          ref_arg     psychnet_arg
#> 1  psychonetrics::panelgvar panel_gvar             data             data
#> 2  psychonetrics::panelgvar panel_gvar             vars             vars
#> 3  psychonetrics::panelgvar panel_gvar            idvar               id
#> 4  psychonetrics::panelgvar panel_gvar          beepvar             wave
#> 5  psychonetrics::panelgvar panel_gvar          missing        na_method
#> 6  psychonetrics::panelgvar panel_gvar            alpha            alpha
#> 7  psychonetrics::panelgvar panel_gvar           adjust           adjust
#> 8  psychonetrics::panelgvar panel_gvar        recursive        recursive
#> 9  psychonetrics::panelgvar panel_gvar        criterion        criterion
#> 10 psychonetrics::panelgvar panel_gvar        estimator           select
#> 11 psychonetrics::panelgvar panel_gvar          nLambda          nlambda
#> 12 psychonetrics::panelgvar panel_gvar lambda_min_ratio lambda_min_ratio
#> 13 psychonetrics::panelgvar panel_gvar      standardize      standardize
#> 14 psychonetrics::panelgvar panel_gvar         datatype                -
#> 15 psychonetrics::panelgvar panel_gvar    within_latent                -
#> 16 psychonetrics::panelgvar panel_gvar   between_latent                -
#> 17 psychonetrics::panelgvar panel_gvar stationary_means                -
#> 18 psychonetrics::panelgvar panel_gvar          covtype                -
#> 19 psychonetrics::panelgvar panel_gvar        optimizer                -
#> 20 psychonetrics::panelgvar panel_gvar                -           select
#> 21 psychonetrics::panelgvar panel_gvar                -           stepup
#> 22 psychonetrics::panelgvar panel_gvar                -            gamma
#> 23 psychonetrics::panelgvar panel_gvar                -            refit
#> 24 psychonetrics::panelgvar panel_gvar                -           labels
#>              status
#> 1         identical
#> 2         identical
#> 3           renamed
#> 4           renamed
#> 5  semantics differ
#> 6         identical
#> 7         identical
#> 8         identical
#> 9         identical
#> 10 semantics differ
#> 11          renamed
#> 12        identical
#> 13 semantics differ
#> 14   reference only
#> 15   reference only
#> 16   reference only
#> 17   reference only
#> 18   reference only
#> 19   reference only
#> 20    psychnet only
#> 21  default differs
#> 22    psychnet only
#> 23    psychnet only
#> 24    psychnet only
#>                                                                                                                                note
#> 1                                                                                            wide (one row per person) or long data
#> 2                                           design matrix (variables x waves); psychnet also takes a character vector for long data
#> 3                                                                                                      person id column (long data)
#> 4                                                                                  wave column (long data); waves ordered by sort()
#> 5                                                                    psychonetrics: FIML when data have NA; psychnet: listwise only
#> 6                                                                                               prune() argument, same default 0.01
#> 7                                                                                             prune() argument, same default "none"
#> 8                                                                                              prune() argument, same default FALSE
#> 9                                                                                             stepup() argument, same default "bic"
#> 10                                             estimator = "PML" (penalized, auto lambda) is select = "ebic"; FIML is not supported
#> 11                                                                                find_penalized_lambda() argument, same default 50
#> 12                                                                              find_penalized_lambda() argument, same default 1e-4
#> 13                                                              "none"/"z"/"z_per_wave" as psychonetrics; "quantile" is not offered
#> 14                                                                                         psychnet infers wide vs long from `vars`
#> 15                                                                                     psychnet fixes "ggm" (the panelgvar default)
#> 16                                                                                     psychnet fixes "ggm" (the panelgvar default)
#> 17                                                                           psychnet: stationary means (the psychonetrics default)
#> 18                                                      psychnet: ML covariance (denominator n), as psychonetrics uses for raw data
#> 19                                                       psychnet: nlminb + Fisher-scoring polish, ML score stored as a certificate
#> 20 "none" (default) = runmodel() full model with every edge's p-value; "prune" = prune() |> stepup(); "ebic" = PML with EBIC lambda
#> 21                                                      stepup() after prune(); psychnet runs it by default (FALSE = prune() alone)
#> 22                                                EBIC gamma for select = "ebic" (psychonetrics encodes it in criterion = "ebic.5")
#> 23                                                                                   unpenalized refit of the EBIC-selected network
#> 24                                                                                                                      node labels
```
