# Node predictability

Reports how well each node is predicted by the others in a fitted
network. For Gaussian graphical models this is the closed-form variance
explained (R-squared) from the precision matrix and needs no data. For
the nodewise models
([`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md),
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md))
it requires the data and reports R-squared for Gaussian nodes and
classification accuracy (`CC`) plus normalized accuracy (`nCC`) for
binary nodes.

## Usage

``` r
net_predict(x, data = NULL, ...)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object.

- data:

  The data the network was estimated from; required for the nodewise
  models (ising / IsingSampler / mgm), ignored for the GGMs.

- ...:

  Unused.

## Value

A tidy `data.frame`, one row per node, with columns `node`, `type`
(`"gaussian"` or `"binary"`), `metric` (`"R2"` or `"nCC"`),
`predictability`, and `accuracy` (classification accuracy for binary
nodes, `NA` for Gaussian).

## Details

Prediction selects the original modeled columns, independently of
display `labels`. Binary MGM categories use the levels recorded when
fitting, including numeric codes such as 1/2; unseen categories are
refused.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853-861.
[doi:10.3758/s13428-017-0910-x](https://doi.org/10.3758/s13428-017-0910-x)

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_predict(ebic_glasso(cor_matrix = S, n = 250))
#>   node     type metric predictability accuracy
#> 1   V1 gaussian     R2      0.1568160       NA
#> 2   V2 gaussian     R2      0.2711166       NA
#> 3   V3 gaussian     R2      0.2711166       NA
#> 4   V4 gaussian     R2      0.2711166       NA
#> 5   V5 gaussian     R2      0.2711166       NA
#> 6   V6 gaussian     R2      0.1568160       NA
# Nodewise models need the data they were fitted to.
b <- dichotomize(SRL_GPT)
net_predict(ising_fit(b), data = b)
#>   node   type metric predictability  accuracy
#> 1  CSU binary    nCC      0.7744361 0.9000000
#> 2   IV binary    nCC      0.7769784 0.8966667
#> 3   SE binary    nCC      0.6641791 0.8500000
#> 4   SR binary    nCC      0.7101449 0.8666667
#> 5   TA binary    nCC      0.0000000 0.5000000
```
