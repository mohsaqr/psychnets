# Split numeric items into binary 0/1 variables

`dichotomize()` splits numeric items into binary 0/1 variables, each
column independently at its median (default), at its mean, or by rank.
Ising network models
([`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md))
describe binary variables, so Likert-type or continuous items must be
split before such a model can be fitted; the resulting network then
describes associations between scoring high and scoring low on the items
rather than between their full response scales.

Missing values stay missing. The result is an integer 0/1 matrix with
the dimensions of the input.

Splitting discards the information in the response scale. To keep it,
estimate a Gaussian graphical model, for example
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
with `cor_method = "auto"` (polychoric correlations for ordinal items),
or a mixed graphical model with
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md).

## Usage

``` r
dichotomize(data, method = c("median", "mean", "rank"))
```

## Arguments

- data:

  A numeric matrix or data frame with one row per observation.

- method:

  Split rule, applied to each column:

  `"median"`

  :   (default) `1` when the value is at or above the column median.

  `"mean"`

  :   `1` when the value is above the column mean.

  `"rank"`

  :   `1` when the value's average rank is above the midpoint. Tied
      values stay together, so the split can be unbalanced but never
      depends on row order.

  On an item with many responses at its floor, the median can equal the
  minimum, so `"median"` puts every observation in the upper group and
  the output column is constant (an Ising fit would drop it); `"mean"`
  or `"rank"` splits such items.

## Value

An integer matrix of 0/1 values with the dimensions and dimnames of
`as.matrix(data)`; missing values stay `NA`.

## Conditions

`psychnet_bad_input` when `data` is not numeric or has no columns, a
column is constant or entirely missing, or `method` is not one of the
rules; `psychnet_constant_variable` (warning) when a split puts every
observation of a column in one group.

## See also

[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
and
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md),
which need binary data.

## Examples

``` r
b <- dichotomize(SRL_GPT, method = "median")
table(b)                          # values are 0/1 only
#> b
#>   0   1 
#> 694 806 
```
