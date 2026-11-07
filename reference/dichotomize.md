# Dichotomize numeric columns to 0/1

Splits each column of a numeric matrix or data frame into a binary 0/1
variable. This is the usual preprocessing step before fitting an Ising
network
([`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md))
to Likert or other ordinal/continuous data, which require binary input.

## Usage

``` r
dichotomize(data, method = c("median", "mean", "rank"))
```

## Arguments

- data:

  Numeric matrix or data frame (rows = observations).

- method:

  Split rule, applied independently to each column:

  `"median"`

  :   (default) `1` if the value is `>=` the column median.

  `"mean"`

  :   `1` if the value is `>` the column mean.

  `"rank"`

  :   `1` for observations whose average rank is above the midpoint.
      Tied values always remain together, so the split may be unbalanced
      rather than depending on row order.

  On a floor- or ceiling-heavy item the median can equal the column
  minimum, so `"median"` (`>=`) puts every observation in the upper
  group and the output column is constant. Such columns are reported by
  a `psychnet_constant_variable` warning (an Ising fit would drop them);
  `method = "mean"` or `"rank"` splits them.

## Value

An integer matrix of `0`/`1` values with the same dimensions and
dimnames as `data`. Missing values stay `NA` in the output (they are not
assigned to either group). A non-numeric `data`, an empty one, or a
constant / all-missing input column is an error of class
`psychnet_bad_input`.

## Examples

``` r
b <- dichotomize(SRL_GPT, method = "median")
table(b)                          # values are 0/1 only
#> b
#>   0   1 
#> 694 806 
```
