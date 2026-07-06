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

  :   `1` for the upper half of the column by rank, giving a balanced
      (~50/50) split that is robust to ties (useful for coarse Likert
      items where a median split is badly unbalanced).

## Value

An integer matrix of `0`/`1` values with the same dimensions and
dimnames as `data`.

## Examples

``` r
b <- dichotomize(SRL_GPT, method = "median")
table(b)                          # values are 0/1 only
#> b
#>   0   1 
#> 694 806 
```
