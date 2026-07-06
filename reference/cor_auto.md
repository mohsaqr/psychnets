# Automatic correlation matrix (polychoric / polyserial / Pearson)

Detects ordinal variables (integer-valued with at most
`ordinal_max_levels` levels) and returns the correlation matrix using a
polychoric correlation for ordinal-ordinal pairs, a polyserial
correlation for ordinal-continuous pairs, and Pearson otherwise,
projected to the nearest positive-definite matrix. The base-R
counterpart of
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html);
this is the correlation `bootnet`/`qgraph` use by default for Likert
data.

## Usage

``` r
cor_auto(data, ordinal_max_levels = 7L, na_method = c("pairwise", "listwise"))
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations).

- ordinal_max_levels:

  Maximum distinct values for a variable to count as ordinal. Default 7.

- na_method:

  `"pairwise"` (default) or `"listwise"`.

## Value

A correlation matrix with the variable names as dimnames.

## Examples

``` r
set.seed(1)
z <- matrix(stats::rnorm(300 * 4), 300, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
x <- apply(z, 2, function(col) as.integer(cut(col, 5)))   # 5-level Likert
cor_auto(x)
#>            V1        V2        V3         V4
#> V1 1.00000000 0.4542199 0.2843910 0.09885485
#> V2 0.45421990 1.0000000 0.4944493 0.32126907
#> V3 0.28439095 0.4944493 1.0000000 0.50448731
#> V4 0.09885485 0.3212691 0.5044873 1.00000000
```
