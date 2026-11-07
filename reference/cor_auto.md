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

  Maximum distinct values for a variable to count as ordinal: a whole
  number \>= 2. Default 7.

- na_method:

  `"pairwise"` (default) or `"listwise"`.

## Value

A correlation matrix with the variable names as dimnames. An invalid
`ordinal_max_levels` is an error of class `psychnet_bad_input`.

## Details

The polychoric correlation is the two-step estimator of Olsson (1979):
thresholds fixed at the marginal normal quantiles, then the
contingency-table likelihood maximised over the correlation; it matches
`psych`/`qgraph` to about 1e-5. The polyserial correlation is the
Olsson, Drasgow & Dorans (1982) *ad-hoc* estimator (the closed form
\\r\_{xy}\\\sigma_y / \sum_k \phi(\tau_k)\\), not the two-step
maximum-likelihood polyserial that
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
obtains through `lavaan`. Both are consistent; on finite samples they
differ by up to about 0.006, so a mixed ordinal / continuous correlation
matrix will not reproduce
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
exactly.

An ordered-factor column of a data frame enters as its integer level
codes and is then ordinal whenever it has at most `ordinal_max_levels`
levels. Unordered factors, character and logical columns are dropped
with a `psychnet_dropped_columns` message; constant columns with a
`psychnet_constant_variable` warning.

## References

Olsson, U. (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. *Psychometrika*, 44(4), 443-460.
[doi:10.1007/BF02296207](https://doi.org/10.1007/BF02296207)

Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial
correlation coefficient. *Psychometrika*, 47(3), 337-347.
[doi:10.1007/BF02294164](https://doi.org/10.1007/BF02294164)

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
