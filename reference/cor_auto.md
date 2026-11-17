# Compute correlations suited to ordinal and continuous variables

`cor_auto()` computes a correlation matrix suited to a mix of ordinal
and continuous variables: polychoric correlations between ordinal items,
polyserial correlations between an ordinal item and a continuous
variable, and Pearson correlations otherwise. The polychoric correlation
assumes that each ordinal item is a coarse categorization of a normally
distributed latent variable and estimates the correlation between the
latent variables (Olsson, 1979); the polyserial correlation does the
same for an ordinal item paired with a continuous variable (Olsson et
al., 1982).

Questionnaire items are usually recorded on a few ordered categories,
such as a five-point Likert scale, and treating these codes as
continuous and computing Pearson correlations understates the
association between the underlying traits. `cor_auto()` decides for each
variable whether it is ordinal (integer values with at most
`ordinal_max_levels` distinct values) and computes a polychoric
correlation for each pair of ordinal variables, a polyserial correlation
for each ordinal-continuous pair and a Pearson correlation otherwise.

It returns the resulting matrix adjusted to the nearest
positive-definite correlation matrix, so that it can be used by the
network estimators.

Network estimators for ordinal data commonly start from this correlation
matrix. The estimators use the same correlations when given
`cor_method = "auto"`, for example
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
and
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md).
The function corresponds to
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html).

## Usage

``` r
cor_auto(data, ordinal_max_levels = 7L, na_method = c("pairwise", "listwise"))
```

## Arguments

- data:

  A numeric data frame or matrix with one row per observation.

- ordinal_max_levels:

  Largest number of distinct values for an integer-valued variable to be
  treated as ordinal: a whole number \>= 2. Default `7`.

- na_method:

  Handling of missing values: `"pairwise"` (default) uses the
  observations available for each pair; `"listwise"` drops incomplete
  rows first.

## Value

A symmetric correlation matrix with the variable names as row and column
names, one row and column per variable kept.

## Details

The polychoric correlation is the two-step estimator of Olsson (1979):
the thresholds are fixed at the normal quantiles of the cumulative
marginal proportions, and the likelihood of the contingency table is
then maximized over the correlation. The polyserial correlation is the
ad hoc estimator of Olsson, Drasgow and Dorans (1982),
\\r\_{xy}\\\sigma_y / \sum_k \phi(\tau_k)\\, rather than the two-step
maximum-likelihood estimator that
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
obtains through `lavaan`. Both are consistent, but on finite samples
they differ slightly (by up to about 0.006 in simulations), so a matrix
that mixes ordinal and continuous variables does not reproduce
[`qgraph::cor_auto()`](https://rdrr.io/pkg/qgraph/man/cor_auto.html)
exactly.

An ordered factor in a data frame enters as its integer level codes and
is ordinal whenever it has at most `ordinal_max_levels` levels.
Unordered factors, character and logical columns are dropped with a
message, and constant columns with a warning. A pair with no estimable
association, such as an item observed at a single level, has correlation
0.

## Conditions

`psychnet_bad_input` when `ordinal_max_levels` or `na_method` is
invalid, `data` has no numeric columns or fewer than two usable
variables, or a pair of variables has fewer than two joint observations;
`psychnet_constant_variable` (warning) when a constant or all-missing
column is dropped. A `psychnet_dropped_columns` message names
non-numeric columns that were dropped.

## References

Olsson, U. (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. *Psychometrika*, 44(4), 443-460.
[doi:10.1007/BF02296207](https://doi.org/10.1007/BF02296207)

Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial
correlation coefficient. *Psychometrika*, 47(3), 337-347.
[doi:10.1007/BF02294164](https://doi.org/10.1007/BF02294164)

## See also

[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md)
and the other estimators, which accept `cor_method = "auto"`.

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
