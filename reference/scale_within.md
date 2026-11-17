# Center and scale variables within persons

`scale_within()` centers and scales variables within persons: it
subtracts each person's mean from each variable and, by default, divides
by that person's standard deviation, as
`mlVAR::mlVAR(scaleWithin = TRUE)` does. Removing each person's own mean
leaves only the within-person part of the associations, so a network
estimated from the result describes how a person's deviations from their
own typical level go together (Epskamp et al., 2018).

When each person is measured on several occasions, the association
between two variables can arise between persons (people who are high on
one tend to be high on the other) or within persons (on occasions when a
person is higher than usual on one, they tend to be higher than usual on
the other). Means and standard deviations use each person's available
values (the standard deviation with denominator n - 1, as
[`stats::sd()`](https://rdrr.io/r/stats/sd.html)). A value is set to
`NA` when its person has fewer than two observed values of the variable
or, with `scale = TRUE`, when the person's values of it are constant,
because there is then no within-person variation to describe.

The result is `data` with the `vars` columns replaced by their
within-person centered (and scaled) values.

The same scaling is built into
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
as `within =`, which also lets
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
pair two occasions by person. For panel data in which within-person
effects from one wave to the next are of interest, use
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md).

## Usage

``` r
scale_within(data, id, vars = NULL, center = TRUE, scale = TRUE)
```

## Arguments

- data:

  A data frame with one row per observation (person by occasion).

- id:

  Name of the column identifying persons; it may not contain missing
  values.

- vars:

  Names of the columns to center and scale. Default `NULL`: every
  numeric column except `id`. Name them explicitly when the data contain
  other numeric columns, such as an occasion index or a covariate, that
  must stay as they are.

- center:

  Subtract each person's mean. Default `TRUE`.

- scale:

  Divide by each person's standard deviation. Default `TRUE`.

## Value

A data frame with the same rows and columns as `data`, in which the
`vars` columns hold the within-person centered (and scaled) values and
every other column is unchanged.

## Details

With exactly two observations per person, each person's centered values
at the two occasions are mirror images (\\x_1 - m = -(x_2 - m)\\), so
the networks of the two occasions are identical by construction and a
comparison between them cannot detect anything; after scaling, every
value is \\\pm 1/\sqrt{2}\\. Within-person comparisons of occasions need
three or more occasions per person.

## Conditions

`psychnet_bad_input` when `data` is not a data frame, `id` does not name
one column or has missing values, `vars` names a missing or non-numeric
column or includes `id`, no numeric column is left to scale, or `center`
or `scale` is not `TRUE` or `FALSE`; `psychnet_mirror_deviations`
(warning) when no person has more than two observations;
`psychnet_constant_within` (warning) when values are set to `NA` for
lack of within-person variation, with their count.

## References

Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018). The
Gaussian graphical model in cross-sectional and time-series data.
*Multivariate Behavioral Research*, 53(4), 453-480.
[doi:10.1080/00273171.2018.1454823](https://doi.org/10.1080/00273171.2018.1454823)

## See also

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `within`;
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md);
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md).

## Examples

``` r
v <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
sw <- scale_within(srl_panel, id = "student", vars = v)
head(sw)
#>   student week   efficacy   planning     effort  enjoyment    anxiety
#> 1    S001    1  1.2036432  0.7538154  0.5293811  0.1346210  0.3484629
#> 2    S001    2 -0.1856365  0.9361336  0.2760021  1.2483039  1.2267804
#> 3    S001    3  0.1976131 -1.0974150  0.6741690 -0.2080506 -1.0453887
#> 4    S001    4 -1.2156197 -0.5925340 -1.4795522 -1.1748743 -0.5298546
#> 5    S002    1 -0.9608572 -0.4997138 -0.9140280 -0.0431398 -0.8965450
#> 6    S002    2  1.3348821  1.0618918  1.1306864 -0.2156990  0.4944581

# The same scaling inside the front door: one within-person network per
# week; compare two weeks with net_compare(nets, c(1, 4)), which pairs the
# students automatically (see ?psychnet).
psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
         within = "student")
#> <psychnet_group> 4 networks by week (method: glasso)
#>  group nodes edges   n
#>      1     5     6 100
#>      2     5     6 100
#>      3     5     7 100
#>      4     5     7 100
```
