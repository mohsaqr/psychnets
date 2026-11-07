# Center and scale variables within persons

Removes each person's own mean from every variable and, by default,
divides by that person's own standard deviation, so that every network
estimated from the result describes *within-person* associations: how a
person's deviations from their own typical level go together. This is
the within-person standardization of `mlVAR::mlVAR(scaleWithin = TRUE)`
(centering per person, then division by the person's SD), offered as a
standalone step so it can precede any analysis. The same scaling is
built into
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
as `within =`, which also lets
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
pair two occasions by person.

## Usage

``` r
scale_within(data, id, vars = NULL, center = TRUE, scale = TRUE)
```

## Arguments

- data:

  A data frame with one row per observation (person x occasion).

- id:

  Name of the column identifying persons. Missing identifiers are an
  error.

- vars:

  Names of the columns to center/scale. Default `NULL`: every numeric
  column except `id`. Name them explicitly when the data carry other
  numeric columns (an occasion index, a covariate) that must stay as
  they are.

- center:

  Subtract each person's mean. Default `TRUE`.

- scale:

  Divide by each person's standard deviation. Default `TRUE`.

## Value

A data frame with the same rows and columns as `data`, the `vars`
columns replaced by their within-person centered (and scaled) values and
every other column unchanged.

Errors of class `psychnet_bad_input` for a missing `id` column, missing
identifiers, or `vars` that are not numeric columns of `data`.

## Details

Means and standard deviations use each person's available values (the SD
has denominator n - 1, as
[`stats::sd()`](https://rdrr.io/r/stats/sd.html)). A value is set to
`NA` when its person has fewer than two observed values of that
variable, or (with `scale = TRUE`) when the person's values of it are
constant: there is no within-person variation to describe. Both cases
raise a `psychnet_constant_within` warning with the count.

**Two occasions per person.** With exactly two observations, each
person's centered deviations at the two occasions are mirror images
(`x1 - m = -(x2 - m)`), so the networks of the two occasions are
identical by construction and a comparison between them cannot detect
anything; after scaling, every value is \\\pm 1/\sqrt{2}\\. When no
person has more than two observations a `psychnet_mirror_deviations`
warning says so. Within-person comparisons of occasions need three or
more occasions per person.

## References

Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018). The
Gaussian graphical model in cross-sectional and time-series data.
*Multivariate Behavioral Research*, 53(4), 453-480.
[doi:10.1080/00273171.2018.1454823](https://doi.org/10.1080/00273171.2018.1454823)

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
