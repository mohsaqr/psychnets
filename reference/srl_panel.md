# Synthetic weekly self-regulated-learning panel

A synthetic panel of 100 students measured in four consecutive weeks on
five self-regulated-learning indicators – a small, fully shareable data
set for the panel graphical VAR
([`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md)).

## Usage

``` r
srl_panel
```

## Format

A data frame with 400 rows (100 students x 4 weeks) and 7 columns:

- `student`:

  Student identifier, `"S001"` to `"S100"`.

- `week`:

  Week of measurement, 1 to 4.

- `efficacy`:

  Self-efficacy.

- `planning`:

  Planning.

- `effort`:

  Effort regulation.

- `enjoyment`:

  Enjoyment.

- `anxiety`:

  Anxiety.

## Source

Simulated; structure derived from the weekly means of the `esm_srl` data
in the idiographic package.

## Details

The data were simulated, once, from a stationary random-intercept VAR(1)
(the panel graphical VAR model) whose structure was taken from the
weekly means of real experience-sampling data (the anonymized `esm_srl`
data of the idiographic package): the contemporaneous partial
correlations, the temporal pattern (autoregressive effects and an effort
-\> efficacy lagged effect, scaled to 0.6 of the estimated values so
trait and state differences are separable in four waves), and the
correlations among students' average levels, with 60% of each variable's
variance between students. Values are on a 0-100 scale (one of the 2,000
values was bounded to the scale) and rounded to one decimal. No real
student is represented.

## Examples

``` r
head(srl_panel)
#>   student week efficacy planning effort enjoyment anxiety
#> 1    S001    1     38.7     54.1   52.1      53.1    44.2
#> 2    S001    2     32.9     55.4   50.7      62.2    48.8
#> 3    S001    3     34.5     40.9   52.9      50.3    36.9
#> 4    S001    4     28.6     44.5   41.0      42.4    39.6
#> 5    S002    1     45.0     54.7   48.6      48.1    38.4
#> 6    S002    2     62.8     67.2   63.7      47.5    44.8
fit <- panel_gvar(srl_panel,
                  vars = c("efficacy", "planning", "effort", "enjoyment",
                           "anxiety"),
                  id = "student", wave = "week")
summary(fit)
#>           network directed nodes edges significant mean_abs_weight
#> 1        temporal     TRUE     5    20           1       0.1046718
#> 2 contemporaneous    FALSE     5    10           5       0.1877681
#> 3         between    FALSE     5    10           2       0.1831649
```
