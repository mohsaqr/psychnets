# Read a moderated mixed graphical model at a value of the moderator

`condition()` reads a moderated mixed graphical model at one value of
the moderator and returns that network as an ordinary `psychnet`;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`print()`](https://rdrr.io/r/base/print.html) read the moderated fit
itself. In a moderated network model, fitted by
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
with `moderators`, the strength of an edge between two variables can
depend on the value of a third variable, the moderator, so each pairwise
association has a main effect and an interaction with the moderator.

The network at a given moderator value, for example in a treatment group
compared with a control group, or at low compared with high stress,
shows how the conditional associations among the other variables differ
across the moderator. For each pair, `condition()` takes the main effect
plus the interaction evaluated at `value`, aggregated over the two
nodewise regressions and combined by the AND or OR rule. Pairs whose
main effect or interaction did not survive the rule at fitting time
contribute nothing. The moderator itself carries no edges, and edge
weights are magnitudes (non-negative), as in
[`mgm::condition()`](https://rdrr.io/pkg/mgm/man/condition.html).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
moderated fit lists the size of the main effect and of the moderation
effect of every pair with a nonzero effect, or, with `value`, the edges
of the conditioned network.

`value` is read on the moderator's own scale: a level of a categorical
moderator as coded in the data (`1` or `2` for a moderator coded 1/2,
`"treat"` for a factor level `"treat"`), or a raw value of a continuous
moderator, so that `value = mean(moderator)` gives the network at the
average moderator. The conditioned network is computed from the fit
rather than estimated, so
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports its certificate as missing.

`condition()` returns an undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

To fit the moderated model, use
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
with `moderators`.

## Usage

``` r
condition(object, value, rule = NULL)

# S3 method for class 'psychnet_moderated'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., value = NULL)

# S3 method for class 'psychnet_moderated'
print(x, ...)
```

## Arguments

- object:

  A `psychnet_moderated` fit from
  [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  with `moderators`.

- value:

  A single moderator value: a level of a categorical moderator (number,
  label, factor or logical, as coded in the data), or a finite number on
  the raw scale of a continuous one. For
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), `NULL`
  (default) gives the table of main and moderation effects.

- rule:

  `"AND"` or `"OR"` combination of the two conditioned nodewise
  estimates of each pair. Default `NULL` uses the rule of the fit.

- x:

  A `psychnet_moderated` fit.

- row.names, optional:

  Ignored.

- ...:

  Unused.

## Value

`condition()` returns an undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md));
the value conditioned on is stored as `mod_value`.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
moderated fit returns, with `value = NULL`, one row per pair of nodes
joined by a main effect or a moderation effect, with columns `from`,
`to`, `main_effect` (magnitude of the pairwise main effect: the mean
absolute nodewise estimate, zero when the pair does not survive the
rule), `interaction` (magnitude of the moderator-by-pair interaction,
`NA` for a pair that includes the moderator), `moderated`
(`interaction > 0`, `NA` for a pair that includes the moderator) and
`involves_moderator`. With a `value`, it returns the edge list (`from`,
`to`, `weight`) of `condition(x, value)`.

[`print()`](https://rdrr.io/r/base/print.html) shows the number of
nodes, the moderator and the fit settings, and returns `x` invisibly.

## Functions

- `as.data.frame(psychnet_moderated)`: Main and moderation effects of
  every pair, or (with `value`) the edge list of the conditioned
  network.

- `print(psychnet_moderated)`: Print the nodes, the moderator and the
  fit settings.

## Conditions

`psychnet_bad_input` when `object` is not a moderated fit, `rule` is not
`"AND"` or `"OR"`, or `value` is missing, not of length one, not a level
of a categorical moderator, or not a finite number for a continuous one.

## See also

[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
for fitting the moderated model;
[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the conditioned network.

## Examples

``` r
set.seed(1)
x1 <- stats::rnorm(400); x2 <- stats::rnorm(400)
mod <- rep(1:2, each = 200)                 # moderator coded 1/2
y <- x1 * (mod == 2) + stats::rnorm(400)   # x1-y edge only when mod == 2
d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = "mod")
fit
#> <psychnet_moderated> 4 nodes, moderator = mod
#>   gamma = 0.25, rule = AND, threshold = LW
#>   moderator levels: 1, 2
#>   as.data.frame(fit) for main and moderation effects; condition(fit, value) for the network at a moderator value.
as.data.frame(fit)            # main and moderation effects per edge
#>   from to main_effect interaction moderated involves_moderator
#> 1   x1  y           0   0.5665183      TRUE              FALSE
condition(fit, value = 1)     # no x1-y edge
#> <psychnet> mgm_moderated network
#>   nodes: 4   edges: 0   (undirected)
condition(fit, value = 2)     # the x1-y edge appears
#> <psychnet> mgm_moderated network
#>   nodes: 4   edges: 1   (undirected)
as.data.frame(fit, value = 2) # the same network as an edge list
#>   from to   weight
#> 1   x1  y 0.646064
```
