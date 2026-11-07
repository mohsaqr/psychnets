# Condition a moderated network at a moderator value

Extracts the effective pairwise network implied by a moderated MGM
([`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
with `moderators`) at a given value of the moderator, mirroring
[`mgm::condition()`](https://rdrr.io/pkg/mgm/man/condition.html): it
applies the AND-rule pre-filter, absorbs the moderator value into the
main-effect coefficients, and re-aggregates the pairwise edges.

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

  A `psychnet_moderated` object from `mgm_fit(..., moderators=)`.

- value:

  A single moderator value to condition on: a level of a categorical
  moderator (number, label, factor or logical, as coded in the data), or
  a finite number on the raw scale of a continuous one.

- rule:

  Symmetrization rule; defaults to the rule used at fit time.

- x:

  A `psychnet_moderated` object.

- row.names, optional:

  Standard
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  arguments; ignored.

- ...:

  Unused.

## Value

`condition()` returns a `psychnet` network object (the moderator node
carries no edges), with `$mod_value` the value conditioned on.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
moderated fit returns, with `value = NULL`, one row per pair of nodes
joined by a main effect or a moderation effect, with columns `from`,
`to`, `main_effect` (aggregated magnitude of the pairwise main effect,
mgm convention: mean absolute nodewise estimate, zero when the pair does
not survive `rule`), `interaction` (aggregated magnitude of the
moderator x pair three-way interaction, `NA` for a pair that includes
the moderator itself), `moderated` (logical: `interaction > 0`) and
`involves_moderator`. With a `value`, it returns the tidy edge list
(`from`, `to`, `weight`) of `condition(x, value)`.

[`print()`](https://rdrr.io/r/base/print.html) on the moderated fit
shows the number of nodes, the moderator and the fit settings, and
returns `x` invisibly.

Errors of class `psychnet_bad_input` when `value` is missing, not length
one, not a level of a categorical moderator, or not a finite number for
a continuous one.

## Details

`value` is always read on the moderator's **own** scale. For a
categorical moderator it is one of the levels observed in the data, in
the coding the data used: `1`/`2` for a moderator coded 1/2, `"treat"`
for a factor level `"treat"` (the internal 0-based level code is never
exposed). For a continuous moderator it is a value on the raw scale; the
moderator is standardized internally (mgm's `scale = TRUE`), and the
value is mapped to `(value - mean) / sd` before conditioning, so
`value = mean(moderator)` gives the network at the average moderator.
([`mgm::condition()`](https://rdrr.io/pkg/mgm/man/condition.html)
instead takes the value on its internal standardized scale:
`condition(fit, v)` here equals
`mgm::condition(ref, list(k = (v - mean) / sd))`.)

## Functions

- `as.data.frame(psychnet_moderated)`: Tidy table of the moderated fit:
  per-edge main and moderation effects, or (with `value`) the
  conditioned edge list.

- `print(psychnet_moderated)`: The moderated fit: nodes, moderator and
  settings.

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
as.data.frame(fit, value = 2) # the same network as a tidy edge list
#>   from to   weight
#> 1   x1  y 0.646064
```
