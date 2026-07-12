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
```

## Arguments

- object:

  A `psychnet_moderated` object from `mgm_fit(..., moderators=)`.

- value:

  Moderator value to condition on (e.g. `0` or `1` for a binary
  moderator, or any numeric for a continuous one).

- rule:

  Symmetrization rule; defaults to the rule used at fit time.

## Value

A `psychnet` network object (the moderator node carries no edges).

## Examples

``` r
set.seed(1)
x1 <- stats::rnorm(400); x2 <- stats::rnorm(400)
mod <- rep(0:1, each = 200)
y <- x1 * (mod == 1) + stats::rnorm(400)   # x1-y edge only when mod == 1
d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
# `moderators=` is the one estimator that needs glmnet (a Suggested package).
if (requireNamespace("glmnet", quietly = TRUE)) {
  fit <- mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4)
  condition(fit, value = 1)
}
#> <psychnet> mgm_moderated network
#>   nodes: 4   edges: 1   (undirected)
```
