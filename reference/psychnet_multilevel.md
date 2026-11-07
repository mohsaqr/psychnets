# Within / between event-data networks

The object returned by
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
for nested event data with `standardize = FALSE`: a pair of Gaussian
graphical networks decomposing the actor-by-action frequency covariance
into a within-actor and a between-actor part.

## Usage

``` r
# S3 method for class 'psychnet_multilevel'
print(x, ...)

# S3 method for class 'psychnet_multilevel'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'psychnet_multilevel'
summary(object, ...)

# S3 method for class 'psychnet_multilevel'
plot(x, ...)
```

## Arguments

- x, object:

  A `psychnet_multilevel` object.

- ...:

  Ignored by `print`, `summary` and `as.data.frame`; passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
  by `plot`.

- row.names, optional:

  Ignored (S3 consistency).

## Value

`x`, invisibly (for `print`), after printing the numbers of actors and
occasions and each network's node and edge counts.

For `as.data.frame`, the two edge lists stacked with a leading `level`
column (`"within"` or `"between"`): columns `level`, `from`, `to`,
`weight`, one row per edge.

For `summary`, a data frame with one row per level (`"within"`,
`"between"`): columns `level`, `nodes`, `edges`, `mean_abs_weight` (0
for a network without edges) and `n` (the level's effective sample size:
pooled within-actor degrees of freedom, or the number of actors).

For `plot`, `x` invisibly, after drawing the within and between networks
side by side with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
(a Suggested package; a `psychnet_bad_input` error says so when it is
not installed).

## Examples

``` r
set.seed(1)
ev <- data.frame(
  Actor   = sample(paste0("s", 1:40), 1200, replace = TRUE),
  Session = sample(1:3, 1200, replace = TRUE),
  Action  = sample(c("read", "quiz", "note", "watch"), 1200, replace = TRUE))
# Actions drawn at random, so neither network is expected to have edges.
ml <- psychnet(ev, actor = "Actor", session = "Session", standardize = FALSE)
ml
#> <psychnet_multilevel> glasso networks (within + between)
#>   actors: 40   occasions: 120
#>   within : 4 nodes, 0 edges
#>   between: 4 nodes, 0 edges
summary(ml)
#>     level nodes edges mean_abs_weight  n
#> 1  within     4     0               0 80
#> 2 between     4     0               0 40
as.data.frame(ml)
#> [1] level  from   to     weight
#> <0 rows> (or 0-length row.names)
if (requireNamespace("cograph", quietly = TRUE)) plot(ml)
```
