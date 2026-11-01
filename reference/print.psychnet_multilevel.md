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
```

## Arguments

- x:

  A `psychnet_multilevel` object.

- ...:

  Ignored.

- row.names, optional:

  Ignored (S3 consistency).

## Value

`x`, invisibly (for `print`), after printing the numbers of actors and
occasions and each network's node and edge counts.

For `as.data.frame`, the two edge lists stacked with a leading `level`
column (`"within"` or `"between"`): columns `level`, `from`, `to`,
`weight`, one row per edge.

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
as.data.frame(ml)
#> [1] level  from   to     weight
#> <0 rows> (or 0-length row.names)
```
