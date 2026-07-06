# Tidy a network bootstrap

Tidy a network bootstrap

## Usage

``` r
# S3 method for class 'psychnet_bootstrap'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., significant = FALSE)
```

## Arguments

- x:

  A `psychnet_bootstrap` object.

- row.names, optional:

  Ignored (S3 consistency).

- ...:

  Unused.

- significant:

  If `TRUE`, return only the edges whose confidence interval excludes
  zero. Default `FALSE` (all edges).

## Value

The tidy `$edges` data frame (one row per edge, with its percentile
interval, inclusion proportion, and `significant` flag).
