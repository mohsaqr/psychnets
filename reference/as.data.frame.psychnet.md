# Tidy edge list for a psychnet network

Tidy edge list for a psychnet network

## Usage

``` r
# S3 method for class 'psychnet'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., include_zero = FALSE)
```

## Arguments

- x:

  A `psychnet` object.

- row.names, optional:

  Ignored (for S3 consistency).

- ...:

  Unused.

- include_zero:

  If TRUE, keep zero-weight (absent) edges. Default FALSE.

## Value

A one-row-per-edge `data.frame` with columns `from`, `to`, `weight`.
