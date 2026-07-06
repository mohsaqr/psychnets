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

`x`, invisibly (for `print`).

For `as.data.frame`, the two edge lists stacked with a `level` column.
