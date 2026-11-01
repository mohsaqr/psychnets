# Print a moderated MGM fit

Print a moderated MGM fit

## Usage

``` r
# S3 method for class 'psychnet_moderated'
print(x, ...)
```

## Arguments

- x:

  A `psychnet_moderated` object.

- ...:

  Unused.

## Value

`x`, invisibly. Prints the number of nodes, the moderator, and the fit
settings.

## Examples

``` r
set.seed(1)
x1 <- stats::rnorm(400); x2 <- stats::rnorm(400)
mod <- rep(0:1, each = 200)
y <- x1 * (mod == 1) + stats::rnorm(400)
d <- data.frame(x1 = x1, x2 = x2, y = y, mod = mod)
print(mgm_fit(d, types = c("g", "g", "g", "c"), moderators = 4))
#> <psychnet_moderated> 4 nodes, moderator = mod
#>   gamma = 0.25, rule = AND, threshold = LW
#>   condition(fit, value) to extract the network at a moderator value.
```
