# Tables of a Network Comparison Test

[`summary()`](https://rdrr.io/r/base/summary.html) gives the tests,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) the
edges.

## Usage

``` r
# S3 method for class 'psychnet_nct'
summary(object, alpha = 0.05, ...)

# S3 method for class 'psychnet_nct'
as.data.frame(x, row.names = NULL, optional = FALSE, alpha = 0.05, ...)
```

## Arguments

- alpha:

  Significance level for the `significant` column. Default 0.05.

- ...:

  Unused.

- x, object:

  A `psychnet_nct` object from
  [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md).

- row.names, optional:

  Ignored (base generic compatibility).

## Value

[`summary()`](https://rdrr.io/r/base/summary.html): a data frame with
one row per global test, columns `test` (`"network structure"`;
`"global strength"`, or `"global expected influence"` when
`abs = FALSE`), `statistic` (`"M"`, `"S"`), `observed`, `p_value` and
`significant` (`p_value < alpha`).

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a data
frame with one row per node pair, columns `from`, `to`, `weight_1`,
`weight_2` (the edge in the first and second group's network; 0/1 when
`weighted = FALSE`), `difference` (absolute difference, statistic E),
`p_value` (permutation p, adjusted by `p_adjust`) and `significant`.

## Examples

``` r
set.seed(1)
a <- matrix(stats::rnorm(100 * 4), 100, 4)
b <- matrix(stats::rnorm(100 * 4), 100, 4)
colnames(a) <- colnames(b) <- paste0("V", 1:4)
fit <- net_compare(a, b, iter = 10)
summary(fit)
#>                test statistic observed p_value significant
#> 1 network structure         M        0       1       FALSE
#> 2   global strength         S        0       1       FALSE
as.data.frame(fit)
#>   from to weight_1 weight_2 difference p_value significant
#> 1   V1 V2        0        0          0       1       FALSE
#> 2   V1 V3        0        0          0       1       FALSE
#> 3   V2 V3        0        0          0       1       FALSE
#> 4   V1 V4        0        0          0       1       FALSE
#> 5   V2 V4        0        0          0       1       FALSE
#> 6   V3 V4        0        0          0       1       FALSE
```
