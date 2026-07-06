# Bootstrap a psychometric network

Resamples observations with replacement, re-estimates the network on
each resample, and summarizes the sampling distribution of every edge
weight and node centrality (mean, percentile confidence interval, and
edge inclusion proportion). An edge is flagged `significant` when its
percentile interval excludes zero. The raw per-resample draws are stored
on the returned object for use by
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md).

## Usage

``` r
net_boot(
  data,
  method = "glasso",
  n_boot = 1000L,
  ci = 0.95,
  measures = c("strength", "expected_influence"),
  centrality_fn = NULL,
  predictability = FALSE,
  threshold = FALSE,
  diff_test = FALSE,
  p_adjust = "none",
  labels = NULL,
  cores = NULL,
  engine = NULL,
  ...
)
```

## Arguments

- data:

  Numeric data frame or matrix (rows = observations).

- method:

  Estimator (see
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)).
  Default `"glasso"`.

- n_boot:

  Number of bootstrap resamples. Default 1000.

- ci:

  Confidence level for percentile intervals. Default 0.95.

- measures:

  Centrality measures to bootstrap. Defaults to the two recommended for
  psychometric networks (`"strength"`, `"expected_influence"`);
  `"betweenness"`/`"closeness"` and custom measures (via
  `centrality_fn`) are also accepted. See
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).

- centrality_fn:

  Optional function supplying any non-built-in `measures` (see
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)).

- predictability:

  Logical; if `TRUE` and the estimator returns a precision matrix (GGM
  family), bootstrap node predictability (R^2) and report its interval.
  Default `FALSE`.

- threshold:

  Logical; if `TRUE`, also return the observed network with every edge
  whose bootstrap interval includes zero set to zero (`$thresholded`).
  Default `FALSE`.

- diff_test:

  Logical; if `TRUE`, also return two-sided bootstrap difference p-value
  matrices for edges (`$edge_diff_p`, `NULL` past 500 edges) and for
  each centrality measure (`$centrality_diff_p`). Default `FALSE`.

- p_adjust:

  Multiple-comparison adjustment applied to the difference p-value
  matrices (any [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html)
  method). Default `"none"`.

- labels:

  Optional node labels.

- cores:

  Number of CPU cores for the resample loop. `NULL` (default) uses two
  thirds of the detected cores; `1` forces a serial run. Parallelism
  uses forking
  ([`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)) and
  falls back to serial on Windows. Because every resample index is drawn
  in the parent process before any fitting, the result is identical for
  any number of cores and reproducible from
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- engine:

  Optional estimator engine forwarded to each resample fit (e.g.
  `"base"`/`"glasso"` for glasso, `"base"`/`"glmnet"` for ising/mgm).
  `NULL` (default) uses the estimator's own default.

- ...:

  Passed to the estimator.

## Value

An object of class `psychnet_bootstrap`: tidy `$edges` (with a
`significant` flag) and `$centrality` data frames, the observed network
in `$observed`, raw resample draws in `$edge_boot`, `$str_boot`,
`$ei_boot`, and the general `$centrality_boot` (named list, one matrix
per measure). Optional `$predictability`, `$thresholded`,
`$edge_diff_p`, `$centrality_diff_p`, plus
`$lambda_path`/`$lambda_selected` when the estimator reports them.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
bs <- net_boot(x, n_boot = 50, cores = 1)
as.data.frame(bs)
#>    from to  observed         mean       lower      upper prop_nonzero
#> 1    V1 V2 0.2647959  0.216657987  0.06460790 0.34930531         0.98
#> 2    V1 V3 0.0000000  0.019208250  0.00000000 0.11096685         0.24
#> 3    V2 V3 0.3114876  0.270731621  0.14521228 0.39115437         1.00
#> 4    V1 V4 0.0000000 -0.024222579 -0.18536844 0.00000000         0.18
#> 5    V2 V4 0.0000000  0.017061134  0.00000000 0.11199844         0.22
#> 6    V3 V4 0.2936317  0.262160092  0.06928587 0.42053424         0.98
#> 7    V1 V5 0.0000000  0.002383375 -0.08504584 0.09770776         0.16
#> 8    V2 V5 0.0000000  0.004825769  0.00000000 0.07379984         0.12
#> 9    V3 V5 0.0000000  0.015522783  0.00000000 0.13699210         0.24
#> 10   V4 V5 0.2807738  0.265568744  0.13497982 0.38097285         1.00
#>    significant
#> 1         TRUE
#> 2        FALSE
#> 3         TRUE
#> 4        FALSE
#> 5        FALSE
#> 6         TRUE
#> 7        FALSE
#> 8        FALSE
#> 9        FALSE
#> 10        TRUE
```
