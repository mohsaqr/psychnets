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
  estimator_args = list(),
  ...
)
```

## Arguments

- data:

  Data frame or matrix (rows = observations). Resampled exactly as given
  – factor columns and incomplete rows included – so each draw is handed
  to the estimator as the full data would be, and the estimator's own
  `na_method` and type handling apply per draw. For `method = "mgm"`
  node types are decided once on the full data and pinned for every
  draw.

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

- estimator_args:

  Named list of estimator arguments. Use this for names consumed by the
  bootstrap itself, such as an estimator's `threshold`. Values here
  override matching arguments in `...`.

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
# method = "pcor" keeps this example fast; the "glasso" default (and
# n_boot >= 1000) is what a real analysis should use.
bs <- net_boot(x, method = "pcor", n_boot = 50, cores = 1)
as.data.frame(bs)
#>    from to    observed         mean       lower      upper prop_nonzero
#> 1    V1 V2  0.31118543  0.298237891  0.15159446 0.40899218            1
#> 2    V1 V3  0.03185023  0.033387300 -0.14933331 0.15793256            1
#> 3    V2 V3  0.33426443  0.337131826  0.21222203 0.48649315            1
#> 4    V1 V4 -0.09321933 -0.099159232 -0.29858901 0.06030023            1
#> 5    V2 V4  0.04470097  0.034035886 -0.08595976 0.19118623            1
#> 6    V3 V4  0.32962441  0.338115084  0.14664746 0.50321968            1
#> 7    V1 V5  0.03435117  0.038882968 -0.11179003 0.19742845            1
#> 8    V2 V5 -0.01153139 -0.006380471 -0.16191974 0.12517682            1
#> 9    V3 V5 -0.01001689 -0.017513993 -0.18367088 0.20086310            1
#> 10   V4 V5  0.33034659  0.360416525  0.21937671 0.46332468            1
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
