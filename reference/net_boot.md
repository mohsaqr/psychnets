# Bootstrap the edges and centralities of a network

`net_boot()` bootstraps a network to show how accurately its edge
weights and node centralities are estimated. The nonparametric bootstrap
approximates the sampling variation of the network, whose edge weights
and centralities would differ in another sample from the same
population, by re-estimating it on samples drawn with replacement from
the observations; the spread of the bootstrapped values shows which
edges, and which differences between them, can be interpreted (Epskamp,
Borsboom & Fried, 2018).

The function resamples the rows of `data`, re-estimates the network with
the chosen estimator on every resample, and summarizes each edge weight
and each node centrality by its bootstrap mean and percentile confidence
interval. For each edge it also reports the proportion of resamples in
which the edge was estimated as nonzero, and flags it as `significant`
when its interval excludes zero. Options add intervals for node
predictability (`predictability`) and a version of the network that
keeps only the edges whose interval excludes zero (`threshold`).

It returns these intervals as tables read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
keeps the bootstrapped values on the object, so that
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
can test differences between edges, or between node centralities,
without resampling again.

The bootstrap describes the accuracy of estimates within one sample. For
the stability of the order of centralities when cases are dropped, use
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md);
for the stability of the edge weights, use
[`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
and
[`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md);
to test whether two groups have different networks, use
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md).

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

# S3 method for class 'psychnet_bootstrap'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("edges", "centrality", "predictability", "thresholded"),
  significant = FALSE
)

# S3 method for class 'psychnet_bootstrap'
print(x, ...)

# S3 method for class 'psychnet_bootstrap'
plot(
  x,
  type = c("edges", "centrality", "edge_diff", "centrality_diff", "predictability"),
  measure = NULL,
  ...
)
```

## Arguments

- data:

  A data frame or matrix with one row per observation, or a
  [psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  estimated from a data set, in which case each group is bootstrapped
  with the group's own estimator settings.

- method:

  The estimator, as in
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md).
  Default `"glasso"`.

- n_boot:

  Number of bootstrap resamples. Default `1000`.

- ci:

  Confidence level of the percentile intervals, a number in (0, 1).
  Default `0.95`.

- measures:

  Centrality measures to bootstrap. Default
  `c("strength", "expected_influence")`; any measure accepted by
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
  can be added.

- centrality_fn:

  Optional function supplying measures that are not built in (see
  [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)).

- predictability:

  If `TRUE`, also bootstrap node predictability (R^2) when the estimator
  returns a precision matrix (the Gaussian graphical models). Default
  `FALSE`.

- threshold:

  If `TRUE`, also keep the observed network with every edge whose
  interval includes zero set to zero. Default `FALSE`.

- diff_test:

  If `TRUE`, also store matrices of pairwise difference p-values on the
  object: for the edges (up to 500 edges) and for each centrality
  measure. Default `FALSE`. The tidy version of these tests is
  [`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md).

- p_adjust:

  Multiplicity adjustment for the matrices stored by `diff_test`, any
  method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"none"`.

- labels:

  Optional node labels, one per node.

- cores:

  Number of CPU cores. `NULL` (default) uses two thirds of the detected
  cores; `1` runs serially. Runs are serial on Windows.

- engine:

  Optional solver for every refit: `"base"` or `"glasso"` for
  `method = "glasso"` or `"huge"`, `"base"` or `"glmnet"` for `"ising"`
  or `"mgm"`. `NULL` (default) keeps the estimator's default.

- estimator_args:

  Named list of arguments for the estimator. Use it for an estimator
  argument that shares a name with an argument of `net_boot()`, e.g.
  `estimator_args = list(threshold = TRUE)`. Values here override those
  in `...`.

- ...:

  Further arguments for the estimator (`net_boot()`); ignored by the
  methods.

- x:

  A `psychnet_bootstrap` object.

- row.names, optional:

  Ignored.

- what:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):
  `"edges"` (default), `"centrality"`, `"predictability"` (needs
  `predictability = TRUE`) or `"thresholded"` (needs
  `threshold = TRUE`).

- significant:

  For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  with `what = "edges"`: if `TRUE`, return only the edges whose interval
  excludes zero. Default `FALSE`.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
  `"edges"` (default; edge weights with their intervals, sorted by
  observed weight), `"centrality"` (one panel per measure),
  `"edge_diff"` or `"centrality_diff"` (the matrix of pairwise
  difference tests for the edges or for one measure), or
  `"predictability"` (needs `predictability = TRUE`).

- measure:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
  `type = "centrality"` or `"centrality_diff"`: the measure(s) to draw.
  Default: all bootstrapped measures, or the first one for
  `"centrality_diff"`.

## Value

`net_boot()` returns an object of class `psychnet_bootstrap`. Its tables
are read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

- `as.data.frame(x)` (`what = "edges"`): one row per edge, with columns
  `from`, `to`, `observed`, `mean` (bootstrap mean), `lower`, `upper`
  (percentile interval), `prop_nonzero` (share of resamples in which the
  edge is nonzero) and `significant` (the interval excludes zero). For a
  directed estimator each ordered pair is a row.

- `what = "centrality"`: one row per node, with `node` and, for each
  measure, its observed value and the columns `<measure>_lower` and
  `<measure>_upper`.

- `what = "predictability"`: one row per node, with `node`, `value` (the
  observed R^2), `lower` and `upper`.

- `what = "thresholded"`: one row per edge, with `from`, `to`,
  `observed`, `weight` (the observed weight, or 0 when its interval
  includes zero) and `retained`.

The object also stores the observed network and the bootstrapped value
of every edge and centrality in every resample, which
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
and the difference plots use.
[`print()`](https://rdrr.io/r/base/print.html) shows the estimator, the
number of resamples, the interval level and the numbers of edges,
significant edges and nodes, and returns `x` invisibly;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns `x`
invisibly. For a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md),
`net_boot()` returns one result per group (see
[psychnet_result_group](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md)).

## Details

Each resample is passed to the estimator exactly as the full data would
be, including incomplete rows and factor columns, so the estimator's own
handling of missing values and variable types applies to every resample.
For `method = "mgm"` the node types are decided once on the full data
and kept for every resample. A resample on which a variable has no
variation keeps that variable as an isolated node.

All resample indices are drawn before any network is fitted, so the
result is reproducible with
[`set.seed()`](https://rdrr.io/r/base/Random.html) and does not depend
on `cores`.

For an estimator whose edge weights are non-negative by construction
(`method = "relimp"`, whose weights are shares of explained variance), a
percentile interval excludes zero almost surely: every edge is then
flagged `significant` and `threshold = TRUE` removes nothing. Pairwise
difference tests with
[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
remain informative.

## Functions

- `as.data.frame(psychnet_bootstrap)`: The edge intervals (default), the
  centrality or predictability intervals, or the thresholded network,
  chosen by `what`.

- `print(psychnet_bootstrap)`: The estimator, the number of resamples,
  the interval level and the edge counts.

- `plot(psychnet_bootstrap)`: Edge, centrality or predictability
  intervals, or a matrix of pairwise difference tests (see `type`).
  Edges whose interval excludes zero are filled points on solid
  segments, the others open points on dashed segments.

## Conditions

`psychnet_bad_input` for an invalid argument, for a group that was not
estimated from a data set, and from
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) or
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) when the
requested table was not computed; `psychnet_estimation_failed` when the
network cannot be estimated on the full data (with the estimator's
message); `psychnet_resample_failed` (warning) when some resamples fail,
naming the first error, and as an error when fewer than two succeed;
`psychnet_constant_variable` (warning) when a constant column is
dropped; `psychnet_ignored_argument` (warning) when
`predictability = TRUE` but the estimator returns no precision matrix;
`psychnet_uninformative_test` (warning) for `threshold = TRUE` or
`as.data.frame(significant = TRUE)` with `method = "relimp"`.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods*, 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## See also

[`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
for pairwise differences;
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
for the measures;
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md),
[`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
and
[`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
for case-dropping and split-half stability;
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
for the estimators.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.4^abs(outer(1:5, 1:5, "-")))
colnames(x) <- paste0("V", 1:5)
# method = "pcor" and n_boot = 50 keep this example fast; a real analysis
# uses the "glasso" default and n_boot >= 1000.
bs <- net_boot(x, method = "pcor", n_boot = 50, cores = 1)
bs
#> <psychnet_bootstrap> pcor, 50 resamples, 95% CI
#>   10 edges (4 significant), 5 nodes, measures: strength, expected_influence
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
as.data.frame(bs, significant = TRUE)
#>   from to  observed      mean     lower     upper prop_nonzero significant
#> 1   V1 V2 0.3111854 0.2982379 0.1515945 0.4089922            1        TRUE
#> 2   V2 V3 0.3342644 0.3371318 0.2122220 0.4864931            1        TRUE
#> 3   V3 V4 0.3296244 0.3381151 0.1466475 0.5032197            1        TRUE
#> 4   V4 V5 0.3303466 0.3604165 0.2193767 0.4633247            1        TRUE
as.data.frame(bs, what = "centrality")
#>   node  strength strength_lower strength_upper expected_influence
#> 1   V1 0.4706062      0.3089894       0.907735          0.2841675
#> 2   V2 0.7016822      0.5439371       1.001855          0.6786194
#> 3   V3 0.7057560      0.6049353       1.116255          0.6857222
#> 4   V4 0.7978913      0.6316038       1.182999          0.6114526
#> 5   V5 0.3862460      0.4024386       0.854440          0.3431495
#>   expected_influence_lower expected_influence_upper
#> 1               0.05695078                0.4525569
#> 2               0.48748685                0.8893284
#> 3               0.48989974                0.9386692
#> 4               0.48030709                0.9109487
#> 5               0.08638032                0.6981966
plot(bs)                       # edge weights with their intervals

plot(bs, type = "centrality")  # centralities with their intervals
```
