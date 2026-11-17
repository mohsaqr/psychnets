# Estimate a psychometric network with any of the package's estimators

`psychnet()` estimates a psychometric network with any of the
cross-sectional estimators, selected by `method`. A psychometric network
represents variables, such as questionnaire items, symptoms or
behaviours, as nodes and their conditional or marginal associations as
edges.

The estimators differ in what an edge means: partial-correlation
networks describe the association between two variables after
controlling for all others, regularized networks retain only the
associations the data support, and Ising and mixed graphical models
apply the same idea to binary and mixed data. Besides `method`, the
arguments select the variables (`vars`), estimate one network per group
(`group`), scale each variable within persons first (`within`), or treat
`data` as an event log (`actor`, `action`, `session`, `time`), which is
converted to action frequencies per occasion before a network is
estimated.

`psychnet()` returns a fitted network; with `group`, one network per
group, and for some event data a pair of within-actor and between-actor
networks.

Each method can also be called directly, which exposes its full set of
options. Panel data measured in a few waves are modelled by
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md),
which `psychnet()` does not route. The method-specific pages describe
each model:
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md),
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md),
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md),
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md),
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md),
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md),
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md),
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md).

## Usage

``` r
psychnet(
  data,
  method = c("glasso", "cor", "pcor", "ising", "mgm", "huge", "ggm", "tmfg", "logo",
    "relimp", "ising_sampler"),
  threshold = 0,
  gamma = NULL,
  labels = NULL,
  vars = NULL,
  group = NULL,
  within = NULL,
  source = c("auto", "data", "eventdata"),
  actor = NULL,
  action = "Action",
  session = NULL,
  time = NULL,
  compute_sessions = TRUE,
  time_threshold = 900,
  id = NULL,
  standardize = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame or matrix with one row per observation and one column per
  variable, or a long event log (one row per event).

- method:

  The estimator, a single name: `"glasso"` (default), `"ggm"`, `"huge"`,
  `"tmfg"`, `"logo"`, `"relimp"`, `"cor"`, `"pcor"`, `"ising"`,
  `"ising_sampler"` or `"mgm"`, or a `qgraph`/`bootnet` alias. Event
  data accept `"glasso"`, `"pcor"`, `"ggm"`, `"huge"`, `"tmfg"` and
  `"logo"`.

- threshold:

  Absolute edge weight below which edges are set to zero, for `"cor"`,
  `"pcor"`, `"glasso"`, `"huge"`, `"ggm"` and `"logo"`: a single number
  \>= 0 (default `0`). For `"mgm"`, give its rule (`"LW"`, `"HW"` or
  `"none"`); the default `0` leaves the
  [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  default in place.

- gamma:

  EBIC hyperparameter for `"glasso"`, `"huge"`, `"ggm"`, `"ising"` and
  `"mgm"`. `NULL` (default) keeps each method's own default.

- labels:

  Optional node labels, one per node.

- vars:

  The variables to use as nodes (default: all). A name range
  (`motivation:regulation`), an index range (`3:9`), names
  (`c(joy, fear)` or `c("joy", "fear")`) or a single name, as in the
  `select` argument of [`subset()`](https://rdrr.io/r/base/subset.html).
  For an event log the actions are the nodes.

- group:

  Optional grouping column(s). One network is estimated per level and a
  [psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  is returned.

- within:

  Optional person-identifier column. Each variable is centred and scaled
  within persons over all rows before estimation (each person's mean
  subtracted and the result divided by the person's standard deviation,
  as `mlVAR::mlVAR(scaleWithin = TRUE)`), and the person column is not a
  node. With `group`, rows are matched by person so that
  [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  can pair two levels (`paired = TRUE`). Without `group`, a person's
  rows are pooled and are not independent. Numeric tables only.

- source:

  Kind of input: `"auto"` (default; an event log when `actor` or `id` is
  given, otherwise a numeric table), `"data"` or `"eventdata"`.

- actor, action, session, time:

  Event-log columns: the actor, the action (default `"Action"`), an
  optional session within actor, and an optional time used only to split
  actors into sessions. Giving `actor` selects event data.

- compute_sessions, time_threshold:

  Whether to split each actor into sessions where the gap in `time`
  exceeds `time_threshold` (defaults `TRUE` and 900, i.e. 15 minutes for
  times in seconds). Has no effect without `time`.

- id:

  Actor column of an event log already reduced to a numeric table with
  one row per occasion.

- standardize:

  For event data with several occasions per actor: `TRUE` (default)
  centres each actor's occasions and fits one network; `FALSE` returns
  separate within-actor and between-actor networks.

- ...:

  Further arguments for the estimator, such as `cor_method` for the
  correlation-based methods, `npn` for `"huge"`, `rule` for the Ising
  methods and `alpha` for `"cor"`, `"pcor"` and `"ising_sampler"`. For
  an event table, `na_method = "listwise"` removes incomplete rows
  before centring; `"pairwise"` (default) uses the available values.

## Value

A fitted network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md));
with `group`, a
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
with one network per level; for event data with several occasions per
actor and `standardize = FALSE`, a
[psychnet_multilevel](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
pair of within-actor and between-actor networks. Edges are read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## Details

`method` accepts the package's short names and the names used by
`qgraph` and `bootnet` (`"EBICglasso"`, `"ggmModSelect"`, `"TMFG"`,
`"LoGo"`, `"IsingFit"`, `"IsingSampler"`), ignoring case and separators;
all resolve to the same estimators, and the fitted network records the
short name.
[`net_crosswalk()`](https://pak.dynasite.org/psychnets/reference/net_crosswalk.md)
maps the arguments of the reference packages to those of this package.

## Conditions

`psychnet_bad_input` when `method` has several values or is unknown,
`vars` names a missing column, `within` does not name one column of a
numeric table, or event data are combined with a method other than the
Gaussian ones; `psychnet_ignored_argument` (warning) when `threshold` or
`gamma` is set for a method that does not use it. With `group`, the
conditions listed in
[psychnet_group](https://pak.dynasite.org/psychnets/reference/psychnet_group.md).
The chosen estimator raises its own conditions.

## See also

[psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
for the fitted object;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
to check a regularized fit;
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
and
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
for analyses of fitted networks;
[`event_frequencies()`](https://pak.dynasite.org/psychnets/reference/event_frequencies.md)
for the event-log conversion.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 5), 200, 5)
colnames(x) <- c("joy", "fear", "calm", "anger", "trust")
psychnet(x, method = "glasso")
#> <psychnet> glasso network
#>   nodes: 5   edges: 0   (undirected)
#>   lambda: 0.06852   gamma: 0.5
#>   optimality (KKT residual): 1.39e-17
psychnet(x, method = "pcor", vars = joy:anger)     # name range
#> <psychnet> pcor network
#>   nodes: 4   edges: 6   (undirected)
psychnet(x, method = "glasso", vars = 1:3)         # index range
#> <psychnet> glasso network
#>   nodes: 3   edges: 0   (undirected)
#>   lambda: 0.06817   gamma: 0.5
#>   optimality (KKT residual): 0.00e+00
psychnet(x, method = "EBICglasso")                 # qgraph alias, same result
#> <psychnet> glasso network
#>   nodes: 5   edges: 0   (undirected)
#>   lambda: 0.06852   gamma: 0.5
#>   optimality (KKT residual): 1.39e-17

ev <- data.frame(
  Actor  = rep(paste0("s", 1:30), each = 20),
  Action = sample(c("read", "quiz", "note", "watch"), 600, replace = TRUE))
psychnet(ev, actor = "Actor", action = "Action")   # event data, auto-detected
#> <psychnet> glasso network
#>   nodes: 4   edges: 6   (undirected)
#>   lambda: 0.004501   gamma: 0.5
#>   optimality (KKT residual): 1.07e-09

d <- data.frame(g = rep(c("A", "B"), each = 100),
                matrix(stats::rnorm(200 * 4), 200, 4,
                       dimnames = list(NULL, paste0("V", 1:4))))
psychnet(d, method = "glasso", group = "g")         # one network per level
#> <psychnet_group> 2 networks by g (method: glasso)
#>  group nodes edges   n
#>      A     4     0 100
#>      B     4     0 100

# Within-person networks, one per week, scaled within each student over all
# four weeks; compare two weeks (paired by student automatically).
nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                 within = "student")
nets
#> <psychnet_group> 4 networks by week (method: glasso)
#>  group nodes edges   n
#>      1     5     6 100
#>      2     5     6 100
#>      3     5     7 100
#>      4     5     7 100
set.seed(1)
net_compare(nets, c(1, 4), iter = 3)  # iter >= 1000 in practice
#> Network Comparison Test: 1 vs 4
#> 3 permutations, paired groups, n = 100 and 100
#> 
#>   Test                  Statistic       p
#>   M  network structure      0.084   1.000 
#>   S  global strength        0.006   1.000 
#>        1 1.516 vs 4 1.510
#>   E  edges differing      0 of 10
#> 
#>   Edges, largest difference first
#>   Edge                       1       4  |diff|       p
#>   planning - effort      0.115   0.199   0.084   0.250 
#>   effort - enjoyment     0.382   0.314   0.068   0.750 
#>   efficacy - anxiety     0.000   0.064   0.064   0.250 
#>   efficacy - planning    0.503   0.440   0.063   1.000 
#>   efficacy - effort      0.312   0.358   0.046   1.000 
#>   planning - enjoyment   0.165   0.125   0.040   0.500 
#>   efficacy - enjoyment   0.039   0.011   0.029   1.000 
#>   effort - anxiety       0.000   0.000   0.000   1.000 
#>   enjoyment - anxiety    0.000   0.000   0.000   1.000 
#>   planning - anxiety     0.000   0.000   0.000   1.000 
#> 
#>   * p < 0.05
```
