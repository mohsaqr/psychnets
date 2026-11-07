# Estimate a psychometric network

The package's main entry point. Supply `data` and a `method`; everything
else is optional fine-grained control. `psychnet()` reads two kinds of
input and picks the right one automatically (`source = "auto"`):

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

  A numeric data frame / matrix (rows = observations), or a long event
  log when `actor` is supplied.

- method:

  Estimator, a single name. One of `"glasso"` (default), `"ggm"`,
  `"tmfg"`, `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`,
  `"huge"`, `"mgm"`, `"cor"`, `"pcor"`; several names at once are an
  error (`psychnet_bad_input`). The `qgraph`/`bootnet` names are
  accepted as aliases (e.g. `"EBICglasso"` -\> `"glasso"`,
  `"ggmModSelect"` -\> `"ggm"`). Event data is restricted to the
  Gaussian graphical methods.

- threshold:

  Absolute-weight threshold below which edges are zeroed for `cor`,
  `pcor`, `glasso`, `huge`, `ggm`, and `logo`: a single number \>= 0
  (validated by the estimator). For `mgm`, supply its character rule
  (`"LW"`, `"HW"`, or `"none"`); the numeric default leaves
  [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)'s
  own default unchanged. A non-zero value for any other method is
  ignored with a `psychnet_ignored_argument` warning.

- gamma:

  EBIC hyperparameter. `NULL` (default) keeps each method's own default
  (0.5 for the regularized Gaussian graphical models, 0 for `ggm`, 0.25
  for `ising`/`mgm`); set it to override. Forwarded only to the
  regularized methods; a value given for another method is ignored with
  a `psychnet_ignored_argument` warning.

- labels:

  Optional node labels.

- vars:

  Which variables to build the network on. Defaults to every variable.
  Selected the tidy way: a name range (`motivation:regulation`), a
  column-index range (`3:9`), a vector of names (`c(joy, fear)` or
  `c("joy", "fear")`), or a single name – the same grammar as
  [`subset()`](https://rdrr.io/r/base/subset.html)'s `select=`. (For an
  event log, actions become the variables, so `vars` selects feature
  columns only when `data` is already a numeric table.)

- group:

  Optional grouping column(s). When supplied, one network is estimated
  per level of `group` and a `psychnet_group` (a named list of networks)
  is returned, which plots as a grid and is iterated per level by the
  framework verbs.

- within:

  Optional name of a person-identifier column. When supplied, every
  network variable is centered and scaled *within persons* first (see
  [`scale_within()`](https://pak.dynasite.org/psychnets/reference/scale_within.md);
  mlVAR's `scaleWithin`), over all rows – so with `group` (e.g.
  occasions) the person means and SDs use every occasion – and the
  person column is not a node. With `group`, the result remembers which
  person each row is, so `net_compare(x, levels, paired = TRUE)` pairs
  the two levels by person. Without `group`, rows of the same person are
  pooled and are not independent. Cross-sectional data only.

- source:

  Input kind: `"auto"` (default; an event log when `actor` is given,
  otherwise a numeric table), `"data"`, or `"eventdata"`.

- actor, action, session, time:

  Event-log columns. `actor` (and, by default, `action = "Action"`) name
  the subject and the event; `session` is an explicit within-actor
  grouping, and `time` is used only to compute sessions from gaps.
  Supplying `actor` selects `source = "eventdata"`.

- compute_sessions, time_threshold:

  Split each actor into sessions from `time` gaps (a new session starts
  when the gap exceeds `time_threshold`, default 900 s).
  `compute_sessions` defaults to `TRUE`; it is a no-op when no `time` is
  supplied.

- id:

  Actor column when an event log has already been reduced to a numeric
  feature table (one row per occasion) rather than raw events.

- standardize:

  For nested event data (several occasions per actor), `TRUE` (default)
  removes the actor clustering by person-centering and fits a single
  network; `FALSE` returns the within/between pair instead.

- ...:

  Passed to the underlying estimator (e.g. `cor_method=` for the
  correlation-based methods, `npn=` for `"huge"`, `rule=` for the Ising
  methods, `alpha=` for the correlation / `"ising_sampler"` methods).
  For event feature tables, `na_method = "listwise"` removes incomplete
  rows before person-centering; `"pairwise"` (default) uses available
  observations for each person's variable means and each correlation.

## Value

A `psychnet` object; a `psychnet_group` (named list of networks) when
`group` is supplied; or, for nested event data with
`standardize = FALSE`, a `psychnet_multilevel` object carrying `$within`
and `$between` networks.

## Details

- **a numeric table** (`source = "data"`) – the ordinary cross-sectional
  case: one row per observation, one column per variable;

- **a long event log** (`source = "eventdata"`) – one row per event,
  with an `actor` and an `action` column (and optionally `session` /
  `time`). The log is converted to action frequencies with
  [`event_frequencies()`](https://pak.dynasite.org/psychnets/reference/event_frequencies.md)
  and, when actors contribute several occasions, decomposed into within-
  and between-actor networks. Passing `actor` switches this on
  automatically.

`method` speaks the package's own short vocabulary – `"glasso"`,
`"ggm"`, `"tmfg"`, `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`,
`"huge"`, `"mgm"`, `"cor"`, `"pcor"`. For interoperability it **also**
accepts the `qgraph`/`bootnet` spellings (`"EBICglasso"`,
`"ggmModSelect"`, `"TMFG"`, `"LoGo"`, `"IsingFit"`, `"IsingSampler"`),
which resolve to the same estimators. Whichever you pass in, the stored
`$method` is the short name.

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
