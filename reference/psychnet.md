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

  Estimator. One of `"glasso"` (default), `"ggm"`, `"tmfg"`, `"logo"`,
  `"relimp"`, `"ising"`, `"ising_sampler"`, `"huge"`, `"mgm"`, `"cor"`,
  `"pcor"`. The `qgraph`/`bootnet` names are accepted as aliases (e.g.
  `"EBICglasso"` -\> `"glasso"`, `"ggmModSelect"` -\> `"ggm"`). Event
  data is restricted to the Gaussian graphical methods.

- threshold:

  Absolute-weight threshold below which edges are zeroed (forwarded only
  to the methods that take it: `cor`, `pcor`, `glasso`, `huge`, `ggm`,
  `logo`).

- gamma:

  EBIC hyperparameter. `NULL` (default) keeps each method's own default
  (0.5 for the regularized Gaussian graphical models, 0 for `ggm`, 0.25
  for `ising`/`mgm`); set it to override. Forwarded only to the
  regularized methods.

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
x <- matrix(stats::rnorm(200 * 5), 200, 5)
colnames(x) <- c("joy", "fear", "calm", "anger", "trust")
psychnet(x, method = "glasso")
#> <psychnet> glasso network
#>   nodes: 5   edges: 0   (undirected)
#>   lambda: 0.1688   gamma: 0.5
#>   optimality (KKT residual): 0.00e+00
psychnet(x, method = "pcor", vars = joy:anger)     # name range
#> <psychnet> pcor network
#>   nodes: 4   edges: 6   (undirected)
psychnet(x, method = "glasso", vars = 1:3)         # index range
#> <psychnet> glasso network
#>   nodes: 3   edges: 0   (undirected)
#>   lambda: 0.06016   gamma: 0.5
#>   optimality (KKT residual): 0.00e+00
psychnet(x, method = "EBICglasso")                 # qgraph alias, same result
#> <psychnet> glasso network
#>   nodes: 5   edges: 0   (undirected)
#>   lambda: 0.1688   gamma: 0.5
#>   optimality (KKT residual): 0.00e+00

ev <- data.frame(
  Actor  = rep(paste0("s", 1:30), each = 20),
  Action = sample(c("read", "quiz", "note", "watch"), 600, replace = TRUE))
psychnet(ev, actor = "Actor", action = "Action")   # event data, auto-detected
#> <psychnet> glasso network
#>   nodes: 4   edges: 6   (undirected)
#>   lambda: 0.005028   gamma: 0.5
#>   optimality (KKT residual): 8.00e-10

d <- data.frame(g = rep(c("A", "B"), each = 100),
                matrix(stats::rnorm(200 * 4), 200, 4,
                       dimnames = list(NULL, paste0("V", 1:4))))
psychnet(d, method = "glasso", group = "g")         # one network per level
#> <psychnet_group> 2 networks by g (method: glasso)
#>  group nodes edges   n
#>      A     4     0 100
#>      B     4     0 100
```
