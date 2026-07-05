# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...),
# extended to read raw event logs (actor / action / time) directly.

#' Estimate a psychometric network
#'
#' The package's main entry point. Supply `data` and a `method`; everything else
#' is optional fine-grained control. `psychnet()` reads two kinds of input and
#' picks the right one automatically (`source = "auto"`):
#'
#' * **a numeric table** (`source = "data"`) -- the ordinary cross-sectional
#'   case: one row per observation, one column per variable;
#' * **a long event log** (`source = "eventdata"`) -- one row per event, with an
#'   `actor` and an `action` column (and optionally `session` / `time`). The log
#'   is converted to action frequencies with [event_frequencies()] and, when
#'   actors contribute several occasions, decomposed into within- and
#'   between-actor networks. Passing `actor` switches this on automatically.
#'
#' `method` speaks the package's own short vocabulary -- `"glasso"`, `"ggm"`,
#' `"tmfg"`, `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`, `"huge"`,
#' `"mgm"`, `"cor"`, `"pcor"`. For interoperability it **also** accepts the
#' `qgraph`/`bootnet` spellings (`"EBICglasso"`, `"ggmModSelect"`, `"TMFG"`,
#' `"LoGo"`, `"IsingFit"`, `"IsingSampler"`), which resolve to the same
#' estimators. Whichever you pass in, the stored `$method` is the short name.
#'
#' @param data A numeric data frame / matrix (rows = observations), or a long
#'   event log when `actor` is supplied.
#' @param method Estimator. One of `"glasso"` (default), `"ggm"`, `"tmfg"`,
#'   `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`, `"huge"`, `"mgm"`,
#'   `"cor"`, `"pcor"`. The `qgraph`/`bootnet` names are accepted as aliases
#'   (e.g. `"EBICglasso"` -> `"glasso"`, `"ggmModSelect"` -> `"ggm"`).
#'   Event data is restricted to the Gaussian graphical methods.
#' @param vars Which variables to build the network on. Defaults to every
#'   variable. Selected the tidy way: a name range (`motivation:regulation`),
#'   a column-index range (`3:9`), a vector of names (`c(joy, fear)` or
#'   `c("joy", "fear")`), or a single name -- the same grammar as `subset()`'s
#'   `select=`. (For an event log, actions become the variables, so `vars`
#'   selects feature columns only when `data` is already a numeric table.)
#' @param group Optional grouping column(s). When supplied, one network is
#'   estimated per level of `group` and a `psychnet_group` (a named list of
#'   networks) is returned, which plots as a grid and is iterated per level by
#'   the framework verbs.
#' @param source Input kind: `"auto"` (default; an event log when `actor` is
#'   given, otherwise a numeric table), `"data"`, or `"eventdata"`.
#' @param actor,action,session,time Event-log columns. `actor` (and, by default,
#'   `action = "Action"`) name the subject and the event; `session` is an
#'   explicit within-actor grouping, and `time` is used only to compute sessions
#'   from gaps. Supplying `actor` selects `source = "eventdata"`.
#' @param compute_sessions,time_threshold Split each actor into sessions from
#'   `time` gaps (a new session starts when the gap exceeds `time_threshold`,
#'   default 900 s). `compute_sessions` defaults to `TRUE`; it is a no-op when no
#'   `time` is supplied.
#' @param id Actor column when an event log has already been reduced to a numeric
#'   feature table (one row per occasion) rather than raw events.
#' @param standardize For nested event data (several occasions per actor),
#'   `TRUE` (default) removes the actor clustering by person-centering and fits a
#'   single network; `FALSE` returns the within/between pair instead.
#' @param threshold Absolute-weight threshold below which edges are zeroed
#'   (forwarded only to the methods that take it: `cor`, `pcor`, `glasso`,
#'   `huge`, `ggm`, `logo`).
#' @param gamma EBIC hyperparameter. `NULL` (default) keeps each method's own
#'   default (0.5 for the regularized Gaussian graphical models, 0 for `ggm`,
#'   0.25 for `ising`/`mgm`); set it to override. Forwarded only to the
#'   regularized methods.
#' @param labels Optional node labels.
#' @param ... Passed to the underlying estimator (e.g. `cor_method=` for the
#'   correlation-based methods, `npn=` for `"huge"`, `rule=` for the Ising
#'   methods, `alpha=` for the correlation / `"ising_sampler"` methods).
#' @return A `psychnet` object; a `psychnet_group` (named list of networks) when
#'   `group` is supplied; or, for nested event data with `standardize = FALSE`, a
#'   `psychnet_multilevel` object carrying `$within` and `$between` networks.
#' @examples
#' x <- matrix(stats::rnorm(200 * 5), 200, 5)
#' colnames(x) <- c("joy", "fear", "calm", "anger", "trust")
#' psychnet(x, method = "glasso")
#' psychnet(x, method = "pcor", vars = joy:anger)     # name range
#' psychnet(x, method = "glasso", vars = 1:3)         # index range
#' psychnet(x, method = "EBICglasso")                 # qgraph alias, same result
#'
#' ev <- data.frame(
#'   Actor  = rep(paste0("s", 1:30), each = 20),
#'   Action = sample(c("read", "quiz", "note", "watch"), 600, replace = TRUE))
#' psychnet(ev, actor = "Actor", action = "Action")   # event data, auto-detected
#'
#' d <- data.frame(g = rep(c("A", "B"), each = 100),
#'                 matrix(stats::rnorm(200 * 4), 200, 4,
#'                        dimnames = list(NULL, paste0("V", 1:4))))
#' psychnet(d, method = "glasso", group = "g")         # one network per level
#' @export
psychnet <- function(data,
                     method = c("glasso", "cor", "pcor", "ising", "mgm",
                                "huge", "ggm", "tmfg", "logo", "relimp",
                                "ising_sampler"),
                     threshold = 0, gamma = NULL, labels = NULL,
                     vars = NULL, group = NULL,
                     source = c("auto", "data", "eventdata"),
                     actor = NULL, action = "Action", session = NULL,
                     time = NULL, compute_sessions = TRUE, time_threshold = 900,
                     id = NULL, standardize = TRUE, ...) {
  vars_q <- substitute(vars)                           # capture for tidy select
  vars_env <- parent.frame()
  source <- match.arg(source)
  method <- .resolve_method(method)

  if (source == "auto")
    source <- if (!is.null(actor) || !is.null(id)) "eventdata" else "data"

  # Group mode: one network per level of `group`, returned as a psychnet_group.
  if (!is.null(group)) {
    return(.psychnet_grouped(data = data, group = group, source = source,
                             vars_q = vars_q, vars_env = vars_env,
                             method = method, actor = actor, action = action,
                             session = session, time = time,
                             compute_sessions = compute_sessions,
                             time_threshold = time_threshold, id = id,
                             standardize = standardize, threshold = threshold,
                             gamma = gamma, labels = labels, dots = list(...)))
  }

  if (source == "eventdata") {
    vars_sel <- if (is.null(vars_q)) NULL else .select_names(data, vars_q, vars_env)
    return(.net_eventdata(data, actor = actor, action = action,
                          session = session, time = time,
                          compute_sessions = compute_sessions,
                          time_threshold = time_threshold, id = id,
                          vars = vars_sel, standardize = standardize,
                          method = method, threshold = threshold,
                          gamma = gamma, labels = labels, ...))
  }

  # Cross-sectional numeric table: optional tidy column selection, then route.
  if (!is.null(vars_q)) data <- .select_cols(data, vars_q, vars_env)

  fn <- switch(
    method,
    cor = cor_network, pcor = pcor_network, glasso = ebic_glasso,
    ising = ising_fit, mgm = mgm_fit, huge = huge_network,
    ggm = ggm_modselect, tmfg = tmfg_network, logo = logo_network,
    relimp = relimp_network, ising_sampler = ising_sampler,
    stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
  # Forward each shared argument only to the methods that accept it, so the
  # callee's own defaults stand otherwise (notably gamma: ising/mgm = 0.25).
  args <- c(list(data = data, labels = labels), list(...))
  if (method %in% c("cor", "pcor", "glasso", "huge", "ggm", "logo"))
    args$threshold <- threshold
  if (!is.null(gamma) &&
      method %in% c("glasso", "huge", "ggm", "ising", "mgm"))
    args$gamma <- gamma
  do.call(fn, args)
}

# Resolve a `vars` expression to the selected column names, using the same
# grammar as subset()'s `select=`: column names are bound to their positions, so
# a name range (a:c), an index range (3:9), c(name, name) / c("a","b") and a bare
# name all evaluate to a set of column indices. NULL means "all columns".
#' @noRd
.select_names <- function(data, vars_q, env) {
  if (is.null(vars_q)) return(NULL)
  cn <- if (is.data.frame(data)) names(data) else colnames(data)
  if (is.null(cn))
    stop("`vars` selection requires named columns.", call. = FALSE)
  nl <- as.list(seq_along(cn))
  names(nl) <- cn
  idx <- eval(vars_q, nl, env)
  if (is.character(idx)) {
    miss <- setdiff(idx, cn)
    if (length(miss))
      stop("Columns not found in data: ", paste(miss, collapse = ", "),
           call. = FALSE)
    return(idx)
  }
  cn[idx]
}

# Tidy column selection returning the subset data (preserving frame/matrix).
#' @noRd
.select_cols <- function(data, vars_q, env) {
  if (is.null(colnames(data)))
    stop("`vars` selection requires named columns.", call. = FALSE)
  sel <- .select_names(data, vars_q, env)
  data[, sel, drop = FALSE]
}

# Resolve a method name to its canonical short form. Accepts the package's own
# short names and the qgraph/bootnet spellings (dual vocabulary); matching is
# case- and separator-insensitive, so "EBICglasso", "ebic_glasso" and "glasso"
# all land on "glasso".
#' @noRd
.resolve_method <- function(method) {
  if (length(method) > 1L) method <- method[1L]
  key <- gsub("[^a-z]", "", tolower(method))         # strip case/underscores
  canon <- c(
    cor = "cor", correlation = "cor",
    pcor = "pcor", partial = "pcor", partialcor = "pcor",
    glasso = "glasso", ebicglasso = "glasso", ebic = "glasso",
    huge = "huge", npn = "huge", nonparanormal = "huge",
    ggm = "ggm", ggms = "ggm", ggmmodselect = "ggm", modselect = "ggm",
    stepwise = "ggm",
    tmfg = "tmfg",
    logo = "logo",
    relimp = "relimp", relativeimportance = "relimp", lmg = "relimp",
    ising = "ising", isingfit = "ising",
    isingsampler = "ising_sampler",
    mgm = "mgm"
  )
  out <- canon[key]
  if (is.na(out)) stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  unname(out)
}
