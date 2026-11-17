# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...),
# extended to read raw event logs (actor / action / time) directly.

#' Estimate a psychometric network with any of the package's estimators
#'
#' @description
#' `psychnet()` estimates a psychometric network with any of the
#' cross-sectional estimators, selected by `method`. A psychometric network
#' represents variables, such as
#' questionnaire items, symptoms or behaviours, as nodes and their conditional
#' or marginal associations as edges.
#'
#' The estimators differ in what an edge means: partial-correlation networks
#' describe the association between two variables after controlling for all
#' others, regularized networks retain only the associations the data support,
#' and Ising and mixed graphical models apply the same idea to binary and
#' mixed data. Besides `method`, the arguments select the variables (`vars`),
#' estimate one network per group (`group`), scale each variable within
#' persons first (`within`), or treat `data` as an event log (`actor`,
#' `action`, `session`, `time`), which is converted to action frequencies per
#' occasion before a network is estimated.
#'
#' `psychnet()` returns a fitted network; with `group`, one network per
#' group, and for some event data a pair of within-actor and between-actor
#' networks.
#'
#' Each method can also be called directly, which exposes its full set of
#' options. Panel data measured in a few waves are modelled by
#' [panel_gvar()], which `psychnet()` does not route. The method-specific
#' pages describe each model: [ebic_glasso()], [ggm_modselect()], [huge_network()], [tmfg_network()],
#' [logo_network()], [relimp_network()], [cor_network()], [pcor_network()],
#' [ising_fit()], [ising_sampler()] and [mgm_fit()].
#'
#' @details
#' `method` accepts the package's short names and the names used by `qgraph`
#' and `bootnet` (`"EBICglasso"`, `"ggmModSelect"`, `"TMFG"`, `"LoGo"`,
#' `"IsingFit"`, `"IsingSampler"`), ignoring case and separators; all resolve
#' to the same estimators, and the fitted network records the short name.
#' [net_crosswalk()] maps the arguments of the reference packages to those of
#' this package.
#'
#' @param data A data frame or matrix with one row per observation and one
#'   column per variable, or a long event log (one row per event).
#' @param method The estimator, a single name: `"glasso"` (default), `"ggm"`,
#'   `"huge"`, `"tmfg"`, `"logo"`, `"relimp"`, `"cor"`, `"pcor"`, `"ising"`,
#'   `"ising_sampler"` or `"mgm"`, or a `qgraph`/`bootnet` alias. Event data
#'   accept `"glasso"`, `"pcor"`, `"ggm"`, `"huge"`, `"tmfg"` and `"logo"`.
#' @param vars The variables to use as nodes (default: all). A name range
#'   (`motivation:regulation`), an index range (`3:9`), names (`c(joy, fear)`
#'   or `c("joy", "fear")`) or a single name, as in the `select` argument of
#'   [subset()]. For an event log the actions are the nodes.
#' @param group Optional grouping column(s). One network is estimated per
#'   level and a [psychnet_group] is returned.
#' @param within Optional person-identifier column. Each variable is centred
#'   and scaled within persons over all rows before estimation (each
#'   person's mean subtracted and the result divided by the person's standard
#'   deviation, as `mlVAR::mlVAR(scaleWithin = TRUE)`), and the person column
#'   is not a node. With `group`,
#'   rows are matched by person so that [net_compare()] can pair two levels
#'   (`paired = TRUE`). Without `group`, a person's rows are pooled and are
#'   not independent. Numeric tables only.
#' @param source Kind of input: `"auto"` (default; an event log when `actor`
#'   or `id` is given, otherwise a numeric table), `"data"` or `"eventdata"`.
#' @param actor,action,session,time Event-log columns: the actor, the action
#'   (default `"Action"`), an optional session within actor, and an optional
#'   time used only to split actors into sessions. Giving `actor` selects
#'   event data.
#' @param compute_sessions,time_threshold Whether to split each actor into
#'   sessions where the gap in `time` exceeds `time_threshold` (defaults
#'   `TRUE` and 900, i.e. 15 minutes for times in seconds). Has no effect
#'   without `time`.
#' @param id Actor column of an event log already reduced to a numeric table
#'   with one row per occasion.
#' @param standardize For event data with several occasions per actor: `TRUE`
#'   (default) centres each actor's occasions and fits one network; `FALSE`
#'   returns separate within-actor and between-actor networks.
#' @param threshold Absolute edge weight below which edges are set to zero, for
#'   `"cor"`, `"pcor"`, `"glasso"`, `"huge"`, `"ggm"` and `"logo"`: a single
#'   number >= 0 (default `0`). For `"mgm"`, give its rule (`"LW"`, `"HW"` or
#'   `"none"`); the default `0` leaves the [mgm_fit()] default in place.
#' @param gamma EBIC hyperparameter for `"glasso"`, `"huge"`, `"ggm"`,
#'   `"ising"` and `"mgm"`. `NULL` (default) keeps each method's own default.
#' @param labels Optional node labels, one per node.
#' @param ... Further arguments for the estimator, such as `cor_method` for
#'   the correlation-based methods, `npn` for `"huge"`, `rule` for the Ising
#'   methods and `alpha` for `"cor"`, `"pcor"` and `"ising_sampler"`. For an
#'   event table, `na_method = "listwise"` removes incomplete rows before
#'   centring; `"pairwise"` (default) uses the available values.
#' @return A fitted network of class `psychnet` ([psychnet-class]); with
#'   `group`, a [psychnet_group] with one network per level; for event data
#'   with several occasions per actor and `standardize = FALSE`, a
#'   [psychnet_multilevel] pair of within-actor and between-actor networks.
#'   Edges are read with `as.data.frame()`.
#' @section Conditions:
#'   `psychnet_bad_input` when `method` has several values or is unknown,
#'   `vars` names a missing column, `within` does not name one column of a
#'   numeric table, or event data are combined with a method other than the
#'   Gaussian ones; `psychnet_ignored_argument` (warning) when `threshold` or
#'   `gamma` is set for a method that does not use it. With `group`, the
#'   conditions listed in [psychnet_group]. The chosen estimator raises its
#'   own conditions.
#' @seealso [psychnet-class] for the fitted object; [certificate()] to check
#'   a regularized fit; [net_centralities()], [net_boot()] and
#'   [net_compare()] for analyses of fitted networks; [event_frequencies()]
#'   for the event-log conversion.
#' @examples
#' set.seed(1)
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
#'
#' # Within-person networks, one per week, scaled within each student over all
#' # four weeks; compare two weeks (paired by student automatically).
#' nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
#'                  within = "student")
#' nets
#' set.seed(1)
#' net_compare(nets, c(1, 4), iter = 3)  # iter >= 1000 in practice
#' @export
psychnet <- function(data,
                     method = c("glasso", "cor", "pcor", "ising", "mgm",
                                "huge", "ggm", "tmfg", "logo", "relimp",
                                "ising_sampler"),
                     threshold = 0, gamma = NULL, labels = NULL,
                     vars = NULL, group = NULL, within = NULL,
                     source = c("auto", "data", "eventdata"),
                     actor = NULL, action = "Action", session = NULL,
                     time = NULL, compute_sessions = TRUE, time_threshold = 900,
                     id = NULL, standardize = TRUE, ...) {
  vars_q <- substitute(vars)                           # capture for tidy select
  vars_env <- parent.frame()
  source <- .psn_match_arg(source)
  if (!missing(method) && length(method) != 1L) {
    .psn_stop(sprintf(paste0(
      "`method` must be a single estimator name, not %d values (%s); fit each ",
      "method with its own call."), length(method),
      paste(method, collapse = ", ")))
  }
  method <- .resolve_method(method)

  if (source == "auto")
    source <- if (!is.null(actor) || !is.null(id)) "eventdata" else "data"

  # Within-person scaling over all rows, before any grouping, so the person
  # means and SDs use every occasion.
  if (!is.null(within)) {
    if (!identical(source, "data"))
      .psn_stop("`within` is for a numeric table with one row per occasion; ",
                "for an event log use `actor` (and `standardize`).")
    if (!is.data.frame(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!(is.character(within) && length(within) == 1L && within %in% names(data)))
      .psn_stop("`within` must name one column of `data` (the person identifier).")
    nodes <- if (is.null(vars_q)) {
      setdiff(names(data)[vapply(data, is.numeric, logical(1))], c(within, group))
    } else {
      setdiff(.select_names(data, vars_q, vars_env), c(within, group))
    }
    data <- .scale_within(data, id = within, vars = nodes)
  }

  # Group mode: one network per level of `group`, returned as a psychnet_group.
  # An ignored threshold / gamma is reported once here; the per-level refits
  # would otherwise repeat the same warning for every group.
  if (!is.null(group)) {
    if (identical(source, "data")) .warn_ignored_args(method, threshold, gamma)
    return(withCallingHandlers(
      .psychnet_grouped(data = data, group = group, source = source,
                        vars_q = vars_q, vars_env = vars_env,
                        method = method, actor = actor, action = action,
                        session = session, time = time,
                        compute_sessions = compute_sessions,
                        time_threshold = time_threshold, id = id,
                        standardize = standardize, threshold = threshold,
                        gamma = gamma, labels = labels, within = within,
                        dots = list(...)),
      psychnet_ignored_argument = function(w) invokeRestart("muffleWarning")))
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
  if (!is.null(within)) data[[within]] <- NULL           # the person is not a node
  if (!is.null(vars_q)) data <- .select_cols(data, vars_q, vars_env)
  .warn_ignored_args(method, threshold, gamma)

  fn <- switch(
    method,
    cor = cor_network, pcor = pcor_network, glasso = ebic_glasso,
    ising = ising_fit, mgm = mgm_fit, huge = huge_network,
    ggm = ggm_modselect, tmfg = tmfg_network, logo = logo_network,
    relimp = relimp_network, ising_sampler = ising_sampler,
    .psn_stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
  # Forward each shared argument only to the methods that accept it, so the
  # callee's own defaults stand otherwise (notably gamma: ising/mgm = 0.25).
  args <- c(list(data = data, labels = labels), list(...))
  if (method %in% .threshold_methods)
    args$threshold <- threshold
  if (identical(method, "mgm") && is.character(threshold))
    args$threshold <- threshold
  if (!is.null(gamma) && method %in% .gamma_methods)
    args$gamma <- gamma
  do.call(fn, args)
}

# The estimators that take the shared front-door arguments. psychnet()
# forwards `threshold` / `gamma` only to these (so every other method keeps its
# own defaults); the sets are shared with .warn_ignored_args() so the routing
# and the warning cannot drift apart.
#' @noRd
.threshold_methods <- c("cor", "pcor", "glasso", "huge", "ggm", "logo")
#' @noRd
.gamma_methods <- c("glasso", "huge", "ggm", "ising", "mgm")

# Warn (class `psychnet_ignored_argument`) when the caller set `threshold` or
# `gamma` for a method that does not use it. Only a non-default value counts:
# `threshold = 0` and `gamma = NULL` are the defaults that internal callers
# (group mode, the resampling verbs) forward unchanged, and they change nothing.
#' @noRd
.warn_ignored_args <- function(method, threshold, gamma) {
  ignored <- character()
  threshold_default <- is.numeric(threshold) &&
    identical(as.numeric(threshold), 0)
  takes_threshold <- method %in% .threshold_methods ||
    (identical(method, "mgm") && is.character(threshold))
  if (!threshold_default && !takes_threshold) ignored <- c(ignored, "threshold")
  if (!is.null(gamma) && !method %in% .gamma_methods)
    ignored <- c(ignored, "gamma")
  if (length(ignored)) {
    .psn_warn(sprintf(paste0(
      "%s ignored: method = \"%s\" does not use %s.%s"),
      paste0("`", ignored, "`", collapse = " and "), method,
      if (length(ignored) > 1L) "them" else "it",
      if ("threshold" %in% ignored && identical(method, "mgm"))
        " For mgm, give its character rule (\"LW\", \"HW\" or \"none\")."
      else ""),
      class = "psychnet_ignored_argument")
  }
  invisible(ignored)
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
    .psn_stop("`vars` selection requires named columns.", call. = FALSE)
  nl <- as.list(seq_along(cn))
  names(nl) <- cn
  idx <- eval(vars_q, nl, env)
  if (is.character(idx)) {
    miss <- setdiff(idx, cn)
    if (length(miss))
      .psn_stop("Columns not found in data: ", paste(miss, collapse = ", "),
           call. = FALSE)
    return(idx)
  }
  cn[idx]
}

# Tidy column selection returning the subset data (preserving frame/matrix).
#' @noRd
.select_cols <- function(data, vars_q, env) {
  if (is.null(colnames(data)))
    .psn_stop("`vars` selection requires named columns.", call. = FALSE)
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
  if (is.na(out)) .psn_stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  unname(out)
}
