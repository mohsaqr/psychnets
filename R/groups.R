# Group mode: estimate one network per level of a grouping variable, mirroring
# Nestimate's `netobject_group`. psychnet(data, group = ) returns a
# `psychnet_group`: a named list of single networks keyed by group level. The
# class is `c("psychnet_group", "netobject_group")` so cograph::splot() lays the
# group out in a grid (it already dispatches on "netobject_group"), and the
# framework verbs (centrality, predict, bootstrap, stability, compare) detect the
# group object and iterate per level automatically -- the same transparent
# `inherits(x, ...)` dispatch Nestimate uses, returning a `*_group` result.
#
# The lean single networks never store their data; the *container* keeps one copy
# of the per-level fit data (`.subsets`) plus the estimation call (`.call`) so the
# data-first verbs (net_boot / net_stability / net_compare) can re-estimate each
# level. That is the re-estimable unit, leaving the per-net objects lean.

#' @noRd
.new_psychnet_group <- function(nets, group_col, source, subsets, call) {
  structure(nets,
            group_col = group_col,
            source    = source,
            subsets   = subsets,
            call      = call,
            class     = c("psychnet_group", "netobject_group"))
}

# Combine tuple identities before creating display labels. Pasting raw values
# merges distinct tuples whenever a value itself contains the separator.
#' @noRd
.tuple_factor <- function(data, cols, sep = ".", lex.order = TRUE) {
  if (length(cols) == 1L) return(droplevels(factor(data[[cols]])))
  codes <- lapply(data[cols], function(x) as.integer(factor(x)))
  key <- interaction(codes, drop = TRUE, lex.order = lex.order)
  first <- match(levels(key), as.character(key))
  display <- vapply(first, function(i)
    paste(vapply(data[cols], function(x) as.character(x[i]), character(1)),
          collapse = sep), character(1))
  levels(key) <- make.unique(display)
  key
}

# Build the per-level networks. `data` is the original (data-frame) input;
# `group` names the grouping column(s); the remaining arguments are forwarded to
# psychnet() unchanged for each level (with the group column removed and any
# `vars` selection already resolved to column names).
#' @noRd
.psychnet_grouped <- function(data, group, source, vars_q, vars_env, method,
                              actor, action, session, time, compute_sessions,
                              time_threshold, id, standardize, threshold,
                              gamma, labels, dots, within = NULL) {
  if (!is.data.frame(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
  miss <- setdiff(group, names(data))
  if (length(miss))
    .psn_stop("Group column(s) not found in data: ", paste(miss, collapse = ", "),
         call. = FALSE)

  # Level order decides group 1 and 2 of net_compare() and the order of every
  # per-group result. A factor grouping keeps its level order (several factor
  # columns: their levels in lexical tuple order); otherwise levels follow
  # their order of appearance.
  all_factor <- all(vapply(data[group], is.factor, logical(1)))
  key_f <- if (length(group) == 1L) NULL
           else .tuple_factor(data, group, sep = "-", lex.order = all_factor)
  key <- if (length(group) == 1L) as.character(data[[group]])
         else as.character(key_f)
  ok   <- !is.na(key)
  if (any(!ok)) {
    .psn_warn(sprintf(
      "%d of %d rows have a missing `group` value and were dropped.",
      sum(!ok), length(ok)), class = "psychnet_rows_dropped")
  }
  if (any(key[ok] == ""))
    .psn_stop("`group` has an empty-string level; give those rows a name or ",
              "set them to NA.")
  levs <- if (!all_factor) unique(key[ok])
          else if (length(group) == 1L) levels(droplevels(data[[group]][ok]))
          else levels(droplevels(key_f[ok]))
  if (length(levs) < 2L)
    .psn_stop("`group` must have at least 2 levels; found ", length(levs), ".",
         call. = FALSE)

  # Resolve a vars selection once (against the full data), then drop the group
  # column(s) from it so the grouping variable is never itself a node.
  vars_sel <- if (is.null(vars_q)) NULL
              else setdiff(.select_names(data, vars_q, vars_env), c(group, within))
  cl <- list(method = method, threshold = threshold, gamma = gamma,
             labels = labels, dots = dots)
  rows_of <- function(lv) data[ok & key == lv, setdiff(names(data), group),
                               drop = FALSE]

  if (source == "eventdata") {
    # Reduce each level's event log to its modelling matrix, fit that, and keep
    # the matrix as the bootstrappable subset -- so an event-data group resamples
    # exactly like a cross-sectional one (no temporal model). A raw log ignores
    # `vars` (actions are the variables); a feature table honors it.
    raw_log <- !is.null(actor)
    des <- lapply(levs, function(lv)
      .event_design(rows_of(lv), actor, action, session, time,
                    compute_sessions, time_threshold, id,
                    if (raw_log) NULL else vars_sel, standardize, labels,
                    na_method = dots$na_method %||% "pairwise"))
    nets <- lapply(des, function(d)
      do.call(.fit_design, c(list(d, method = method, threshold = threshold,
                                  gamma = gamma), dots)))
    # A level reduces to a single design matrix unless it is a within/between
    # split (standardize = FALSE with several occasions); only all-single groups
    # are bootstrappable, so mark the container source accordingly.
    subsets <- lapply(des, function(d) d$matrix)
    # The ordinary resamplers draw independent rows. Centered repeated-person
    # rows remain dependent and their effective n is N-J, so do not advertise
    # these designs as cross-sectional bootstrappable subsets.
    independent <- vapply(des, function(d)
      !is.null(d$matrix) && identical(as.integer(d$n), nrow(d$matrix)), logical(1))
    src <- if (all(independent)) "data" else "eventdata"
  } else {
    # Cross-sectional table: keep the selected columns and re-estimate per level.
    cols_of <- function(lv) {
      sub <- rows_of(lv)
      if (!is.null(within)) sub[[within]] <- NULL       # the person is not a node
      if (is.null(vars_sel)) sub else sub[, intersect(vars_sel, names(sub)),
                                          drop = FALSE]
    }
    subsets <- lapply(levs, cols_of)
    # Which person each row is, per level: lets net_compare(paired = TRUE)
    # pair two levels by person rather than by row position.
    persons <- if (is.null(within)) NULL
               else stats::setNames(lapply(levs, function(lv) rows_of(lv)[[within]]),
                                    levs)
    .check_group_sizes(subsets, levs)
    nets <- lapply(subsets, function(s)
      do.call(psychnet, c(list(data = s, method = method, threshold = threshold,
                               gamma = gamma, labels = labels), dots)))
    src <- "data"
  }
  if (source != "eventdata") {
    # Each network remembers its level, rows, persons and fitting settings, so
    # net_compare(nets$A, nets$B) can compare two members directly.
    nets <- lapply(seq_along(levs), function(k)
      structure(nets[[k]], group_member = list(
        level = levs[k], rows = subsets[[k]],
        persons = if (is.null(within)) NULL else persons[[k]],
        within = within, call = cl)))
  }
  names(nets) <- levs
  names(subsets) <- levs
  out <- .new_psychnet_group(nets, group_col = group, source = src,
                             subsets = subsets, call = cl)
  if (source != "eventdata" && !is.null(within)) {
    attr(out, "persons") <- persons
    attr(out, "within") <- within
  }
  out
}

# A network is not estimable from fewer than 3 rows, and is saturated -- or
# singular -- from no more rows than variables: refuse the first, warn on the
# second, naming the levels.
#' @noRd
.check_group_sizes <- function(subsets, levs) {
  rows <- vapply(subsets, nrow, integer(1))
  vars <- vapply(subsets, function(s) sum(vapply(s, is.numeric, logical(1))),
                 integer(1))
  tiny <- rows < 3L
  if (any(tiny))
    .psn_stop(sprintf(
      "Every group level needs at least 3 rows to estimate a network; %s.",
      paste(sprintf("level \"%s\" has %d", levs[tiny], rows[tiny]),
            collapse = ", ")))
  small <- rows <= vars
  if (any(small))
    .psn_warn(sprintf(paste0(
      "Group level(s) with no more rows than variables: %s. Their networks ",
      "rest on a saturated or singular correlation matrix and are unreliable."),
      paste(sprintf("\"%s\" (%d rows, %d variables)", levs[small],
                    rows[small], vars[small]), collapse = ", ")),
      class = "psychnet_small_sample")
  invisible(NULL)
}

#' Per-group psychometric networks
#'
#' The object returned by [psychnet()] when `group` is supplied: a named list of
#' [psychnet] networks, one per level of the grouping variable. It plots as a
#' grid with `cograph::splot()` and is consumed per level by the framework verbs
#' ([net_centralities()], [net_predict()], [net_boot()], [net_stability()],
#' [net_compare()]), each returning a `*_group` result.
#' Data-based resampling and comparison require independent cross-sectional
#' rows. Groups built from several occasions per actor need resampling at the
#' actor level and are refused by these row-based diagnostics.
#'
#' Networks are listed in the grouping variable's level order when it is a
#' factor (so its first level is group 1 of [net_compare()]), otherwise in
#' order of first appearance. Rows with a missing group value are dropped
#' with a `psychnet_rows_dropped` warning; an empty-string level, or a
#' cross-sectional level with fewer than 3 rows, is a `psychnet_bad_input`
#' error; a level with no more rows than variables raises a
#' `psychnet_small_sample` warning.
#'
#' @param x A `psychnet_group` object.
#' @param ... Ignored.
#' @return `x`, invisibly (for `print`), after printing one row per group with
#'   its node and edge counts and sample size.
#' @examples
#' d <- rbind(data.frame(SRL_GPT, source = "GPT"),
#'            data.frame(SRL_Claude, source = "Claude"))
#' nets <- psychnet(d, method = "glasso", group = "source")
#' nets
#' summary(nets)
#' as.data.frame(nets)
#' @name psychnet_group
#' @aliases psychnet_group
#' @export
print.psychnet_group <- function(x, ...) {
  gc <- paste(attr(x, "group_col"), collapse = ", ")
  cat(sprintf("<psychnet_group> %d networks by %s (method: %s)\n",
              length(x), gc, attr(x, "call")$method))
  tab <- .group_network_summary(x)
  tab$mean_abs_weight <- NULL
  print(tab, row.names = FALSE)
  invisible(x)
}

#' @rdname psychnet_group
#' @param row.names,optional Ignored (S3 consistency).
#' @return For `as.data.frame`, the per-group edge lists stacked with a leading
#'   `group` column: one row per edge per group, columns `group`, `from`,
#'   `to`, `weight`.
#' @export
as.data.frame.psychnet_group <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  .stack_by_group(lapply(x, as.data.frame), names(x))
}

#' @rdname psychnet_group
#' @param object A `psychnet_group` object.
#' @return For `summary`, a data frame with node/edge counts and mean absolute
#'   edge weight: one row per group, or per group and `level` when a group
#'   contains within/between networks.
#' @export
summary.psychnet_group <- function(object, ...) {
  tab <- .group_network_summary(object)
  tab$n <- NULL
  tab
}

# A group may hold a within/between pair; keep both levels visible in its
# summary instead of treating the pair as if it were a single network.
#' @noRd
.group_network_summary <- function(x) {
  nested <- any(vapply(x, inherits, logical(1), "psychnet_multilevel"))
  parts <- lapply(seq_along(x), function(i) {
    nets <- if (inherits(x[[i]], "psychnet_multilevel")) x[[i]] else list(x[[i]])
    lev <- if (inherits(x[[i]], "psychnet_multilevel")) names(nets) else "network"
    rows <- lapply(seq_along(nets), function(j) {
      net <- nets[[j]]
      row <- data.frame(group = names(x)[i], nodes = nrow(net$nodes),
        edges = nrow(net$edges), n = as.integer(net$n %||% NA_integer_),
        mean_abs_weight = if (nrow(net$edges)) mean(abs(net$edges$weight)) else 0,
        stringsAsFactors = FALSE)
      if (nested) row <- cbind(row["group"], level = lev[j], row[-1L])
      row
    })
    do.call(rbind, rows)
  })
  do.call(rbind, parts)
}

# ---- generic result-group plumbing (shared by every verb's group result) -----

`%||%` <- function(a, b) if (is.null(a)) b else a

# Stack a list of per-group data frames into one, prefixing a `group` column.
#' @noRd
.stack_by_group <- function(dfs, levels) {
  parts <- Map(function(df, lv) {
    if (is.null(df) || !nrow(df))
      return(data.frame(group = character(0), df))
    data.frame(group = rep(lv, nrow(df)), df, row.names = NULL,
               stringsAsFactors = FALSE)
  }, dfs, levels)
  out <- do.call(rbind, unname(parts))   # unnamed: rbind would build "1.1" row names
  rownames(out) <- NULL
  out
}

# Object-first verbs (centrality, predict): apply `fn` to each fitted network.
#' @noRd
# Object-first group verbs (net_centralities, net_predict) operate per network;
# a within/between (psychnet_multilevel) level has two networks, not one, so
# reject it with a clear message instead of crashing inside the bare-matrix path.
.reject_multilevel_group <- function(x, verb) {
  if (any(vapply(x, inherits, logical(1), "psychnet_multilevel")))
    .psn_stop(sprintf(paste0("%s() does not support within/between (multilevel) ",
         "groups. Rebuild the group with standardize = TRUE for one ",
         "person-centered network per level, or read the within and between ",
         "edges with as.data.frame()."), verb))
  invisible(NULL)
}

.group_obj_apply <- function(x, fn, out_class, ...) {
  res <- lapply(x, fn, ...)
  names(res) <- names(x)
  attr(res, "group_col") <- attr(x, "group_col")
  class(res) <- c(out_class, "psychnet_result_group")
  res
}

# Data-first group dispatch shared by the resampling verbs (net_boot,
# net_stability, net_casedrop_reliability, net_split_reliability): re-run `verb` on
# each level's stored cross-sectional subset, reproducing the same estimator
# configuration the group networks were built with through estimator_args.
# `args` carries the verb's own diagnostic options.
.group_data_apply <- function(data, verb, verb_name, out_class, args) {
  if (!identical(attr(data, "source"), "data"))
    .psn_stop(sprintf("%s() supports group mode for cross-sectional data only.",
                 verb_name), call. = FALSE)
  subs <- attr(data, "subsets"); cl <- attr(data, "call")
  estimator_args <- cl$dots
  if (is.null(estimator_args)) estimator_args <- list()
  if (!is.null(cl$gamma)) estimator_args$gamma <- cl$gamma
  threshold_method <- cl$method %in% c("cor", "pcor", "glasso", "huge",
                                        "ggm", "logo") ||
    (identical(cl$method, "mgm") && is.character(cl$threshold))
  if (threshold_method && !is.null(cl$threshold))
    estimator_args$threshold <- cl$threshold
  supplied <- args$estimator_args
  if (length(supplied)) estimator_args[names(supplied)] <- supplied
  args$estimator_args <- estimator_args
  fixed <- c(list(method = cl$method, labels = cl$labels), args)
  res <- lapply(names(subs), function(lv)
    do.call(verb, c(list(subs[[lv]]), fixed)))
  names(res) <- names(subs)
  attr(res, "group_col") <- attr(data, "group_col")
  class(res) <- c(out_class, "psychnet_result_group")
  res
}

#' Per-group framework results
#'
#' The list returned by a framework verb ([net_centralities()], [net_predict()],
#' [net_boot()], [net_stability()]) applied to a [psychnet_group]: one result per
#' group level, keyed by level. `as.data.frame()` stacks the per-group tables
#' with a leading `group` column.
#'
#' @param x A `psychnet_result_group` object.
#' @param ... Ignored.
#' @return `x`, invisibly (for `print`), after printing each group's result
#'   under a group heading.
#' @examples
#' d <- rbind(data.frame(SRL_GPT, source = "GPT"),
#'            data.frame(SRL_Claude, source = "Claude"))
#' nets <- psychnet(d, method = "glasso", group = "source")
#' cent <- net_centralities(nets)
#' cent
#' as.data.frame(cent)
#' @name psychnet_result_group
#' @aliases psychnet_result_group
#' @export
print.psychnet_result_group <- function(x, ...) {
  cat(sprintf("<%s> %d groups: %s\n", class(x)[1L], length(x),
              paste(names(x), collapse = ", ")))
  lapply(names(x), function(nm) {
    cat(sprintf("\n--- %s ---\n", nm))
    print(x[[nm]])
  })
  invisible(x)
}

#' @rdname psychnet_result_group
#' @param row.names,optional Ignored (S3 consistency).
#' @return For `as.data.frame`, the per-group tables stacked with a `group`
#'   column: each group's result converted by its own `as.data.frame()`
#'   method, with `...` (e.g. `what =`) passed on to it. An error converting
#'   any group's result is raised, not dropped.
#' @export
as.data.frame.psychnet_result_group <- function(x, row.names = NULL,
                                                optional = FALSE, ...) {
  dfs <- lapply(x, function(el) as.data.frame(el, ...))
  .stack_by_group(dfs, names(x))
}
