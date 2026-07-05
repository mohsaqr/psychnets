# Networks from actor / time / event data, clean-room base R. Mirrors the
# Nestimate event pipeline -- prepare the log (actor, action, time, session),
# convert to action frequencies, then estimate -- reimplemented with no compiled
# dependency.
#
# The structure of the data decides the model:
#   * one occasion per actor (only `actor`)                  -> a single GGM on
#     the actor x action frequency matrix (cross-sectional network);
#   * several occasions per actor (explicit `session`, or sessions COMPUTED from
#     time gaps) -> the observations are nested in actors, so the covariance is
#     decomposed into a WITHIN-actor and a BETWEEN-actor network -- or, with
#     `standardize = TRUE`, the clustering is removed by person-centering and a
#     single within network is fit.
#
# A `time` column is used only to COMPUTE sessions from gaps (no temporal /
# lagged modelling): `compute_sessions = TRUE` starts a new session whenever the
# gap between consecutive events exceeds `time_threshold`.

# Coerce a time column to numeric seconds for gap computation: numeric as-is,
# date/time via as.numeric, character parsed as a timestamp (ISO8601 etc.).
#' @noRd
.as_time_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.numeric(x))
  parsed <- suppressWarnings(as.numeric(as.POSIXct(as.character(x), tz = "UTC")))
  if (all(is.na(parsed))) {
    stop("`time` could not be parsed as numeric or a timestamp.", call. = FALSE)
  }
  parsed
}

# Long event log -> (occasions x actions) count matrix + the actor of each
# occasion. An occasion is one (actor[, session]) cell; with compute_sessions it
# is one (actor[, session], inferred-session) cell. Returns features + actor so
# .net_eventdata can cluster occasions by actor without column juggling.
#' @noRd
.event_freq <- function(data, actor, action, session, time,
                        compute_sessions, time_threshold) {
  acts <- as.character(data[[action]])
  keep <- !is.na(acts) & acts != ""
  data <- data[keep, , drop = FALSE]; acts <- acts[keep]
  levels_a <- sort(unique(acts))
  if (length(levels_a) < 2L) {
    stop("`action` has fewer than 2 distinct values; nothing to relate.",
         call. = FALSE)
  }
  base_group <- if (length(c(actor, session)) == 1L) factor(data[[actor]])
                else interaction(data[c(actor, session)], drop = TRUE, lex.order = TRUE)

  if (isTRUE(compute_sessions) && !is.null(time)) {
    tnum <- .as_time_numeric(data[[time]])
    ord  <- order(base_group, tnum)
    data <- data[ord, , drop = FALSE]; acts <- acts[ord]
    base_group <- base_group[ord, drop = TRUE]; tnum <- tnum[ord]
    # --- session detection copied from Nestimate::prepare() (time gaps) ---
    inferred <- stats::ave(tnum, base_group, FUN = function(t) {
      gaps <- c(NA_real_, diff(t))
      new_session <- is.na(gaps) | gaps > time_threshold
      cumsum(new_session)
    })
    # ---------------------------------------------------------------------
    occ <- interaction(base_group, inferred, drop = TRUE, lex.order = TRUE)
  } else {
    occ <- base_group
  }
  occ <- factor(occ)

  counts <- unclass(table(occ, factor(acts, levels = levels_a)))
  feats  <- matrix(as.numeric(counts), nrow(counts), ncol(counts),
                   dimnames = list(NULL, levels_a))
  actor_key <- if (length(actor) == 1L) as.character(data[[actor]])
               else as.character(interaction(data[actor], drop = TRUE))
  actor_per_occ <- vapply(split(actor_key, occ), `[`, character(1), 1L)
  list(features = feats, actor = unname(actor_per_occ[levels(occ)]),
       occasion = levels(occ))
}

#' Action frequencies from an event log
#'
#' Converts a long event log -- one row per event, with an actor, an action, and
#' optionally a session and time -- into a frequency table: one row per occasion,
#' one column per distinct action, holding that action's count. The conversion
#' mirrors the "frequency" format of the Nestimate event pipeline. An occasion is
#' one (actor, session) cell; with `compute_sessions = TRUE` the `time` column is
#' used to split each actor into sessions wherever the gap between consecutive
#' events exceeds `time_threshold`. This is the frequency input that
#' [psychnet()] builds internally for event data (`source = "eventdata"`).
#'
#' @param data A long event log (data frame), one row per event.
#' @param actor Column(s) naming the actor / subject. Default `"Actor"`.
#' @param action Column naming the action / state. Default `"Action"`.
#' @param session Optional column(s) naming an explicit session within an actor.
#' @param time Optional column used to compute sessions from gaps (see
#'   `compute_sessions`); it is not used for any temporal model.
#' @param compute_sessions If `TRUE`, split each actor into sessions from the
#'   `time` gaps (a new session starts when the gap exceeds `time_threshold`).
#'   Default `TRUE`.
#' @param time_threshold Maximum gap (in the units of `time`, seconds for a
#'   timestamp) between consecutive events before a new session begins. Default
#'   `900` (15 minutes), as in Nestimate.
#' @return A `data.frame` with an `actor` column, a `session` index, and one
#'   integer count column per action (one row per occasion).
#' @examples
#' ev <- data.frame(
#'   Actor   = rep(c("a", "b"), each = 6),
#'   Session = rep(rep(1:2, each = 3), 2),
#'   Action  = c("read","quiz","read", "quiz","read","note",
#'               "note","note","read", "read","quiz","quiz"))
#' event_frequencies(ev, session = "Session")
#' @export
event_frequencies <- function(data, actor = "Actor", action = "Action",
                              session = NULL, time = NULL,
                              compute_sessions = TRUE, time_threshold = 900) {
  stopifnot(is.data.frame(data), action %in% names(data))
  miss <- setdiff(c(actor, session, time), names(data))
  if (length(miss)) {
    stop("Columns not found in data: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  ef <- .event_freq(data, actor, action, session, time,
                    compute_sessions, time_threshold)
  data.frame(actor = ef$actor, session = ef$occasion, ef$features,
             row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
}

# Split a numeric matrix into between-actor (per-actor means) and within-actor
# (deviations from each actor's mean) components, with n_between = #actors and
# n_within = #rows - #actors (pooled within df).
#' @noRd
.within_between_decompose <- function(V, ids) {
  f  <- factor(ids)
  J  <- nlevels(f)
  if (J < 2L) stop("Need at least 2 actors for a between-actor network.",
                   call. = FALSE)
  nj <- tabulate(f)
  means  <- rowsum(V, f) / nj
  within <- V - means[as.integer(f), , drop = FALSE]
  list(between = means, within = within,
       n_between = J, n_within = nrow(V) - J)
}

# Resolve and call a GGM estimator on a correlation matrix at sample size n.
#' @noRd
.fit_ggm <- function(mat, n, method, ...) {
  fn <- switch(method, glasso = ebic_glasso, pcor = pcor_network,
               ggm = ggm_modselect, huge = huge_network,
               tmfg = tmfg_network, logo = logo_network,
               stop(sprintf("Unsupported method '%s' for net_multilevel().", method),
                    call. = FALSE))
  # A feature with no variance at this level (e.g. an action that is constant
  # across occasions, or constant within every actor) makes the correlation
  # undefined; drop it rather than fail deep in the estimator.
  sds <- apply(mat, 2L, stats::sd)
  bad <- !is.finite(sds) | sds < 1e-12
  if (any(bad)) {
    warning(sprintf("Dropping %d zero-variance feature(s) at this level: %s.",
                    sum(bad), paste(colnames(mat)[bad], collapse = ", ")),
            call. = FALSE)
    mat <- mat[, !bad, drop = FALSE]
  }
  if (ncol(mat) < 2L) {
    stop("Fewer than 2 features with variance remain to relate.", call. = FALSE)
  }
  # Build the association matrix the way the requested estimator would, so
  # data-level options forwarded through `...` are honored: cor_method
  # (spearman/kendall/auto) for the correlation-based GGMs, and the
  # nonparanormal transform (npn) for `huge`. The level's effective sample size
  # (n_within / n_between, or n) is preserved by passing it explicitly alongside
  # the precomputed matrix rather than letting the estimator recount rows.
  dots <- list(...)
  cor_method <- dots$cor_method %||% "pearson"
  na_method  <- dots$na_method  %||% "pairwise"
  S <- if (method == "huge")
    .npn_cor(mat, dots$npn %||% "shrinkage", na_method)$S
  else
    .cor_input(mat, method = cor_method, na_method = na_method)$S
  dots$cor_method <- NULL; dots$na_method <- NULL; dots$npn <- NULL
  args <- c(list(cor_matrix = S), dots)
  # Only the estimators with an `n` formal take a sample size (TMFG does not).
  if ("n" %in% names(formals(fn)) && !"n" %in% names(args)) args$n <- n
  do.call(fn, args)
}

# Build the numeric modelling matrix from an event log / feature table. Returns a
# list describing the design the network is fit on:
#   * one occasion per actor, or standardize=TRUE -> $matrix (+ $n): one matrix
#     (raw frequencies, or person-centered deviations) fit as a single GGM;
#   * several occasions per actor & standardize=FALSE -> $within / $between
#     (+ counts): the covariance split into within- and between-actor parts.
# Keeping this separate from the fit lets group mode reuse the very matrix the
# network is built on, so event-data groups bootstrap like cross-sectional ones.
#' @noRd
.event_design <- function(data, actor, action, session, time,
                          compute_sessions, time_threshold, id, vars,
                          standardize, labels) {
  if (!is.null(actor)) {                                # event-log mode
    miss <- setdiff(c(actor, session, time), names(data))
    if (length(miss)) {
      stop("Columns not found in data: ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    ef  <- .event_freq(data, actor, action, session, time,
                       compute_sessions, time_threshold)
    V   <- ef$features
    ids <- ef$actor
  } else {                                              # numeric feature table
    if (is.null(id) || !all(id %in% names(data))) {
      stop("Supply `actor` for an event log, or `id` naming the actor column of a feature table.",
           call. = FALSE)
    }
    if (is.null(vars)) {
      vars <- setdiff(names(data), id)
      vars <- vars[vapply(data[vars], is.numeric, logical(1))]
    }
    V   <- as.matrix(data[, vars, drop = FALSE])
    ids <- if (length(id) == 1L) data[[id]] else interaction(data[id], drop = TRUE)
  }
  if (ncol(V) < 2L) stop("Need at least 2 feature columns.", call. = FALSE)

  # Honor user labels: rename the feature columns so the labels flow through the
  # fit (and survive any zero-variance column drop in .fit_ggm).
  if (!is.null(labels)) {
    if (length(labels) != ncol(V)) {
      stop(sprintf("`labels` length (%d) must match the number of features (%d).",
                   length(labels), ncol(V)), call. = FALSE)
    }
    colnames(V) <- labels
  }

  occ_per_actor <- tabulate(factor(ids))
  if (max(occ_per_actor) == 1L) return(list(matrix = V, n = nrow(V)))
  dec <- .within_between_decompose(V, ids)
  if (isTRUE(standardize)) return(list(matrix = dec$within, n = dec$n_within))
  list(within = dec$within, between = dec$between,
       n_within = dec$n_within, n_between = dec$n_between, n_obs = nrow(V))
}

# Fit a GGM (or within/between pair) from an .event_design() result. `method` /
# `threshold` / `gamma` are already resolved by psychnet(); `...` carries extra
# estimator arguments (e.g. cor_method).
#' @noRd
.fit_design <- function(des, method, threshold, gamma, ...) {
  extra <- list(...)
  if (method %in% c("pcor", "glasso", "huge", "ggm", "logo"))
    extra$threshold <- threshold
  if (!is.null(gamma) && method %in% c("glasso", "huge", "ggm"))
    extra$gamma <- gamma
  fit <- function(mat, n) do.call(.fit_ggm, c(list(mat, n, method), extra))
  if (!is.null(des$matrix)) return(fit(des$matrix, des$n))
  # A two-network group (within + between): the list elements are the only
  # contents so cograph::splot() lays them out as a grid (it dispatches on
  # "netobject_group", exactly like psychnet_group); run-level metadata lives in
  # attributes, not as list elements.
  out <- list(within  = fit(des$within,  des$n_within),
              between = fit(des$between, des$n_between))
  attr(out, "method") <- method
  attr(out, "n_actors") <- des$n_between
  attr(out, "n_obs") <- des$n_obs
  attr(out, "group_col") <- "level"
  class(out) <- c("psychnet_multilevel", "netobject_group")
  out
}

# Networks from actor / event data, dispatched from psychnet(source="eventdata").
# Adapts to the data's structure: one occasion per actor -> a single GGM; several
# occasions per actor -> a person-centered single network (standardize=TRUE) or a
# within/between pair (standardize=FALSE). A `time` column is used only to compute
# sessions; there is no temporal model.
#' @noRd
.net_eventdata <- function(data, actor = NULL, action = "Action",
                           session = NULL, time = NULL,
                           compute_sessions = TRUE, time_threshold = 900,
                           id = NULL, vars = NULL, standardize = TRUE,
                           method = "glasso", threshold = 0, gamma = NULL,
                           labels = NULL, ...) {
  stopifnot(is.data.frame(data), is.logical(standardize), length(standardize) == 1L)
  if (!method %in% c("glasso", "pcor", "ggm", "huge", "tmfg", "logo")) {
    stop(sprintf("Event data supports the Gaussian graphical methods only; '%s' is not one of them.",
                 method), call. = FALSE)
  }
  des <- .event_design(data, actor, action, session, time, compute_sessions,
                       time_threshold, id, vars, standardize, labels)
  .fit_design(des, method, threshold, gamma, ...)
}

#' Within / between event-data networks
#'
#' The object returned by [psychnet()] for nested event data with
#' `standardize = FALSE`: a pair of Gaussian graphical networks decomposing the
#' actor-by-action frequency covariance into a within-actor and a between-actor
#' part.
#'
#' @param x A `psychnet_multilevel` object.
#' @param ... Ignored.
#' @return `x`, invisibly (for `print`).
#' @export
print.psychnet_multilevel <- function(x, ...) {
  cat(sprintf("<psychnet_multilevel> %s networks (within + between)\n",
              attr(x, "method")))
  cat(sprintf("  actors: %d   occasions: %d\n",
              attr(x, "n_actors"), attr(x, "n_obs")))
  cat(sprintf("  within : %d nodes, %d edges\n",
              nrow(x$within$nodes), nrow(x$within$edges)))
  cat(sprintf("  between: %d nodes, %d edges\n",
              nrow(x$between$nodes), nrow(x$between$edges)))
  invisible(x)
}

#' @rdname print.psychnet_multilevel
#' @param row.names,optional Ignored (S3 consistency).
#' @return For `as.data.frame`, the two edge lists stacked with a `level` column.
#' @export
as.data.frame.psychnet_multilevel <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  tag <- function(df, lv) data.frame(level = rep(lv, nrow(df)), df,
                                      row.names = NULL, stringsAsFactors = FALSE)
  rbind(tag(as.data.frame(x$within), "within"),
        tag(as.data.frame(x$between), "between"))
}
