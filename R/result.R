# The shared result object returned by every psychnet estimator. It mirrors the
# Nestimate / cograph "netobject" so a fitted psychnet network is plottable by
# cograph::splot() and comparable by Nestimate, while staying lean: it carries
# only the estimated content (the weighted network, the tidy node and edge
# tables, the sample size, the estimator name, and method-specific extras such
# as the precision matrix, the EBIC penalty, and the KKT certificate).

# Tidy one-row-per-edge table from a weighted adjacency matrix, with LABEL
# from/to (the user-facing form returned by as.data.frame()).
#' @noRd
.edges_df <- function(weights, labels, directed, include_zero = FALSE) {
  idx <- if (isTRUE(directed)) which(row(weights) != col(weights), arr.ind = TRUE)
         else which(upper.tri(weights), arr.ind = TRUE)
  w <- weights[idx]
  keep <- if (include_zero) rep(TRUE, length(w)) else abs(w) > 1e-12
  data.frame(from = labels[idx[keep, 1L]], to = labels[idx[keep, 2L]],
             weight = w[keep], stringsAsFactors = FALSE, row.names = NULL)
}

# Edge table with INTEGER node-index from/to -- the format cograph's
# cograph_network requires, so the object plots directly via cograph::splot().
#' @noRd
.edges_index_df <- function(weights, directed, include_zero = FALSE) {
  idx <- if (isTRUE(directed)) which(row(weights) != col(weights), arr.ind = TRUE)
         else which(upper.tri(weights), arr.ind = TRUE)
  w <- weights[idx]
  keep <- if (include_zero) rep(TRUE, length(w)) else abs(w) > 1e-12
  data.frame(from = as.integer(idx[keep, 1L]), to = as.integer(idx[keep, 2L]),
             weight = w[keep], row.names = NULL)
}

#' Construct a psychnet network object
#'
#' @param graph p x p weighted adjacency matrix (partial correlations for the
#'   Gaussian models, regression-weight matrix for Ising/mgm).
#' @param labels Character node labels.
#' @param method Estimator name.
#' @param directed Logical; TRUE only for inherently directed estimators.
#' @param n_obs Sample size used.
#' @param extra Named list of method-specific fields (e.g. precision, lambda).
#' @param data Fit data, passed through only so nodewise predictability (ising /
#'   mgm) can be computed and stored on the node table at fit time. `NULL` for
#'   the GGM / structural methods, whose predictability is closed-form or stored.
#' @return An object of class `c("psychnet", "cograph_network")`.
#' @noRd
.new_psychnet <- function(graph, labels, method, directed, n_obs,
                          extra = list(), data = NULL) {
  if (length(labels) != ncol(graph)) {
    .psn_stop(sprintf(paste0("labels length (%d) does not match the network ",
                        "dimension (%d); a non-numeric or zero-variance column ",
                        "may have been dropped from the data."),
                 length(labels), ncol(graph)), call. = FALSE)
  }
  if (anyNA(labels) || anyDuplicated(labels))
    .psn_stop("`labels` must be unique and non-missing.", call. = FALSE)
  weights <- graph
  diag(weights) <- 0
  dimnames(weights) <- list(labels, labels)
  # Build a genuine cograph_network: nodes carry x/y layout slots, edges use
  # integer node indices, and meta/data/node_groups are present -- so the object
  # plots directly with cograph::splot(net) (not splot(net$weights)).
  source_names <- extra$nodewise$source_names
  if (is.null(source_names) && !is.null(data) && ncol(data) == length(labels))
    source_names <- colnames(data)
  if (is.null(source_names)) source_names <- labels
  if (!is.null(extra$nodewise)) extra$nodewise$source_names <- source_names
  nodes <- data.frame(id = seq_along(labels), label = labels, name = source_names,
                      x = NA_real_, y = NA_real_, stringsAsFactors = FALSE)
  edges <- .edges_index_df(weights, directed)
  net <- structure(
    c(list(nodes = nodes, edges = edges, directed = directed, weights = weights,
           data = NULL, meta = list(source = "psychnet"), node_groups = NULL,
           method = method, n = n_obs), extra),
    class = c("psychnet", "cograph_network")
  )
  # Store node predictability + the default-draw flag so cograph::splot() can
  # render the predictability ring straight from the object.
  .attach_predictability(net, data = data)
}

#' Fitted psychometric network
#'
#' @description
#' A `psychnet` object holds a fitted psychometric network; every estimator in
#' the package returns its network in this form. The network consists of
#' nodes, which are the variables,
#' and weighted edges, which are the estimated associations: partial
#' correlations for the Gaussian graphical models, correlations for
#' [cor_network()], symmetrized regression weights for the Ising and mixed
#' graphical models, and directed importance shares for [relimp_network()].
#'
#' `print()` reports the estimator, the numbers of nodes and edges, the
#' selected penalty and the certificate; `summary()` adds the range of the edge
#' weights; `as.data.frame()` returns the edges as a table; and `plot()` draws
#' the network with the `cograph` package. The certificate is read with
#' [certificate()].
#'
#' Analyses of a fitted network, such as [net_centralities()], [net_boot()]
#' and [net_stability()], return their own result classes, documented on their
#' pages. Networks estimated per group are collected in a [psychnet_group].
#'
#' @details
#' The object is a list that follows the network object of the `cograph` and
#' `Nestimate` packages, so it can be passed to `cograph::splot()` directly.
#' It stores the weighted adjacency matrix with a zero diagonal (`weights`),
#' the node table (`nodes`, with `id`, `label` and `name`), the edge table in
#' the cograph format with integer node indices (`edges`), whether the network
#' is directed (`directed`), the estimator (`method`) and the sample size
#' (`n`), together with fields specific to each estimator, such as the
#' precision matrix (`precision`), the selected penalty (`lambda`), the EBIC
#' hyperparameter (`gamma`), the correlation matrix (`cor_matrix`), the
#' certificate (`kkt`), the Ising thresholds (`thresholds`) and the node types
#' of a mixed model (`types`). Field access with `$` also recognises the names
#' used by earlier versions: `graph` (the weights), `labels`, `n_nodes`,
#' `n_edges` and `n_obs`.
#'
#' @param x,object A `psychnet` object.
#' @param name For `$`: a field name, or one of the earlier names listed in
#'   Details.
#' @param include_zero For `as.data.frame()`: keep absent (zero-weight) edges.
#'   Default `FALSE`.
#' @param row.names,optional Ignored.
#' @param ... For `plot()`, passed to `cograph::splot()`; ignored by the other
#'   methods.
#' @return `as.data.frame()` returns a data frame with one row per edge
#'   (each node pair once for an undirected network, each ordered pair for a
#'   directed one) and columns `from`, `to` (node labels) and `weight`.
#'   `print()` returns `x` invisibly. `summary()` prints the network with its
#'   edge-weight range and returns the edge table invisibly. `plot()` returns
#'   the value of `cograph::splot()` invisibly. `$` returns the requested
#'   field, or `NULL` when the name is neither a field nor an earlier name.
#' @section Conditions:
#'   `psychnet_bad_input` from `plot()` when the `cograph` package is not
#'   installed.
#' @seealso [psychnet()] and the estimators it routes to; [certificate()];
#'   [psychnet_group] for networks per group.
#' @examples
#' fit <- pcor_network(SRL_GPT)
#' fit
#' summary(fit)
#' as.data.frame(fit)
#' as.data.frame(fit, include_zero = TRUE)
#' if (requireNamespace("cograph", quietly = TRUE)) plot(fit)
#' @name psychnet-class
#' @aliases psychnet-class
NULL

#' @describeIn psychnet-class Field access, including the legacy aliases.
#' @method $ psychnet
#' @export
`$.psychnet` <- function(x, name) {
  fld <- .subset2(x, name)
  if (!is.null(fld)) return(fld)
  switch(name,
    graph   = .subset2(x, "weights"),
    labels  = .subset2(x, "nodes")[["label"]],
    n_nodes = nrow(.subset2(x, "nodes")),
    n_edges = nrow(.subset2(x, "edges")),
    n_obs   = .subset2(x, "n"),
    NULL)
}

#' @describeIn psychnet-class The estimator, the node and edge counts, the
#'   selected penalty and the certificate (see [certificate()]).
#' @export
print.psychnet <- function(x, ...) {
  cat(sprintf("<psychnet> %s network\n", x$method))
  cat(sprintf("  nodes: %d   edges: %d   (%s)\n",
              nrow(x$nodes), nrow(x$edges),
              if (isTRUE(x$directed)) "directed" else "undirected"))
  if (!is.null(x$lambda)) {
    cat(sprintf("  lambda: %.4g   gamma: %.2g\n", x$lambda, x$gamma))
  }
  # Structural certificates (TMFG, relimp) check a defining identity, not an
  # optimum, so they must not be labelled as a KKT residual.
  if (identical(.certificate_kind(x$method), "structural")) {
    value <- if (identical(x$method, "tmfg") && !is.null(x$adjacency))
      tmfg_certificate(x) else x$kkt
    if (!is.null(value)) {
      cat(sprintf("  structural certificate (%s): %s\n",
                  if (identical(x$method, "tmfg"))
                    "violated TMFG invariants" else "decomposition residual",
                  if (identical(x$method, "tmfg")) format(value)
                  else sprintf("%.2e", value)))
    }
  } else if (!is.null(x$kkt)) {
    if (is.finite(x$kkt)) {
      cat(sprintf("  optimality (KKT residual): %.2e\n", x$kkt))
    } else if (identical(x$certificate_target, "pre_threshold_fit")) {
      cat("  optimality: unavailable for post-threshold graph",
          sprintf(" (pre-threshold residual %.2e)\n", x$fit_kkt))
    } else {
      cat("  optimality: unavailable\n")
    }
  }
  invisible(x)
}

#' @describeIn psychnet-class The tidy edge list, one row per edge.
#' @export
as.data.frame.psychnet <- function(x, row.names = NULL, optional = FALSE, ...,
                                   include_zero = FALSE) {
  # The stored $edges uses integer node indices (cograph format); the tidy
  # accessor returns the label-based edge list from the weight matrix.
  .edges_df(x$weights, x$nodes$label, x$directed, include_zero = include_zero)
}

#' @describeIn psychnet-class The printout plus the edge-weight range.
#' @export
summary.psychnet <- function(object, ...) {
  print(object)
  edges <- as.data.frame(object)
  ew <- edges$weight
  if (length(ew)) {
    cat(sprintf("  edge weight: range [%.3f, %.3f], mean %.3f\n",
                min(ew), max(ew), mean(ew)))
  }
  invisible(edges)
}
