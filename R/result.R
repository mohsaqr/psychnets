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

#' Psychometric network objects
#'
#' Every estimator ([cor_network()], [ebic_glasso()], [ising_fit()],
#' [mgm_fit()], ..., and the front door [psychnet()]) returns an object of class
#' `c("psychnet", "cograph_network")`. Read it with `print()`, `summary()` and
#' `as.data.frame()`, and draw it with `plot()`, which hands it to
#' `cograph::splot()`. Diagnostics computed from it ([net_centralities()],
#' [net_boot()], [net_stability()], ...) have their own result classes, with
#' methods documented on those verbs' pages.
#'
#' The stored fields follow the cograph/Nestimate network object: `weights`
#' (the p x p weighted adjacency matrix with zero diagonal; partial
#' correlations for the Gaussian models, symmetrized regression weights for
#' Ising and mixed models), `nodes` (data frame `id`, `label`, `name`), `edges`
#' (integer node indices `from`, `to` and `weight`, the cograph format; use
#' `as.data.frame()` for labelled edges), `directed`, `method` and `n` (sample
#' size), plus method-specific fields such as `precision`, `lambda`, `gamma`,
#' `cor_matrix`, `kkt`, `thresholds`, `types` and `nodewise`. Read the
#' optimality certificate with [certificate()].
#'
#' @param x,object A `psychnet` object.
#' @param name For `$`: a field name. Besides the stored fields, the legacy
#'   aliases `graph` (= `weights`), `labels` (the node labels), `n_nodes`,
#'   `n_edges` and `n_obs` (= `n`) are recognised.
#' @param include_zero For `as.data.frame()`: if `TRUE`, keep zero-weight
#'   (absent) edges. Default `FALSE`.
#' @param row.names,optional Ignored (S3 consistency).
#' @param ... For `plot()`, passed to [cograph::splot()]; ignored by the other
#'   methods.
#' @return `print()` returns `x` invisibly. `as.data.frame()` returns a
#'   data frame with one row per edge and columns `from`, `to` (node labels)
#'   and `weight`. `summary()` prints the network and its edge-weight range and
#'   returns the same edge table invisibly. `plot()` returns the value of
#'   [cograph::splot()] invisibly. `$` returns the requested field, or `NULL`
#'   when the name is neither a stored field nor a known alias.
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

#' @describeIn psychnet-class The estimator, node and edge counts, the
#'   selected penalty and the certificate: the optimality (KKT) residual for
#'   the regularized / constrained estimators, the structural certificate for
#'   TMFG (number of violated invariants) and relimp (decomposition residual).
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
