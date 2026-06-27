# The shared result object returned by every psychnet estimator. It mirrors the
# Nestimate / cograph "netobject" so a fitted psychnet network is plottable by
# cograph::splot() and comparable by Nestimate, while staying lean: it carries
# only the estimated content (the weighted network, the tidy node and edge
# tables, the sample size, the estimator name, and method-specific extras such
# as the precision matrix, the EBIC penalty, and the KKT certificate).

# Tidy one-row-per-edge table from a weighted adjacency matrix.
#' @noRd
.edges_df <- function(weights, labels, directed, include_zero = FALSE) {
  idx <- if (isTRUE(directed)) which(row(weights) != col(weights), arr.ind = TRUE)
         else which(upper.tri(weights), arr.ind = TRUE)
  w <- weights[idx]
  keep <- if (include_zero) rep(TRUE, length(w)) else abs(w) > 1e-12
  data.frame(from = labels[idx[keep, 1L]], to = labels[idx[keep, 2L]],
             weight = w[keep], stringsAsFactors = FALSE, row.names = NULL)
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
#' @return An object of class `c("psychnet", "cograph_network")`.
#' @noRd
.new_psychnet <- function(graph, labels, method, directed, n_obs,
                          extra = list()) {
  if (length(labels) != ncol(graph)) {
    stop(sprintf(paste0("labels length (%d) does not match the network ",
                        "dimension (%d); a non-numeric or zero-variance column ",
                        "may have been dropped from the data."),
                 length(labels), ncol(graph)), call. = FALSE)
  }
  weights <- graph
  diag(weights) <- 0
  dimnames(weights) <- list(labels, labels)
  nodes <- data.frame(id = seq_along(labels), label = labels, name = labels,
                      stringsAsFactors = FALSE)
  edges <- .edges_df(weights, labels, directed)
  structure(
    c(list(weights = weights, nodes = nodes, edges = edges,
           directed = directed, method = method, n = n_obs), extra),
    class = c("psychnet", "cograph_network")
  )
}

#' Back-compatible field access for a psychnet object
#'
#' The canonical (`str`-visible) fields are the lean netobject set; this method
#' adds virtual aliases so older/external accessors keep working without storing
#' redundant fields.
#'
#' @param x A `psychnet` object.
#' @param name Field name. Canonical fields plus the legacy aliases `graph`
#'   (= `weights`), `labels` (= `nodes$label`), `n_nodes`, `n_edges`, and
#'   `n_obs` (= `n`).
#' @return The requested field, or `NULL` if neither a canonical field nor a
#'   known alias.
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

#' Print a psychnet network
#'
#' @param x A `psychnet` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet <- function(x, ...) {
  cat(sprintf("<psychnet> %s network\n", x$method))
  cat(sprintf("  nodes: %d   edges: %d   (%s)\n",
              nrow(x$nodes), nrow(x$edges),
              if (isTRUE(x$directed)) "directed" else "undirected"))
  if (!is.null(x$lambda)) {
    cat(sprintf("  lambda: %.4g   gamma: %.2g\n", x$lambda, x$gamma))
  }
  if (!is.null(x$kkt)) {
    cat(sprintf("  optimality (KKT residual): %.2e\n", x$kkt))
  }
  invisible(x)
}

#' Tidy edge list for a psychnet network
#'
#' @param x A `psychnet` object.
#' @param ... Unused.
#' @param include_zero If TRUE, keep zero-weight (absent) edges. Default FALSE.
#' @param row.names,optional Ignored (for S3 consistency).
#' @return A one-row-per-edge `data.frame` with columns `from`, `to`, `weight`.
#' @export
as.data.frame.psychnet <- function(x, row.names = NULL, optional = FALSE, ...,
                                   include_zero = FALSE) {
  if (include_zero) {
    return(.edges_df(x$weights, x$nodes$label, x$directed, include_zero = TRUE))
  }
  x$edges
}

#' Summarize a psychnet network
#'
#' @param object A `psychnet` object.
#' @param ... Unused.
#' @return The tidy edge list (invisibly); prints a summary as a side effect.
#' @export
summary.psychnet <- function(object, ...) {
  print(object)
  ew <- object$edges$weight
  if (length(ew)) {
    cat(sprintf("  edge weight: range [%.3f, %.3f], mean %.3f\n",
                min(ew), max(ew), mean(ew)))
  }
  invisible(object$edges)
}
