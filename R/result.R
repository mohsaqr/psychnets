# The shared result object returned by every psychnet estimator.

#' Construct a psychnet network object
#'
#' @param graph p x p weighted adjacency matrix (partial correlations for the
#'   Gaussian models, regression-weight matrix for Ising/mgm).
#' @param labels Character node labels.
#' @param method Estimator name.
#' @param directed Logical; TRUE only for inherently directed estimators.
#' @param n_obs Sample size used.
#' @param extra Named list of method-specific fields (e.g. precision, lambda).
#' @return An object of class `psychnet`.
#' @noRd
.new_psychnet <- function(graph, labels, method, directed, n_obs,
                          extra = list()) {
  if (length(labels) != ncol(graph)) {
    stop(sprintf(paste0("labels length (%d) does not match the network ",
                        "dimension (%d); a non-numeric or zero-variance column ",
                        "may have been dropped from the data."),
                 length(labels), ncol(graph)), call. = FALSE)
  }
  g <- graph
  diag(g) <- 0
  n_edges <- if (directed) {
    sum(abs(g) > 1e-12)
  } else {
    sum(abs(g[upper.tri(g)]) > 1e-12)
  }
  dimnames(g) <- list(labels, labels)
  structure(
    c(list(graph = g, labels = labels, method = method,
           directed = directed, n_nodes = length(labels),
           n_edges = n_edges, n_obs = n_obs),
      extra),
    class = "psychnet"
  )
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
              x$n_nodes, x$n_edges,
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
  g <- x$graph
  labs <- x$labels
  if (isTRUE(x$directed)) {
    idx <- which(row(g) != col(g), arr.ind = TRUE)
  } else {
    idx <- which(upper.tri(g), arr.ind = TRUE)
  }
  w <- g[idx]
  keep <- if (include_zero) rep(TRUE, length(w)) else abs(w) > 1e-12
  data.frame(
    from   = labs[idx[keep, 1L]],
    to     = labs[idx[keep, 2L]],
    weight = w[keep],
    stringsAsFactors = FALSE
  )
}

#' Summarize a psychnet network
#'
#' @param object A `psychnet` object.
#' @param ... Unused.
#' @return The tidy edge list (invisibly); prints a summary as a side effect.
#' @export
summary.psychnet <- function(object, ...) {
  print(object)
  ew <- as.data.frame(object)$weight
  if (length(ew)) {
    cat(sprintf("  edge weight: range [%.3f, %.3f], mean %.3f\n",
                min(ew), max(ew), mean(ew)))
  }
  invisible(as.data.frame(object))
}
