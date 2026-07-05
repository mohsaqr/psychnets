# Node centrality for psychometric networks (pure base R). The default measures
# are strength and expected influence (Robinaugh, Millner & McNally 2016) -- the
# two recommended for psychometric networks because, unlike closeness and
# betweenness, they are well-defined on signed, weighted graphs. Betweenness and
# closeness are available on request (computed on the absolute, inverted weight
# graph), and arbitrary measures via `centrality_fn`, so a resampling caller
# (net_boot) can report whatever a downstream analysis needs without re-deriving
# the shortest-path machinery.

# Weighted shortest-path lengths and shortest-path counts by Floyd-Warshall
# (igraph-free; ported from Nestimate). `invert = TRUE` maps an association
# weight to a distance via 1/w, so a stronger edge is a shorter step.
#' @noRd
.psn_floyd_warshall <- function(W, invert = TRUE) {
  n   <- nrow(W)
  pos <- W > 0
  D            <- matrix(Inf, n, n)
  D[pos]       <- if (invert) 1 / W[pos] else W[pos]
  diag(D)      <- 0
  sigma        <- matrix(0L, n, n)
  sigma[pos]   <- 1L
  diag(sigma)  <- 1L
  Reduce(function(s, k) {
    D <- s$D; sigma <- s$sigma
    new_d <- outer(D[, k], D[k, ], "+")
    new_s <- outer(sigma[, k], sigma[k, ], "*")
    shorter      <- new_d < D & is.finite(new_d)
    equal        <- (new_d == D) & is.finite(new_d) & new_d > 0
    shorter[k, ] <- FALSE; shorter[, k] <- FALSE
    equal[k, ]   <- FALSE; equal[, k]   <- FALSE
    sigma[shorter] <- new_s[shorter]
    sigma[equal]   <- sigma[equal] + new_s[equal]
    new_D <- D; new_D[shorter] <- new_d[shorter]
    list(D = new_D, sigma = sigma)
  }, seq_len(n), list(D = D, sigma = sigma))
}

# Undirected weighted betweenness on the absolute-weight graph.
#' @noRd
.psn_betweenness <- function(g) {
  W <- abs(g); diag(W) <- 0
  n <- nrow(W)
  if (n < 3L) return(stats::setNames(numeric(n), rownames(W)))
  sp <- .psn_floyd_warshall(W, invert = TRUE)
  D <- sp$D; sg <- sp$sigma
  btw <- vapply(seq_len(n), function(v) {
    idx <- seq_len(n)[-v]
    d_svt   <- outer(D[idx, v], D[v, idx], "+")
    d_st    <- D[idx, idx]
    on_path <- is.finite(d_st) & sg[idx, idx] > 0L & abs(d_svt - d_st) < 1e-10
    diag(on_path) <- FALSE
    sum((outer(sg[idx, v], sg[v, idx], "*") / sg[idx, idx])[on_path], na.rm = TRUE)
  }, numeric(1))
  norm <- (n - 1) * (n - 2) / 2
  if (norm > 0) btw <- btw / norm
  stats::setNames(btw, rownames(W))
}

# Edge betweenness: for each present edge a->b, the fraction of shortest paths
# (over all source/target pairs) that traverse it. Geodesics use inverse
# absolute weights (strong edges = short distance), matching the node measures.
.psn_edge_betweenness <- function(W, invert = TRUE) {
  n <- nrow(W)
  EB <- matrix(0, n, n, dimnames = dimnames(W))
  if (n < 2L) return(EB)
  sp <- .psn_floyd_warshall(W, invert)
  D <- sp$D; sg <- sp$sigma
  pos <- W > 0
  len <- matrix(Inf, n, n); len[pos] <- if (invert) 1 / W[pos] else W[pos]
  edges <- which(pos, arr.ind = TRUE)
  vals <- vapply(seq_len(nrow(edges)), function(e) {
    a <- edges[e, 1L]; b <- edges[e, 2L]
    # A shortest s->t path uses edge a->b iff d(s,a) + len(a,b) + d(b,t) = d(s,t);
    # such paths number sigma(s,a)*sigma(b,t), as a fraction of sigma(s,t).
    through <- outer(D[, a], D[b, ], "+") + len[a, b]
    on_path <- is.finite(D) & sg > 0L & abs(through - D) < 1e-9
    diag(on_path) <- FALSE
    sum((outer(sg[, a], sg[b, ]) / sg)[on_path])
  }, numeric(1))
  EB[edges] <- vals
  EB
}

#' Edge betweenness centrality
#'
#' For each edge, the share of weighted shortest paths (across all node pairs)
#' that pass through it - a high value marks an edge that bridges otherwise
#' distant parts of the network. Geodesics are computed on inverse absolute
#' weights, so strong edges count as short, matching [net_centralities()]'s node
#' betweenness/closeness.
#'
#' @param x A `psychnet` object or a square weighted adjacency matrix.
#' @param invert If `TRUE` (default) edge weights are inverted to distances
#'   (strong association = short path). Set `FALSE` to treat weights as distances.
#' @param labels Optional node labels (used when `x` is a bare matrix).
#' @return A tidy `data.frame`, one row per edge: `from`, `to`,
#'   `edge_betweenness`. Undirected networks give one row per unordered edge.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' net_edge_betweenness(ebic_glasso(cor_matrix = S, n = 400))
#' @export
net_edge_betweenness <- function(x, invert = TRUE, labels = NULL) {
  if (inherits(x, "psychnet")) {
    g <- x$weights; labs <- x$nodes$label; directed <- isTRUE(x$directed)
  } else {
    if (!is.matrix(x) && !is.data.frame(x))
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g))
      stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    labs <- if (!is.null(labels)) labels else colnames(g)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(g)))
    directed <- FALSE
  }
  W <- abs(g); diag(W) <- 0
  rownames(W) <- colnames(W) <- labs
  EB <- .psn_edge_betweenness(W, invert = invert)
  mask <- if (directed) (W > 0) else (W > 0 & upper.tri(W))
  ij <- which(mask, arr.ind = TRUE)
  data.frame(from = labs[ij[, 1L]], to = labs[ij[, 2L]],
             edge_betweenness = EB[ij],
             stringsAsFactors = FALSE, row.names = NULL)
}

# Undirected weighted closeness on the absolute-weight graph: (#reachable) /
# sum(distances). Isolated nodes score 0.
#' @noRd
.psn_closeness <- function(g) {
  W <- abs(g); diag(W) <- 0
  n <- nrow(W)
  D <- .psn_floyd_warshall(W, invert = TRUE)$D
  cl <- vapply(seq_len(n), function(v) {
    d <- D[v, ]; d <- d[is.finite(d) & d > 0]
    if (length(d) == 0L) 0 else length(d) / sum(d)
  }, numeric(1))
  stats::setNames(cl, rownames(W))
}

#' Node centrality
#'
#' @param x A [psychnet] object or a weighted adjacency matrix.
#' @param measures Character vector of measures to return. Any of `"strength"`,
#'   `"expected_influence"` (the defaults, recommended for psychometric
#'   networks), `"betweenness"`, `"closeness"`, plus any names supplied via
#'   `centrality_fn`. Betweenness and closeness are computed on the absolute,
#'   inverted-weight graph and are not generally meaningful on signed networks --
#'   request them only when a downstream comparison needs them.
#' @param centrality_fn Optional function taking the weighted adjacency matrix
#'   and returning a named list of node-centrality vectors, used to supply any
#'   `measures` not built in.
#' @param ... Unused.
#' @return A tidy `data.frame`, one row per node, with a `node` column and one
#'   column per requested measure (`strength` = sum of absolute edge weights,
#'   `expected_influence` = sum of signed edge weights, by default).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' net_centralities(ebic_glasso(cor_matrix = S, n = 250))
#' @export
net_centralities <- function(x, measures = c("strength", "expected_influence"),
                             centrality_fn = NULL, ...) {
  if (inherits(x, "psychnet_group")) {
    .reject_multilevel_group(x, "net_centralities")
    return(.group_obj_apply(x, net_centralities, "psychnet_centrality_group",
                            measures = measures, centrality_fn = centrality_fn,
                            ...))
  }
  if (inherits(x, "psychnet")) {
    g <- x$weights
    labs <- x$nodes$label
  } else {
    if (!is.matrix(x) && !is.data.frame(x)) {
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    }
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g)) {
      stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    }
    labs <- colnames(g)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(g)))
  }
  diag(g) <- 0
  rownames(g) <- colnames(g) <- labs

  builtin <- c("strength", "expected_influence", "expected_influence_2step",
               "betweenness", "closeness")
  external <- setdiff(measures, builtin)
  if (length(external) && is.null(centrality_fn)) {
    stop("`centrality_fn` is required for measures: ",
         paste(external, collapse = ", "), ".", call. = FALSE)
  }
  custom <- if (length(external)) {
    cv <- centrality_fn(g)
    if (!is.list(cv)) stop("`centrality_fn` must return a named list.",
                           call. = FALSE)
    cv
  } else NULL

  one <- function(m) {
    switch(m,
      strength            = rowSums(abs(g)),
      expected_influence  = rowSums(g),
      # 2-step EI (Robinaugh et al. 2016): a node's own 1-step EI plus the
      # 1-step EI it transmits through its neighbours, ei2 = ei1 + g %*% ei1.
      expected_influence_2step = { ei1 <- rowSums(g); ei1 + as.vector(g %*% ei1) },
      betweenness         = .psn_betweenness(g),
      closeness           = .psn_closeness(g),
      {
        v <- custom[[m]]
        if (is.null(v)) stop("`centrality_fn` did not return '", m, "'.",
                             call. = FALSE)
        stats::setNames(as.numeric(v), labs)
      })
  }
  cols <- lapply(measures, one)
  names(cols) <- measures
  out <- data.frame(node = labs, cols, row.names = NULL,
                    stringsAsFactors = FALSE, check.names = FALSE)
  # A lightweight subclass over `data.frame` so `plot()` finds a method while
  # every data-frame operation (printing, `[`, `$`, `as.data.frame`) is unchanged.
  class(out) <- c("psychnet_centrality", "data.frame")
  out
}
