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
    # idx x idx contains both s->t and t->s for an undirected graph. Count each
    # unordered source-target pair once before applying the usual normalization.
    sum((outer(sg[idx, v], sg[v, idx], "*") / sg[idx, idx])[on_path],
        na.rm = TRUE) / 2
  }, numeric(1))
  norm <- (n - 1) * (n - 2) / 2
  if (norm > 0) btw <- btw / norm
  stats::setNames(btw, rownames(W))
}

# Edge betweenness: for each present edge a->b, the number of shortest paths
# (over all ordered source/target pairs, each pair's sigma(s,t) geodesics
# weighted 1/sigma(s,t)) that traverse it. Geodesics use inverse
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
#' For each edge, the number of weighted shortest paths (across all node
#' pairs) that pass through it - a high value marks an edge that bridges
#' otherwise distant parts of the network. When a pair is joined by several
#' equally short paths each counts as its share (`1 / number of geodesics`), so
#' the value is the raw, unnormalized edge betweenness (igraph's
#' `edge_betweenness()`), not a proportion: on the path `a-b-c-d-e` the four
#' edges score 4, 6, 6, 4. Geodesics are computed on inverse absolute
#' weights, so strong edges count as short, matching [net_centralities()]'s node
#' betweenness/closeness.
#'
#' @param x A `psychnet` object or a square weighted adjacency matrix. A
#'   `psychnet` is directed when its `directed` field says so; a bare matrix is
#'   treated as directed when it is not symmetric (`isSymmetric()`), with
#'   `x[i, j]` the edge `i -> j`.
#' @param invert If `TRUE` (default) edge weights are inverted to distances
#'   (strong association = short path). Set `FALSE` to treat weights as distances.
#' @param labels Optional node labels (used when `x` is a bare matrix): one
#'   unique label per node.
#' @return A tidy `data.frame`, one row per edge: `from`, `to`,
#'   `edge_betweenness` (count of shortest paths through the edge). Undirected
#'   networks give one row per unordered edge (counting each unordered node
#'   pair once); directed networks one row per arc, counting ordered pairs.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' net_edge_betweenness(ebic_glasso(cor_matrix = S, n = 400))
#' @export
net_edge_betweenness <- function(x, invert = TRUE, labels = NULL) {
  if (inherits(x, "psychnet")) {
    g <- x$weights; labs <- x$nodes$label; directed <- isTRUE(x$directed)
  } else {
    if (!is.matrix(x) && !is.data.frame(x))
      .psn_stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g))
      .psn_stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    .psn_check_labels(labels, ncol(g))
    labs <- if (!is.null(labels)) as.character(labels) else colnames(g)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(g)))
    if (any(!is.finite(g)))
      .psn_stop("`x` must contain only finite edge weights.")
    # A bare matrix carries no `directed` flag: an asymmetric one is a directed
    # network (reporting only its upper triangle would drop half the arcs).
    directed <- !isSymmetric(unname(g))
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

# Align one `centrality_fn` result to the node labels: a named vector by name
# (exactly the node labels, any order), an unnamed one by position. Anything
# else would silently attach values to the wrong nodes.
#' @noRd
.psn_align_custom <- function(v, m, labs) {
  if (!is.numeric(v) && !is.logical(v))
    .psn_stop(sprintf("`centrality_fn` measure '%s' must be numeric.", m))
  if (length(v) != length(labs))
    .psn_stop(sprintf(
      "`centrality_fn` measure '%s' has %d values but the network has %d nodes.",
      m, length(v), length(labs)))
  nm <- names(v)
  if (is.null(nm)) return(stats::setNames(as.numeric(v), labs))
  missing_n <- setdiff(labs, nm)
  extra_n   <- setdiff(nm, labs)
  if (length(missing_n) || length(extra_n) || anyDuplicated(nm))
    .psn_stop(sprintf(paste0(
      "`centrality_fn` measure '%s' is named but its names are not the node ",
      "labels (missing: %s; unknown: %s%s)."), m,
      if (length(missing_n)) paste(missing_n, collapse = ", ") else "none",
      if (length(extra_n)) paste(extra_n, collapse = ", ") else "none",
      if (anyDuplicated(nm)) "; duplicated names" else ""))
  stats::setNames(as.numeric(v[labs]), labs)
}

#' Node centrality
#'
#' Summarizes each node's position in a fitted network as a tidy table. For
#' psychometric networks the recommended measures are strength (the sum of
#' absolute edge weights) and expected influence (the signed sum, so negative
#' edges offset positive ones; Robinaugh, Millner & McNally 2016). Betweenness
#' and closeness follow Opsahl, Agneessens and Skvoretz (2010) on the
#' inverted absolute weights. For a directed network (e.g. [relimp_network()])
#' strength and expected influence sum the outgoing edges.
#'
#' **Scaling of betweenness and closeness.** Both are *normalized*, as igraph's
#' `betweenness(normalized = TRUE)` / `closeness(normalized = TRUE)`:
#' betweenness is divided by the number of node pairs not involving the node,
#' `(p - 1)(p - 2) / 2`, and closeness is the number of reachable nodes over
#' the sum of distances to them (the inverse mean distance). They therefore
#' differ from `qgraph::centrality_auto()`, which reports raw betweenness
#' (path counts) and closeness as `1 / sum(distances)`.
#'
#' @param x A [psychnet] object or a weighted adjacency matrix. For `plot()`,
#'   the `psychnet_centrality` table returned by `net_centralities()`.
#' @param measures Character vector of measures to return (for `plot()`: which
#'   of the returned measure columns to draw; default all). Any of `"strength"`,
#'   `"expected_influence"` (the defaults, recommended for psychometric
#'   networks), `"expected_influence_2step"` (a node's one-step expected
#'   influence plus the one-step expected influence it transmits through its
#'   neighbours), `"betweenness"`, `"closeness"`, plus any names supplied via
#'   `centrality_fn`. Betweenness and closeness are computed on the absolute,
#'   inverted-weight graph and are not generally meaningful on signed networks --
#'   request them only when a downstream comparison needs them.
#' @param centrality_fn Optional function taking the weighted adjacency matrix
#'   (with node-label dimnames) and returning a named list of node-centrality
#'   vectors, used to supply any `measures` not built in. Each vector must have
#'   one numeric value per node. A *named* vector is aligned to the nodes by
#'   label (its names must be exactly the node labels, in any order); an
#'   unnamed vector is taken in node order. A wrong length, missing or extra
#'   names, or a non-numeric value is an error of class `psychnet_bad_input`.
#' @param type For `plot()`: `"bar"` (default; one sorted lollipop panel per
#'   measure) or `"line"` (the qgraph/bootnet centrality plot: one panel per
#'   measure, nodes on a shared vertical axis).
#' @param scale For `plot(type = "line")`: the per-measure transform, `"raw"`
#'   (default; each panel keeps its own axis), `"z"` (z-score per measure) or
#'   `"relative"` (min-max scaled to \[0, 1\]). `type = "bar"` always shows
#'   raw values.
#' @param ... Unused.
#' @return A tidy `data.frame` of class `c("psychnet_centrality",
#'   "data.frame")`, one row per node, with a `node` column and one column per
#'   requested measure (`strength` = sum of absolute edge weights,
#'   `expected_influence` = sum of signed edge weights, by default;
#'   `betweenness` and `closeness` are normalized, see Details); `plot()`
#'   draws it and returns it invisibly. For a [psychnet_group] the result is a
#'   per-group list (see [psychnet_result_group]).
#' @references Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016).
#'   Identifying highly influential nodes in the complicated grief network.
#'   *Journal of Abnormal Psychology*, 125(6), 747-757.
#'   \doi{10.1037/abn0000181}
#'
#'   Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in
#'   weighted networks: Generalizing degree and shortest paths. *Social
#'   Networks*, 32(3), 245-251. \doi{10.1016/j.socnet.2010.03.006}
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' cent <- net_centralities(ebic_glasso(cor_matrix = S, n = 250))
#' cent
#' plot(cent)
#' plot(cent, type = "line")
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
      .psn_stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    }
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g)) {
      .psn_stop("`x` must be a square numeric weighted adjacency matrix.",
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
    .psn_stop("`centrality_fn` is required for measures: ",
         paste(external, collapse = ", "), ".", call. = FALSE)
  }
  custom <- if (length(external)) {
    if (!is.function(centrality_fn))
      .psn_stop("`centrality_fn` must be a function.")
    cv <- centrality_fn(g)
    if (!is.list(cv) || is.null(names(cv)))
      .psn_stop("`centrality_fn` must return a named list.", call. = FALSE)
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
        if (is.null(v)) .psn_stop("`centrality_fn` did not return '", m, "'.",
                             call. = FALSE)
        .psn_align_custom(v, m, labs)
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
