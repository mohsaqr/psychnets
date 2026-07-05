# Bridge centrality (Jones, Ma & McNally 2021), ported from networktools::bridge
# for the undirected weighted case (the psychometric GGM setting). Bridge metrics
# measure how strongly a node connects to communities OTHER than its own:
#   * bridge_strength      - sum of |edges| to other communities
#   * bridge_betweenness   - times the node lies on a shortest path between two
#                            other-community nodes (positive edges, inverse-weight)
#   * bridge_closeness     - inverse mean distance to other-community nodes
#   * bridge_ei1 / bridge_ei2 - 1- and 2-step expected influence restricted to
#                            other communities (signed)
# Communities must be supplied (no automatic detection); this matches the way
# bridge centrality is reported in applied work.

#' Bridge centrality
#'
#' Computes bridge centrality for an undirected weighted network: how strongly
#' each node connects to communities other than its own. You supply the
#' community membership; psychnets does not detect it.
#'
#' @param x A `psychnet` object or a square weighted adjacency matrix.
#' @param communities Community membership, one entry per node: a vector aligned
#'   to the node order, or a named vector / list keyed by node label.
#' @param normalize If `TRUE`, divide each metric by the number of available
#'   other-community nodes (comparable across differently sized networks).
#'   Default `FALSE`.
#' @param labels Optional node labels (used when `x` is a bare matrix).
#' @return A tidy `data.frame` (class `psychnet_bridge`), one row per node, with
#'   columns `node`, `community`, `bridge_strength`, `bridge_betweenness`,
#'   `bridge_closeness`, `bridge_ei1`, `bridge_ei2`. Visualise with
#'   [plot.psychnet_bridge()].
#' @references Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality.
#'   *Multivariate Behavioral Research*, 56(2), 353-367.
#' @examples
#' S <- 0.3^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 400)
#' net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2))
#' @export
net_bridge <- function(x, communities, normalize = FALSE, labels = NULL) {
  if (inherits(x, "psychnet")) { g <- x$weights; labs <- x$nodes$label }
  else {
    if (!is.matrix(x) && !is.data.frame(x))
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g))
      stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    labs <- if (!is.null(labels)) labels else colnames(g)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(g)))
  }
  p <- nrow(g)
  diag(g) <- 0
  rownames(g) <- colnames(g) <- labs

  # Resolve communities to a character vector aligned to `labs`.
  if (is.list(communities)) {
    stacked <- utils::stack(communities)
    comm <- stats::setNames(as.character(stacked$ind), as.character(stacked$values))
    communities <- comm[labs]
  } else if (!is.null(names(communities))) {
    communities <- communities[labs]
  }
  communities <- as.character(communities)
  if (length(communities) != p)
    stop("`communities` must have one entry per node.", call. = FALSE)
  if (anyNA(communities))
    stop("`communities` could not be aligned to every node label.", call. = FALSE)
  if (length(unique(communities)) < 2L)
    stop("Bridge centrality needs at least two communities.", call. = FALSE)

  other <- outer(communities, communities, "!=")   # TRUE where comms differ

  # 1-step quantities: cross-community absolute strength and signed EI.
  bridge_strength <- rowSums(abs(g) * other)
  ei1 <- rowSums(g * other)

  # Shortest paths on positive edges with inverse weights (networktools' g2).
  sp <- .psn_floyd_warshall(g, invert = TRUE)
  D <- sp$D; sg <- sp$sigma

  # Bridge betweenness: # shortest cross-community paths each node lies inside.
  btw <- vapply(seq_len(p), function(i) {
    s <- seq_len(p)[-i]
    d_si <- D[s, i]; d_it <- D[i, s]
    d_st <- D[s, s]
    through <- outer(d_si, d_it, "+")
    on_path <- is.finite(d_st) & d_st > 0 & abs(through - d_st) < 1e-10
    cnt <- outer(sg[s, i], sg[i, s], "*")
    cross <- other[s, s]
    sum(cnt[on_path & cross])
  }, numeric(1)) / 2          # undirected: each unordered pair counted twice

  # Bridge closeness: inverse mean distance to other-community nodes.
  clo <- vapply(seq_len(p), function(i) {
    d <- D[i, other[i, ]]
    d <- d[is.finite(d) & d > 0]
    if (length(d) == 0L) 0 else 1 / mean(d)
  }, numeric(1))

  # Bridge 2-step EI: ei1 plus influence routed through every other community.
  ucomm <- unique(communities)
  infcomm <- lapply(ucomm, function(cc) as.vector(g %*% (communities == cc)))
  names(infcomm) <- ucomm
  ei2 <- vapply(seq_len(p), function(i) {
    others_c <- setdiff(ucomm, communities[i])
    add <- sum(vapply(others_c, function(cc) sum(g[i, ] * infcomm[[cc]]),
                      numeric(1)))
    ei1[i] + add
  }, numeric(1))

  if (isTRUE(normalize)) {
    np <- rowSums(other)                      # other-community node count per node
    k  <- length(ucomm)
    bridge_strength <- bridge_strength / np
    ei1 <- ei1 / np
    btw <- btw / ((length(communities) - 1) * np)
    ei2 <- ei2 / (np + (k - 1) * (np - np / (k - 1)))
  }

  out <- data.frame(node = labs, community = communities,
                    bridge_strength = bridge_strength,
                    bridge_betweenness = btw, bridge_closeness = clo,
                    bridge_ei1 = ei1, bridge_ei2 = ei2,
                    stringsAsFactors = FALSE, row.names = NULL)
  class(out) <- c("psychnet_bridge", "data.frame")
  out
}

#' Plot bridge centrality
#'
#' One sorted horizontal panel per bridge measure (nodes coloured by community).
#'
#' @param x A `psychnet_bridge` data frame from [net_bridge()].
#' @param measures Which bridge columns to draw. Default: all five.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the plot it draws.
#' @examples
#' S <- 0.3^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 400)
#' plot(net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2)))
#' @export
plot.psychnet_bridge <- function(x, measures = NULL, ...) {
  cols <- if (is.null(measures))
    c("bridge_strength", "bridge_betweenness", "bridge_closeness",
      "bridge_ei1", "bridge_ei2") else measures
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mfrow = c(1, length(cols)),
                mar = c(4, .psn_left_margin(x$node), 3, 1))
  for (m in cols)
    .psn_pointrange(x[[m]], x[[m]], x[[m]], x$node, main = m, xlab = "",
                    ref = numeric(0), sort = TRUE)
  invisible(x)
}
