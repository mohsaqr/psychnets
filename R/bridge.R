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

# Resolve a community specification to a character vector aligned to `labs`
# (NA where a label was given no community). Accepted forms: an unnamed vector
# with one entry per label, in label order; a vector named by label; a list of
# label vectors, one element per community. A label claimed by two communities
# and a name/label that is not in `labs` are errors (psychnet_bad_input) -- a
# silent first-match or a silently ignored label would misassign nodes.
# `what` ("node" / "item") words the messages.
#' @noRd
.psn_resolve_communities <- function(communities, labs, what = "node") {
  if (missing(communities) || is.null(communities) || !length(communities))
    .psn_stop("`communities` is required: one community per ", what, ".")
  labs <- as.character(labs)
  if (is.list(communities)) {
    ok <- vapply(communities, function(v) is.atomic(v) && !is.null(v),
                 logical(1))
    if (!all(ok))
      .psn_stop("A list `communities` must hold one vector of ", what,
                " labels per community.")
    cnames <- names(communities) %||% as.character(seq_along(communities))
    blank <- is.na(cnames) | !nzchar(cnames)
    cnames[blank] <- as.character(seq_along(communities))[blank]
    members <- as.character(unlist(communities, use.names = FALSE))
    owner   <- rep(cnames, lengths(communities))
    keyed   <- stats::setNames(owner, members)
  } else if (!is.null(names(communities))) {
    keyed <- stats::setNames(as.character(communities), names(communities))
  } else {
    if (length(communities) != length(labs))
      .psn_stop(sprintf(paste0(
        "`communities` must have one entry per %s: it has %d entries but there ",
        "are %d %ss (or pass a vector named by %s label, or a list of label ",
        "vectors)."), what, length(communities), length(labs), what, what))
    return(stats::setNames(as.character(communities), labs))
  }
  dup <- unique(names(keyed)[duplicated(names(keyed))])
  if (length(dup))
    .psn_stop("`communities` assigns ", what, "(s) to more than one community: ",
              paste(dup, collapse = ", "), ".")
  unknown <- setdiff(names(keyed), labs)
  if (length(unknown))
    .psn_stop("`communities` names ", what, "(s) that are not in the network: ",
              paste(unknown, collapse = ", "), ".")
  stats::setNames(unname(keyed[labs]), labs)
}

#' Compute bridge centrality between communities of nodes
#'
#' @description
#' `net_bridge()` computes the bridge centrality of each node in a network whose
#' nodes belong to known communities, such as the symptoms of two disorders or
#' the items of different scales. Bridge centrality (Jones, Ma & McNally, 2021)
#' describes how strongly a node is connected to the communities other than its
#' own; nodes with high bridge centrality are those through which the
#' communities are linked, for example symptoms that may connect two disorders
#' in comorbidity.
#'
#' The function computes five measures for each node of an undirected network:
#' bridge strength, the sum of the absolute weights of its edges to other
#' communities; one-step bridge expected influence (`bridge_ei1`), the sum of
#' the signed weights of those edges; two-step bridge expected influence
#' (`bridge_ei2`), which adds the signed weights of the paths of two edges that
#' end in another community; bridge betweenness, the number of shortest paths
#' between nodes of different communities that pass through the node; and
#' bridge closeness, the inverse of the mean distance from the node to the
#' nodes of other communities. Shortest paths use the positive edges, with
#' weights inverted to distances. The communities are given by the user;
#' `normalize` divides the measures by the number of nodes they can involve.
#'
#' It returns one row per node, with its community and the five measures.
#'
#' For centrality without regard to communities, use [net_centralities()]. To
#' merge each community into a single node and estimate the network between
#' communities, use [net_aggregate()].
#'
#' @param x An undirected fitted network ([psychnet-class]) or a symmetric
#'   weighted adjacency matrix. For `plot()`, the table returned by
#'   `net_bridge()`.
#' @param communities The community of each node, in one of three forms: an
#'   unnamed vector with one entry per node, in node order; a vector named by
#'   node label; or a list of node-label vectors, one per community, whose
#'   names are the community names, e.g.
#'   `list(anxiety = c("A1", "A2"), mood = c("M1", "M2"))`. Every node must
#'   belong to exactly one community, and there must be at least two.
#' @param normalize If `TRUE`, divide bridge strength and one-step bridge
#'   expected influence by the number `q` of nodes in other communities,
#'   bridge betweenness by `(p - 1) q`, and two-step bridge expected influence
#'   by `q + (p - 1)(q - q / (p - 1))`, where `p` is the number of nodes, as
#'   `networktools::bridge(normalize = TRUE)` does. Bridge closeness is
#'   unchanged. Default `FALSE`.
#' @param labels Optional node labels for a bare matrix, one unique label per
#'   node.
#' @param measures For `plot()`: the bridge columns to draw. Default: all
#'   five.
#' @param ... Unused.
#' @return A data frame of class `psychnet_bridge`, one row per node, with
#'   columns `node`, `community`, `bridge_strength`, `bridge_betweenness`,
#'   `bridge_closeness`, `bridge_ei1` and `bridge_ei2`. `plot()` draws one
#'   sorted panel per measure, each node's point coloured and shaped by its
#'   community, with a legend of the communities, and returns `x` invisibly.
#' @section Conditions:
#'   `psychnet_bad_input` when `x` is directed, is not a finite symmetric
#'   square matrix, or `labels` is invalid; when `communities` lists a node
#'   twice, names a node that is not in the network, leaves a node without a
#'   community, has the wrong length, or defines fewer than two communities.
#'   `plot()` raises `psychnet_bad_input` for a `measures` entry that is not a
#'   bridge column, `psychnet_plot_skipped` (warning) when it skips a measure
#'   without finite values, and `psychnet_nothing_to_plot` (an error) when no
#'   requested measure has one.
#' @seealso [net_centralities()] for centrality without communities;
#'   [net_aggregate()] for the network between communities.
#' @references Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality.
#'   *Multivariate Behavioral Research*, 56(2), 353-367.
#' @examples
#' S <- 0.3^abs(outer(1:6, 1:6, "-"))
#' fit <- ebic_glasso(cor_matrix = S, n = 400)
#' br <- net_bridge(fit, communities = c(1, 1, 1, 2, 2, 2))
#' br
#' plot(br)
#' @export
net_bridge <- function(x, communities, normalize = FALSE, labels = NULL) {
  if (inherits(x, "psychnet")) {
    if (isTRUE(x$directed))
      .psn_stop("net_bridge() requires an undirected network.", call. = FALSE)
    g <- x$weights; labs <- x$nodes$label
  }
  else {
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
  }
  if (any(!is.finite(g)) || any(abs(g - t(g)) > 1e-8))
    .psn_stop("`x` must be a finite symmetric (undirected) matrix.", call. = FALSE)
  p <- nrow(g)
  diag(g) <- 0
  rownames(g) <- colnames(g) <- labs

  communities <- unname(.psn_resolve_communities(communities, labs, "node"))
  if (anyNA(communities))
    .psn_stop("`communities` gives no community for node(s): ",
              paste(labs[is.na(communities)], collapse = ", "), ".")
  .psn_check_flag(normalize, "normalize")
  if (length(unique(communities)) < 2L)
    .psn_stop("Bridge centrality needs at least two communities.", call. = FALSE)

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
    # networktools::bridge(normalize = TRUE): q = number of other-community
    # nodes; the 2-step divisor uses the number of NODES p (not of communities):
    # q + (p - 1) * (q - q / (p - 1)). p >= 2 is guaranteed by >= 2 communities.
    np <- rowSums(other)                      # other-community node count per node
    bridge_strength <- bridge_strength / np
    ei1 <- ei1 / np
    btw <- btw / ((p - 1) * np)
    ei2 <- ei2 / (np + (p - 1) * (np - np / (p - 1)))
  }

  out <- data.frame(node = labs, community = communities,
                    bridge_strength = bridge_strength,
                    bridge_betweenness = btw, bridge_closeness = clo,
                    bridge_ei1 = ei1, bridge_ei2 = ei2,
                    stringsAsFactors = FALSE, row.names = NULL)
  class(out) <- c("psychnet_bridge", "data.frame")
  out
}

#' @describeIn net_bridge One sorted panel per bridge measure, with points
#'   coloured and shaped by community.
#' @export
plot.psychnet_bridge <- function(x, measures = NULL, ...) {
  avail <- c("bridge_strength", "bridge_betweenness", "bridge_closeness",
             "bridge_ei1", "bridge_ei2")
  cols <- .psn_plot_measures(x, measures %||% avail, avail, "plot.psychnet_bridge")
  comm  <- as.character(x$community)
  ucomm <- unique(comm)
  key   <- .psn_group_key(length(ucomm))
  gi    <- match(comm, ucomm)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  # One extra bottom line holds the community legend shared by every panel.
  graphics::par(mfrow = c(1, length(cols)), oma = c(2.5, 0, 0, 0),
                mar = c(4, .psn_left_margin(x$node), 3, 1))
  invisible(lapply(cols, function(m)
    .psn_pointrange(x[[m]], x[[m]], x[[m]], x$node, main = m, xlab = "",
                    ref = numeric(0), sort = TRUE,
                    point_col = key$col[gi], point_pch = key$pch[gi])))
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
                new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", legend = ucomm, col = key$col, pch = key$pch,
                   horiz = TRUE, bty = "n", cex = 0.8, title = "community",
                   xpd = NA)
  invisible(x)
}
