# Triangulated Maximally Filtered Graph (Massara, Di Matteo & Aste 2016),
# clean-room base R. A planar information-filtering network: greedily insert
# each vertex into the triangular face that maximizes its connectivity, keeping
# the graph chordal and yielding exactly 3(p-2) edges. The greedy vertex-into-
# face insertion is genuinely sequential (each step depends on the current face
# set), so it is a literal loop. The certificate is structural -- there is no
# convex objective -- and is exposed via [tmfg_certificate()].

# Greedy TMFG construction from a similarity matrix. Returns the adjacency and
# the 4-clique / 3-clique-separator decomposition (consumed by LoGo).
#' @noRd
.tmfg_build <- function(W) {
  p <- ncol(W)
  # Seed tetrahedron: the 4 vertices maximizing the sum of their ABOVE-MEAN
  # absolute associations (Massara et al. 2016), the canonical rule used by
  # NetworkToolbox::TMFG. The threshold/mean is taken over the full matrix
  # (diagonal included, as in the reference) before the diagonal is zeroed.
  seed <- order(rowSums(W * (W > mean(W))), decreasing = TRUE)[1:4]
  diag(W) <- 0
  adj <- matrix(FALSE, p, p)
  for (a in 1:3) for (b in (a + 1):4) {
    adj[seed[a], seed[b]] <- adj[seed[b], seed[a]] <- TRUE
  }
  cliques <- list(seed)
  separators <- list()
  faces <- list(seed[c(1, 2, 3)], seed[c(1, 2, 4)],
                seed[c(1, 3, 4)], seed[c(2, 3, 4)])
  remaining <- setdiff(seq_len(p), seed)

  while (length(remaining) > 0L) {
    best_gain <- -Inf; best_v <- NA_integer_; best_f <- NA_integer_
    for (v in remaining) {
      for (fi in seq_along(faces)) {
        f <- faces[[fi]]
        g <- W[v, f[1]] + W[v, f[2]] + W[v, f[3]]
        if (g > best_gain) { best_gain <- g; best_v <- v; best_f <- fi }
      }
    }
    f <- faces[[best_f]]; v <- best_v
    adj[v, f] <- adj[f, v] <- TRUE
    cliques[[length(cliques) + 1L]] <- c(v, f)
    separators[[length(separators) + 1L]] <- f
    faces[[best_f]] <- NULL
    faces <- c(faces, list(c(v, f[1], f[2]), c(v, f[1], f[3]), c(v, f[2], f[3])))
    remaining <- setdiff(remaining, v)
  }
  list(adj = adj, cliques = cliques, separators = separators)
}

# Chordality test by maximum-cardinality search (Tarjan & Yannakakis 1984).
#' @noRd
.is_chordal <- function(adj) {
  p <- ncol(adj)
  wt <- rep(0, p); visited <- rep(FALSE, p); ord <- integer(p)
  for (i in p:1) {
    cand <- which(!visited)
    v <- cand[which.max(wt[cand])]
    ord[i] <- v; visited[v] <- TRUE
    nb <- which(adj[v, ] & !visited); wt[nb] <- wt[nb] + 1
  }
  pos <- integer(p); pos[ord] <- seq_len(p)
  for (i in seq_len(p)) {
    v <- ord[i]
    later <- which(adj[v, ] & pos > pos[v])
    if (length(later) >= 2L) {
      u <- later[which.min(pos[later])]
      others <- setdiff(later, u)
      if (length(others) && !all(adj[u, others])) return(FALSE)
    }
  }
  TRUE
}

#' Structural certificate for a TMFG network
#'
#' A TMFG has no convex objective, so its correctness is certified structurally:
#' a valid TMFG on `p >= 3` nodes has exactly `3(p - 2)` edges, is connected, and
#' is chordal (every cycle of length >= 4 has a chord). Returns a non-negative
#' score that is `0` for a valid TMFG.
#'
#' @param x A [psychnet] object produced by [tmfg_network()].
#' @return Scalar; `0` certifies a valid TMFG (correct edge count, connected,
#'   chordal), otherwise a positive integer counting the violated invariants.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 6), 200, 6)
#' tmfg_certificate(tmfg_network(x))
#' @export
tmfg_certificate <- function(x) {
  if (!inherits(x, "psychnet") || is.null(x$adjacency)) {
    stop("`x` must be a TMFG network from tmfg_network().", call. = FALSE)
  }
  adj <- x$adjacency
  p <- ncol(adj)
  n_edges <- sum(adj[upper.tri(adj)])
  edges_ok <- if (p >= 3L) n_edges == 3L * (p - 2L) else TRUE
  chordal <- .is_chordal(adj)
  # connectivity by breadth-first reach from node 1
  reached <- 1L; frontier <- 1L
  while (length(frontier) > 0L) {
    nb <- which(apply(adj[frontier, , drop = FALSE], 2L, any))
    new <- setdiff(nb, reached)
    reached <- c(reached, new); frontier <- new
  }
  connected <- length(reached) == p
  (!edges_ok) + (!chordal) + (!connected)
}

#' Triangulated Maximally Filtered Graph (TMFG)
#'
#' Builds a sparse, planar, chordal association network by greedily retaining
#' the `3(p - 2)` most informative edges (Massara et al. 2016). Equivalent in
#' purpose to `NetworkToolbox::TMFG()` / `bootnet`'s `"TMFG"` default, pure base
#' R; correctness is certified structurally by [tmfg_certificate()].
#'
#' @param data Numeric data frame or matrix (rows = observations). Optional if
#'   `cor_matrix` is supplied.
#' @param cor_matrix Optional correlation matrix.
#' @param n Accepted and ignored. TMFG is a structural filter and needs no
#'   sample size; the argument exists only so a uniform `(cor_matrix=, n=)` call
#'   shared with the other estimators does not partial-match `na_method`.
#' @param cor_method Correlation when `data` is supplied: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric/polyserial; see [cor_auto()]).
#' @param na_method Missing-data handling when `data` is supplied: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels.
#' @return A `psychnet` object whose `$weights` is the filtered (signed)
#'   correlation matrix on the retained edges, with `$adjacency`, `$cliques`,
#'   `$separators` (the chordal decomposition used by [logo_network()]), and
#'   `$cor_matrix`.
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 6), 200, 6)
#' tmfg_network(x)
#' @export
tmfg_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- match.arg(cor_method)
  na_method <- match.arg(na_method)
  if (is.null(cor_matrix)) {
    ci <- .cor_input(data, method = cor_method, na_method = na_method)
    S <- ci$S; n_obs <- ci$n
    if (is.null(labels)) labels <- ci$labels
  } else {
    S <- .check_cor_matrix(cor_matrix)
    if (is.null(labels)) {
      labels <- colnames(S)
      if (is.null(labels)) labels <- paste0("V", seq_len(ncol(S)))
    }
    n_obs <- NA_integer_
  }
  p <- ncol(S)
  if (p < 4L) stop("TMFG requires at least 4 variables.", call. = FALSE)

  built <- .tmfg_build(abs(S))
  g <- S * built$adj
  diag(g) <- 0
  dimnames(g) <- dimnames(built$adj) <- dimnames(S) <- list(labels, labels)

  .new_psychnet(
    graph = g, labels = labels, method = "tmfg",
    directed = FALSE, n_obs = n_obs,
    extra = list(adjacency = built$adj, cliques = built$cliques,
                 separators = built$separators, cor_matrix = S)
  )
}
