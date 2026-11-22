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

# Planarity of an undirected graph by the left-right planarity test
# (de Fraysseix & Rosenstiehl 1985; Brandes 2009), the linear-time test that
# NetworkX's check_planarity() implements; it decides the same property as the
# Boyer & Myrvold (2004) test. `adj` is a symmetric logical adjacency matrix.
# Both depth-first phases are written iteratively (an explicit vertex stack),
# so the depth of the search is not bounded by R's expression nesting. The
# phases are sequential by nature (each step reads the state the previous one
# left), so they are `while` loops, and the roots are visited in turn.
#' @noRd
.is_planar <- function(adj) {
  p <- nrow(adj)
  m <- sum(adj[upper.tri(adj)])
  # Every graph with fewer than five nodes or fewer than nine edges (K3,3 has
  # nine) is planar; a planar graph has at most 3p - 6 edges (Euler).
  if (p < 5L || m < 9L) return(TRUE)
  if (m > 3L * p - 6L) return(FALSE)
  nbrs <- lapply(seq_len(p), function(v) which(adj[v, ]))
  st <- new.env(parent = emptyenv())
  st$from <- st$to <- st$lowpt <- st$lowpt2 <- st$nesting <- integer(m)
  st$ref <- st$lowpt_edge <- st$stack_bottom <- integer(m)
  st$ne <- 0L
  st$oriented <- matrix(FALSE, p, p)
  st$height <- rep(NA_integer_, p)
  st$parent_edge <- integer(p)

  # --- phase 1: orientation, heights, lowpoints and nesting depths ---------
  finish_edge <- function(ep, v) {
    st$nesting[ep] <- 2L * st$lowpt[ep] + as.integer(st$lowpt2[ep] < st$height[v])
    e <- st$parent_edge[v]
    if (e > 0L) {
      if (st$lowpt[ep] < st$lowpt[e]) {
        st$lowpt2[e] <- min(st$lowpt[e], st$lowpt2[ep])
        st$lowpt[e] <- st$lowpt[ep]
      } else if (st$lowpt[ep] > st$lowpt[e]) {
        st$lowpt2[e] <- min(st$lowpt2[e], st$lowpt[ep])
      } else {
        st$lowpt2[e] <- min(st$lowpt2[e], st$lowpt2[ep])
      }
    }
  }
  roots <- integer()
  next_nb <- rep(1L, p)
  while (anyNA(st$height)) {
    root <- which(is.na(st$height))[1L]
    st$height[root] <- 0L
    roots <- c(roots, root)
    stack <- root
    while (length(stack)) {
      v <- stack[length(stack)]
      if (next_nb[v] > length(nbrs[[v]])) {
        stack <- stack[-length(stack)]
        e <- st$parent_edge[v]
        if (e > 0L) finish_edge(e, st$from[e])
        next
      }
      w <- nbrs[[v]][next_nb[v]]
      next_nb[v] <- next_nb[v] + 1L
      if (st$oriented[v, w]) next
      st$oriented[v, w] <- st$oriented[w, v] <- TRUE
      st$ne <- st$ne + 1L
      ep <- st$ne
      st$from[ep] <- v; st$to[ep] <- w
      st$lowpt[ep] <- st$lowpt2[ep] <- st$height[v]
      if (is.na(st$height[w])) {
        st$parent_edge[w] <- ep
        st$height[w] <- st$height[v] + 1L
        stack <- c(stack, w)
      } else {
        st$lowpt[ep] <- st$height[w]
        finish_edge(ep, v)
      }
    }
  }

  # --- phase 2: testing, with a stack of conflict pairs --------------------
  # A conflict pair holds a left and a right interval of return edges, each
  # (low, high) with 0 for "none"; pairs are numbered, the stack holds numbers.
  out_edges <- split(seq_len(m), factor(st$from, levels = seq_len(p)))
  ordered <- lapply(out_edges, function(ids) ids[order(st$nesting[ids], ids)])
  st$Ll <- st$Lh <- st$Rl <- st$Rh <- integer()
  st$np <- 0L
  st$S <- integer()
  top_S <- function() if (length(st$S)) st$S[[length(st$S)]] else 0L
  pop_S <- function() {
    k <- st$S[[length(st$S)]]
    st$S <- st$S[-length(st$S)]
    k
  }
  push_S <- function(k) st$S <- c(st$S, k)
  new_pair <- function(rl = 0L, rh = 0L) {
    st$np <- st$np + 1L
    k <- st$np
    st$Ll[k] <- 0L; st$Lh[k] <- 0L; st$Rl[k] <- rl; st$Rh[k] <- rh
    k
  }
  swap <- function(k) {
    ll <- st$Ll[k]; lh <- st$Lh[k]
    st$Ll[k] <- st$Rl[k]; st$Lh[k] <- st$Rh[k]
    st$Rl[k] <- ll; st$Rh[k] <- lh
  }
  l_empty <- function(k) st$Ll[k] == 0L && st$Lh[k] == 0L
  r_empty <- function(k) st$Rl[k] == 0L && st$Rh[k] == 0L
  l_conflicting <- function(k, b) !l_empty(k) && st$lowpt[st$Lh[k]] > st$lowpt[b]
  r_conflicting <- function(k, b) !r_empty(k) && st$lowpt[st$Rh[k]] > st$lowpt[b]
  lowest <- function(k) {
    if (l_empty(k)) return(st$lowpt[st$Rl[k]])
    if (r_empty(k)) return(st$lowpt[st$Ll[k]])
    min(st$lowpt[st$Ll[k]], st$lowpt[st$Rl[k]])
  }
  set_ref <- function(edge, value) if (edge > 0L) st$ref[edge] <- value

  add_constraints <- function(ei, e) {
    P <- new_pair()
    # merge the return edges of ei into the right interval of P
    repeat {
      Q <- pop_S()
      if (!l_empty(Q)) swap(Q)
      if (!l_empty(Q)) return(FALSE)
      if (st$lowpt[st$Rl[Q]] > st$lowpt[e]) {
        if (r_empty(P)) st$Rh[P] <- st$Rh[Q] else set_ref(st$Rl[P], st$Rh[Q])
        st$Rl[P] <- st$Rl[Q]
      } else {
        set_ref(st$Rl[Q], st$lowpt_edge[e])
      }
      if (top_S() == st$stack_bottom[ei]) break
    }
    # merge the conflicting return edges of the earlier siblings into the left
    repeat {
      t <- top_S()
      if (t == 0L || !(l_conflicting(t, ei) || r_conflicting(t, ei))) break
      Q <- pop_S()
      if (r_conflicting(Q, ei)) swap(Q)
      if (r_conflicting(Q, ei)) return(FALSE)
      set_ref(st$Rl[P], st$Rh[Q])
      if (st$Rl[Q] != 0L) st$Rl[P] <- st$Rl[Q]
      if (l_empty(P)) st$Lh[P] <- st$Lh[Q] else set_ref(st$Ll[P], st$Lh[Q])
      st$Ll[P] <- st$Ll[Q]
    }
    if (!(l_empty(P) && r_empty(P))) push_S(P)
    TRUE
  }

  remove_back_edges <- function(e) {
    u <- st$from[e]
    while (length(st$S) && lowest(top_S()) == st$height[u]) pop_S()
    if (length(st$S)) {
      P <- pop_S()
      while (st$Lh[P] != 0L && st$to[st$Lh[P]] == u) st$Lh[P] <- st$ref[st$Lh[P]]
      if (st$Lh[P] == 0L && st$Ll[P] != 0L) {
        set_ref(st$Ll[P], st$Rl[P])
        st$Ll[P] <- 0L
      }
      while (st$Rh[P] != 0L && st$to[st$Rh[P]] == u) st$Rh[P] <- st$ref[st$Rh[P]]
      if (st$Rh[P] == 0L && st$Rl[P] != 0L) {
        set_ref(st$Rl[P], st$Ll[P])
        st$Rl[P] <- 0L
      }
      push_S(P)
    }
    # the reference of e is a highest return edge
    if (st$lowpt[e] < st$height[u] && length(st$S)) {
      t <- top_S()
      hl <- st$Lh[t]; hr <- st$Rh[t]
      st$ref[e] <- if (hl != 0L && (hr == 0L || st$lowpt[hl] > st$lowpt[hr])) hl else hr
    }
  }

  next_edge <- rep(1L, p)
  started <- logical(m)
  planar <- TRUE
  for (root in roots) {
    stack <- root
    while (length(stack) && planar) {
      v <- stack[length(stack)]
      if (next_edge[v] > length(ordered[[v]])) {
        stack <- stack[-length(stack)]
        e <- st$parent_edge[v]
        if (e > 0L) remove_back_edges(e)
        next
      }
      ei <- ordered[[v]][next_edge[v]]
      w <- st$to[ei]
      if (!started[ei]) {
        started[ei] <- TRUE
        st$stack_bottom[ei] <- top_S()
        if (st$parent_edge[w] == ei) {            # a tree edge: descend first
          stack <- c(stack, w)
          next
        }
        st$lowpt_edge[ei] <- ei                   # a back edge
        push_S(new_pair(rl = ei, rh = ei))
      }
      # integrate the return edges of ei
      next_edge[v] <- next_edge[v] + 1L
      if (st$lowpt[ei] < st$height[v]) {
        e <- st$parent_edge[v]
        if (identical(ei, ordered[[v]][1L])) {
          st$lowpt_edge[e] <- st$lowpt_edge[ei]
        } else if (!add_constraints(ei, e)) {
          planar <- FALSE
        }
      }
    }
    if (!planar) break
  }
  planar
}

#' Check that a network has the structure of a TMFG
#'
#' @description
#' `tmfg_certificate()` checks that a network from [tmfg_network()] has the
#' structure of a triangulated maximally filtered graph (TMFG). A TMFG is
#' defined by its structure rather than as the optimum of an objective: for
#' \eqn{p} variables it is a connected planar graph with exactly
#' \eqn{3(p - 2)} edges that is chordal (every cycle of four or more nodes has
#' a chord) and is built from four-variable cliques (Massara et al., 2016).
#'
#' These properties are what make the graph a TMFG and what the LoGo model
#' ([logo_network()]) relies on, so checking them shows that a returned
#' network is a valid TMFG, for instance that no edge was lost or added.
#' `tmfg_certificate()` counts how many of the defining properties the network
#' violates.
#'
#' It returns that count: 0 confirms a valid TMFG.
#'
#' [certificate()] reports the same value for a fitted TMFG network. The
#' check does not assess which correlations the construction chose to retain.
#'
#' @details
#' Four properties are checked on the stored graph: the number of edges is
#' \eqn{3(p - 2)}; the graph is chordal, tested by maximum-cardinality search;
#' it is connected; and it is planar. Planarity is decided for any graph by
#' the left-right planarity test (de Fraysseix & Rosenstiehl, 1985; Brandes,
#' 2009), the linear-time test of NetworkX's `check_planarity()`, which decides
#' the same property as the test of Boyer and Myrvold (2004). Each property is
#' tested on its own, so a graph that breaks several of them, such as
#' \eqn{K_{3,3}} (too few edges, not chordal, not planar), counts each one.
#'
#' @param x A TMFG network from [tmfg_network()] ([psychnet-class]).
#' @return A single integer from 0 to 4, the number of violated properties;
#'   0 confirms a valid TMFG.
#' @section Conditions:
#'   `psychnet_bad_input` when `x` is not a network from [tmfg_network()].
#' @seealso [tmfg_network()]; [certificate()]; [lmg_certificate()], the other
#'   structural check.
#' @references Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network
#'   filtering for big data: Triangulated maximally filtered graph. *Journal of
#'   Complex Networks*, 5(2), 161-178. \doi{10.1093/comnet/cnw015}
#'
#'   de Fraysseix, H., & Rosenstiehl, P. (1985). A characterization of planar
#'   graphs by Tremaux orders. *Combinatorica*, 5(2), 127-135.
#'   \doi{10.1007/BF02579375}
#'
#'   Brandes, U. (2009). *The left-right planarity test*. Manuscript,
#'   University of Konstanz.
#'
#'   Boyer, J. M., & Myrvold, W. J. (2004). On the cutting edge: Simplified
#'   O(n) planarity by edge addition. *Journal of Graph Algorithms and
#'   Applications*, 8(3), 241-273. \doi{10.7155/jgaa.00091}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 6), 200, 6)
#' tmfg_certificate(tmfg_network(x))
#' @export
tmfg_certificate <- function(x) {
  if (!inherits(x, "psychnet") || is.null(x$adjacency)) {
    .psn_stop("`x` must be a TMFG network from tmfg_network().", call. = FALSE)
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
  planar <- .is_planar(adj)
  (!edges_ok) + (!chordal) + (!connected) + (!planar)
}

#' Filter a correlation network with the triangulated maximally filtered graph
#'
#' @description
#' `tmfg_network()` filters a correlation matrix to a sparse backbone with the
#' triangulated maximally filtered graph (TMFG). The TMFG retains, for
#' \eqn{p} variables, exactly \eqn{3(p - 2)} edges, chosen greedily to carry
#' as much of the absolute correlation as possible while the graph stays
#' planar and built from tetrahedra (four-variable cliques) (Massara et al.,
#' 2016).
#'
#' Like other information-filtering networks, the TMFG keeps the strongest
#' associations while imposing a fixed, interpretable structure. It shows
#' which associations form the core of the correlation structure, without a
#' sample-size-dependent selection rule. `tmfg_network()` builds the TMFG from
#' raw data or a correlation matrix.
#'
#' It returns the graph with the signed correlation of each retained pair as
#' the edge weight; the edges are marginal correlations, not partial
#' correlations. The structure that defines the graph is checked by
#' [certificate()].
#'
#' For a Gaussian graphical model with partial correlations on the TMFG graph
#' use [logo_network()]; for an unfiltered correlation network use
#' [cor_network()]. The graph is the one built by `NetworkToolbox::TMFG()`.
#'
#' @details
#' The construction starts from the four variables with the largest sums of
#' above-average absolute correlations, joined as a tetrahedron. It then adds
#' the remaining variables one at a time, each time choosing the variable and
#' the triangular face of the current graph for which the variable's summed
#' absolute correlation with the face's three variables is largest, and
#' joining the variable to those three. The result is a chordal planar graph
#' with \eqn{3(p - 2)} edges, which [tmfg_certificate()] verifies.
#'
#' The graph always has \eqn{3(p - 2)} edges, but a retained pair whose
#' correlation is exactly zero has edge weight 0 and is therefore absent from
#' `as.data.frame()` and from the printed edge count. With continuous data this
#' happens only for a correlation matrix with exact zeros, such as a
#' constructed or block-diagonal `cor_matrix`.
#'
#' @param data A numeric data frame or matrix with one row per observation.
#'   Not needed when `cor_matrix` is given.
#' @param cor_matrix A correlation matrix to use instead of `data`; a
#'   covariance matrix is rescaled to a correlation matrix.
#' @param n Accepted and ignored: the TMFG needs no sample size. It exists so
#'   that a call written for the other estimators, `(cor_matrix = , n = )`,
#'   also works here.
#' @param cor_method Correlation computed from `data`: `"pearson"` (default),
#'   `"spearman"`, `"kendall"`, or `"auto"` (polychoric and polyserial
#'   correlations for ordinal items; see [cor_auto()]).
#' @param na_method Handling of missing values in `data`: `"pairwise"`
#'   (default) or `"listwise"`. See [ebic_glasso()].
#' @param labels Optional node labels, one unique label per node of the
#'   network.
#' @return A fitted network of class `psychnet` ([psychnet-class]) whose edge
#'   weights are the signed correlations on the retained edges.
#'   `as.data.frame()` returns one row per edge with columns `from`, `to` and
#'   `weight`; [certificate()] reports the structural check. The object also
#'   stores the TMFG graph as a logical matrix, including any retained edge of
#'   weight 0 (`adjacency`), its four-variable cliques and three-variable
#'   separators (`cliques`, `separators`; used by [logo_network()]) and the
#'   correlation matrix (`cor_matrix`). The sample size is recorded when the
#'   network is built from `data` and is `NA` with `cor_matrix`.
#' @section Conditions:
#'   `psychnet_bad_input` when there are fewer than four variables,
#'   `cor_matrix` is not a finite, symmetric, positive semi-definite matrix
#'   with a positive diagonal, an argument is invalid, or `labels` do not
#'   match the nodes; `psychnet_constant_variable` (warning) when a constant
#'   or all-missing column is dropped.
#' @seealso [tmfg_certificate()] and [certificate()] for the check;
#'   [logo_network()] for the Gaussian graphical model on this graph;
#'   [psychnet()], which routes `method = "tmfg"` here; [psychnet-class] for
#'   the result.
#' @references Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network
#'   filtering for big data: Triangulated maximally filtered graph. *Journal of
#'   Complex Networks*, 5(2), 161-178. \doi{10.1093/comnet/cnw015}
#' @examples
#' set.seed(1)
#' x <- matrix(stats::rnorm(200 * 6), 200, 6)
#' tmfg_network(x)
#' @export
tmfg_network <- function(data = NULL, cor_matrix = NULL, n = NULL,
                         cor_method = c("pearson", "spearman", "kendall", "auto"),
                         na_method = c("pairwise", "listwise"), labels = NULL) {
  cor_method <- .psn_match_arg(cor_method)
  na_method <- .psn_match_arg(na_method)
  inp <- .network_input(data, cor_matrix, n, cor_method, na_method, labels,
                        n_role = "ignored")
  S <- inp$S; labels <- inp$labels
  n_obs <- if (is.null(cor_matrix)) inp$n else NA_integer_
  p <- ncol(S)
  if (p < 4L) .psn_stop("TMFG requires at least 4 variables.", call. = FALSE)

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
