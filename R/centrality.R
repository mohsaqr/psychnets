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

  builtin <- c("strength", "expected_influence", "betweenness", "closeness")
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
  data.frame(node = labs, cols, row.names = NULL, stringsAsFactors = FALSE,
             check.names = FALSE)
}
