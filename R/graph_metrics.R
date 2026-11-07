# Graph-level metrics ported from qgraph: weighted clustering coefficients
# (Watts-Strogatz, Zhang-Horvath, Onnela, Barrat) and the small-world index.
# Pure base-R matrix algebra + a degree-preserving edge-swap null for the
# small-world reference graphs.

# Diagonal of a triple matrix product, A %*% A %*% A, without forming it twice.
.psn_diag_cube <- function(A) rowSums((A %*% A) * t(A))

#' Weighted clustering coefficients
#'
#' Per-node clustering coefficients for a weighted (optionally signed) network:
#' Watts-Strogatz, Zhang-Horvath, Onnela, and Barrat. For signed networks the
#' signed variants are also returned (a closed triangle of three negative or one
#' negative edge lowers the signed coefficient). Ported from
#' `qgraph::clustcoef_auto`.
#'
#' @param x A `psychnet` object or a square weighted adjacency matrix.
#' @param labels Optional node labels (used when `x` is a bare matrix).
#' @return A tidy `data.frame` (class `psychnet_clustering`), one row per node:
#'   `node`, `clustWS`, `signed_clustWS`, `clustZhang`, `signed_clustZhang`,
#'   `clustOnnela`, `signed_clustOnnela`, `clustBarrat`.
#' @references Costantini, G., & Perugini, M. (2014). Generalization of
#'   clustering coefficients to signed correlation networks. *PLoS ONE*, 9(2),
#'   e88669. \doi{10.1371/journal.pone.0088669}
#'
#'   Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#'   'small-world' networks. *Nature*, 393, 440-442. \doi{10.1038/30918}
#'
#'   Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
#'   co-expression network analysis. *Statistical Applications in Genetics and
#'   Molecular Biology*, 4(1), Article 17. \doi{10.2202/1544-6115.1128}
#'
#'   Onnela, J.-P., Saramaki, J., Kertesz, J., & Kaski, K. (2005). Intensity and
#'   coherence of motifs in weighted complex networks. *Physical Review E*,
#'   71(6), 065103. \doi{10.1103/PhysRevE.71.065103}
#'
#'   Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004).
#'   The architecture of complex weighted networks. *Proceedings of the National
#'   Academy of Sciences*, 101(11), 3747-3752. \doi{10.1073/pnas.0400087101}
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' net_clustering(ebic_glasso(cor_matrix = S, n = 400))
#' @export
net_clustering <- function(x, labels = NULL) {
  if (inherits(x, "psychnet")) {
    if (isTRUE(x$directed))
      .psn_stop("net_clustering() requires an undirected network.", call. = FALSE)
    W <- x$weights; labs <- x$nodes$label
  }
  else {
    if (!is.matrix(x) && !is.data.frame(x))
      .psn_stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    W <- as.matrix(x)
    if (!is.numeric(W) || nrow(W) != ncol(W))
      .psn_stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    .psn_check_labels(labels, ncol(W))
    labs <- if (!is.null(labels)) as.character(labels) else colnames(W)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(W)))
  }
  if (any(!is.finite(W)) || any(abs(W - t(W)) > 1e-8))
    .psn_stop("`x` must be a finite symmetric (undirected) matrix.", call. = FALSE)
  diag(W) <- 0
  if (min(W, na.rm = TRUE) < -1 || max(W, na.rm = TRUE) > 1) W <- W / max(abs(W))
  aW <- abs(W)

  # Watts-Strogatz on the sign matrix.
  A <- sign(W)
  aA <- abs(A)
  ki <- colSums(aA); den_ws <- ki * (ki - 1)
  clustWS        <- .psn_diag_cube(aA) / den_ws
  signed_clustWS <- .psn_diag_cube(A)  / den_ws

  # Zhang-Horvath on the weights.
  den_z <- colSums(aW)^2 - colSums(W^2)
  clustZhang        <- .psn_diag_cube(aW) / den_z
  signed_clustZhang <- .psn_diag_cube(W)  / den_z

  # Onnela: geometric mean of the three edge weights (sign-preserving cube root).
  aW13 <- aW^(1 / 3)
  W13  <- sign(W) * aW13
  Aon  <- (aW > 0) * 1
  kio  <- colSums(Aon); den_o <- kio * (kio - 1)
  clustOnnela        <- .psn_diag_cube(aW13) / den_o
  signed_clustOnnela <- .psn_diag_cube(W13)  / den_o

  # Barrat: strength-weighted, on absolute weights.
  Ab <- (aW > 0) * 1
  s  <- colSums(aW); kb <- colSums(Ab)
  num_b <- rowSums(aW * (Ab %*% Ab))
  clustBarrat <- num_b / (s * (kb - 1))

  out <- data.frame(node = labs,
                    clustWS = clustWS, signed_clustWS = signed_clustWS,
                    clustZhang = clustZhang, signed_clustZhang = signed_clustZhang,
                    clustOnnela = clustOnnela, signed_clustOnnela = signed_clustOnnela,
                    clustBarrat = clustBarrat,
                    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  # A node with fewer than two neighbours has a 0/0 coefficient: report 0, as
  # qgraph::clustcoef_auto does (an isolate closes no triangle).
  coef_cols <- setdiff(names(out), "node")
  out[coef_cols] <- lapply(out[coef_cols], function(v) { v[is.na(v)] <- 0; v })
  class(out) <- c("psychnet_clustering", "data.frame")
  out
}

# Unweighted global transitivity (3 * triangles / connected triads).
.psn_transitivity <- function(A) {
  diag(A) <- 0
  tri  <- sum(diag(A %*% A %*% A))                 # 6 * triangles
  k    <- rowSums(A)
  trip <- sum(k * (k - 1))                         # 2 * paths of length two
  if (trip < 0.5) 0 else tri / trip                # trip is a whole count
}

# Mean finite unweighted shortest-path length (BFS hop counts via repeated
# squaring of reachability is overkill; Floyd on unit weights is simplest).
.psn_aspl <- function(A) {
  diag(A) <- 0
  D <- .psn_floyd_warshall((A > 0) * 1, invert = FALSE)$D
  d <- D[upper.tri(D)]
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0L) NA_real_ else mean(d)
}

# One degree-preserving random graph by double-edge swaps on the edge list.
.psn_rewire <- function(A, swaps) {
  el <- which(A > 0 & upper.tri(A), arr.ind = TRUE)
  m <- nrow(el)
  if (m < 2L) return(A)
  # A literal for loop (documented exception to the no-loop rule): each swap is
  # proposed against the graph the previous swaps produced, so the chain is
  # genuinely sequential and cannot be vectorized.
  # Use a fixed number of proposals, retaining rejected moves. Stopping after
  # a fixed number of accepted swaps instead samples the jump chain, which
  # overweights graphs with more available swaps.
  for (attempt in seq_len(swaps)) {
    ij <- sample.int(m, 2L)
    a <- el[ij[1], 1]; b <- el[ij[1], 2]
    c <- el[ij[2], 1]; d <- el[ij[2], 2]
    # An undirected edge has no distinguished source endpoint. Randomize one
    # orientation so both reconnections can occur; otherwise the upper-triangle
    # edge list preserves arbitrary source/target degrees and biases the null.
    if (sample.int(2L, 1L) == 2L) { tmp <- c; c <- d; d <- tmp }
    if (length(unique(c(a, b, c, d))) < 4L) next         # shared endpoint
    if (A[a, d] || A[c, b]) next                         # would duplicate
    A[a, b] <- A[b, a] <- 0; A[c, d] <- A[d, c] <- 0
    A[a, d] <- A[d, a] <- 1; A[c, b] <- A[b, c] <- 1
    el[ij[1], ] <- c(a, d); el[ij[2], ] <- c(c, b)
  }
  A
}

#' Small-world index
#'
#' The Humphries & Gurney (2008) small-world index `sigma`: observed transitivity
#' and average shortest-path length compared with degree-preserving random graphs
#' (`sigma = (C/C_rand) / (L/L_rand)`; `sigma > 1` indicates small-worldness).
#' Computed on the binarised network. Ported from `qgraph::smallworldness`.
#'
#' @param x A `psychnet` object or a square weighted adjacency matrix.
#' @param n_rand Number of degree-preserving random graphs. Default 100.
#' @param seed Optional integer for reproducibility of the random graphs.
#' @return A tidy one-row `data.frame` (class `psychnet_smallworld`):
#'   `smallworldness`, `transitivity`, `aspl`, `transitivity_rand`, `aspl_rand`.
#'   When the random reference graphs have no triangles or no finite paths the
#'   index is undefined: `smallworldness` is `NA` and a warning of class
#'   `psychnet_undefined_statistic` is raised.
#' @references Humphries, M. D., & Gurney, K. (2008). Network
#'   'small-world-ness': A quantitative method for determining canonical network
#'   equivalence. *PLoS ONE*, 3(4), e0002051.
#'   \doi{10.1371/journal.pone.0002051}
#' @examples
#' S <- 0.4^abs(outer(1:8, 1:8, "-"))
#' net_smallworld(ebic_glasso(cor_matrix = S, n = 400), n_rand = 50, seed = 1)
#' @export
net_smallworld <- function(x, n_rand = 100L, seed = NULL) {
  if (inherits(x, "psychnet")) {
    if (isTRUE(x$directed))
      .psn_stop("net_smallworld() requires an undirected network.", call. = FALSE)
    W <- x$weights
  }
  else {
    if (!is.matrix(x) && !is.data.frame(x))
      .psn_stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    W <- as.matrix(x)
  }
  if (!is.numeric(W) || nrow(W) != ncol(W) || nrow(W) < 2L)
    .psn_stop("`x` must be a square numeric weighted adjacency matrix.", call. = FALSE)
  if (any(!is.finite(W)) || any(abs(W - t(W)) > 1e-8))
    .psn_stop("`x` must be a finite symmetric (undirected) matrix.", call. = FALSE)
  if (!is.numeric(n_rand) || length(n_rand) != 1L || !is.finite(n_rand) ||
      n_rand < 1 || n_rand != as.integer(n_rand))
    .psn_stop("`n_rand` must be a positive whole number.", call. = FALSE)
  n_rand <- as.integer(n_rand)
  # Seed the rewiring draws without leaving the caller's RNG stream disturbed.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  A <- (abs(W) > 0) * 1; diag(A) <- 0
  m <- sum(A) / 2
  C <- .psn_transitivity(A)
  L <- .psn_aspl(A)
  swaps <- max(10L, as.integer(10 * m))
  rand <- lapply(seq_len(as.integer(n_rand)), function(.) .psn_rewire(A, swaps))
  C_r <- mean(vapply(rand, .psn_transitivity, numeric(1)), na.rm = TRUE)
  L_r <- mean(vapply(rand, .psn_aspl, numeric(1)), na.rm = TRUE)
  # The index is undefined when the reference graphs have no triangles
  # (C_r = 0) or no finite paths (L_r = 0) -- e.g. a (near-)acyclic binarised
  # network. Return NA with a warning rather than Inf/NaN.
  tol <- sqrt(.Machine$double.eps)
  if (!is.finite(C_r) || C_r < tol || !is.finite(L_r) || L_r < tol) {
    .psn_warn("Small-world index undefined: random reference graphs have no ",
              "triangles or no finite paths; returning NA.",
              class = "psychnet_undefined_statistic")
    sigma <- NA_real_
  } else {
    sigma <- (C / C_r) / (L / L_r)
  }
  out <- data.frame(smallworldness = sigma, transitivity = C, aspl = L,
                    transitivity_rand = C_r, aspl_rand = L_r,
                    row.names = NULL)
  class(out) <- c("psychnet_smallworld", "data.frame")
  out
}
