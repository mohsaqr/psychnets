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
#' @references Costantini, G., & Perugini, M. (2014); Watts & Strogatz (1998);
#'   Zhang & Horvath (2005); Onnela et al. (2005); Barrat et al. (2004).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' net_clustering(ebic_glasso(cor_matrix = S, n = 400))
#' @export
net_clustering <- function(x, labels = NULL) {
  if (inherits(x, "psychnet")) { W <- x$weights; labs <- x$nodes$label }
  else {
    if (!is.matrix(x) && !is.data.frame(x))
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    W <- as.matrix(x)
    if (!is.numeric(W) || nrow(W) != ncol(W))
      stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    labs <- if (!is.null(labels)) labels else colnames(W)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(W)))
  }
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
  for (nm in names(out)[-1]) out[[nm]][is.na(out[[nm]])] <- 0  # isolates -> 0
  class(out) <- c("psychnet_clustering", "data.frame")
  out
}

# Unweighted global transitivity (3 * triangles / connected triads).
.psn_transitivity <- function(A) {
  diag(A) <- 0
  tri  <- sum(diag(A %*% A %*% A))                 # 6 * triangles
  k    <- rowSums(A)
  trip <- sum(k * (k - 1))                         # 2 * paths of length two
  if (trip == 0) 0 else tri / trip
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
  for (s in seq_len(swaps)) {
    ij <- sample.int(m, 2L)
    a <- el[ij[1], 1]; b <- el[ij[1], 2]
    c <- el[ij[2], 1]; d <- el[ij[2], 2]
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
#' @references Humphries, M. D., & Gurney, K. (2008). *PLoS ONE*, 3(4), e0002051.
#' @examples
#' S <- 0.4^abs(outer(1:8, 1:8, "-"))
#' net_smallworld(ebic_glasso(cor_matrix = S, n = 400), n_rand = 50, seed = 1)
#' @export
net_smallworld <- function(x, n_rand = 100L, seed = NULL) {
  if (inherits(x, "psychnet")) W <- x$weights
  else {
    if (!is.matrix(x) && !is.data.frame(x))
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    W <- as.matrix(x)
  }
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
  if (!is.finite(C_r) || C_r == 0 || !is.finite(L_r) || L_r == 0) {
    warning("Small-world index undefined: random reference graphs have no ",
            "triangles or no finite paths; returning NA.", call. = FALSE)
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
