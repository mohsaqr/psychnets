# Reference-side code shared by the equivalence tests, in base R. Nothing here calls psychnets: each
# function shapes a reference's output, or builds an input both sides read.

#' The upper triangle of a square matrix as one row per node pair
#'
#' psychnets' as.data.frame() names an undirected edge by its row node (`from`) and column node (`to`)
#' in the upper triangle, column by column; this writes a reference matrix the same way.
#'
#' @param m A square matrix with dimnames.
#' @param name The name of the value column.
#' @return A data frame `from`, `to`, `<name>`.
ref_edges <- function(m, name = "weight") {
  m <- as.matrix(m)
  stopifnot("`m` must be square with names" = nrow(m) == ncol(m) && !is.null(colnames(m)))
  if (is.null(rownames(m))) rownames(m) <- colnames(m)
  upper <- upper.tri(m)
  out <- data.frame(from = rownames(m)[row(m)[upper]], to = colnames(m)[col(m)[upper]],
                    value = m[upper], stringsAsFactors = FALSE)
  names(out)[[3L]] <- name
  out
}

#' Several square matrices side by side, one row per node pair (upper triangle)
#'
#' @param ... Named matrices of one shape.
#' @return A data frame `from`, `to` and one column per matrix.
ref_edge_table <- function(...) {
  tables <- list(...)
  Reduce(\(left, right) merge(left, right, by = c("from", "to"), sort = FALSE),
         Map(ref_edges, tables, names(tables)))
}

#' All ordered node pairs of a square matrix (for a directed network)
#'
#' @param m A square matrix with dimnames; rows are `from`, columns are `to`.
#' @param name The name of the value column.
#' @return A data frame `from`, `to`, `<name>`, the diagonal included.
ref_cells <- function(m, name) {
  m <- as.matrix(m)
  out <- data.frame(from = rownames(m)[as.vector(row(m))], to = colnames(m)[as.vector(col(m))],
                    value = as.vector(m), stringsAsFactors = FALSE)
  names(out)[[3L]] <- name
  out
}

#' A weighted adjacency matrix read from a tab-separated file with row names
#'
#' @param path The file.
#' @return A numeric matrix with dimnames.
ref_read_matrix <- function(path) {
  as.matrix(utils::read.delim(path, row.names = 1L, check.names = FALSE))
}

#' Two groups of simulated multivariate data for the network comparison test
#'
#' Group 1 has one common factor (loading 0.4); group 2 a weaker factor plus a V1 -> V2 path. The
#' generator of the package's own NCT tests (tests/testthat/helper-nct.R), copied here.
#'
#' @param seed The seed.
#' @param n,p Rows per group and variables.
#' @return A list of two numeric matrices `a` and `b`.
ref_nct_groups <- function(seed, n = 120, p = 5) {
  set.seed(seed)
  a <- matrix(stats::rnorm(n * p), n, p)
  a <- a + 0.4 * matrix(rowSums(a), n, p)
  b <- matrix(stats::rnorm(n * p), n, p)
  b <- b + 0.15 * matrix(rowSums(b), n, p)
  b[, 2] <- b[, 2] + 0.5 * b[, 1]
  colnames(a) <- colnames(b) <- paste0("V", seq_len(p))
  list(a = a, b = b)
}

#' A simulated panel: a stationary random-intercept VAR(1), wide form with its design matrix
#'
#' y_it = mu + b_i + w_it, w_t = B w_{t-1} + z_t; the generator of the package's panel tests
#' (tests/testthat/helper-panelgvar.R), copied here. The RNG state is restored on exit.
#'
#' @param n,p,waves Persons, variables, waves.
#' @param seed The seed.
#' @return A list: `data` (one row per person, columns `id` and `V<j>_w<t>`) and `design` (a p x waves
#'   matrix of those column names, as psychonetrics::panelgvar(vars =) takes it).
ref_panel <- function(n = 400, p = 3, waves = 4, seed = 1) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  b <- diag(0.4, p)
  b[2, 1] <- 0.25
  b[3, 2] <- -0.2
  sz <- diag(p)
  sz[1, 2] <- sz[2, 1] <- 0.3
  sb <- diag(0.6, p)
  sb[1, 3] <- sb[3, 1] <- 0.25
  s0 <- matrix(solve(diag(p * p) - b %x% b, as.vector(sz)), p, p)
  trait <- matrix(stats::rnorm(n * p), n) %*% chol(sb)
  w <- matrix(stats::rnorm(n * p), n) %*% chol(s0)
  states <- Reduce(\(prev, t) prev %*% t(b) + matrix(stats::rnorm(n * p), n) %*% chol(sz),
                   seq_len(waves - 1L), w, accumulate = TRUE)
  vars <- paste0("V", seq_len(p))
  columns <- lapply(seq_len(waves), \(t) {
    y <- sweep(states[[t]] + trait, 2L, seq_len(p), "+")
    colnames(y) <- paste0(vars, "_w", t)
    y
  })
  design <- matrix(paste0(rep(vars, waves), "_w", rep(seq_len(waves), each = p)), nrow = p,
                   dimnames = list(vars, NULL))
  list(data = data.frame(id = seq_len(n), do.call(cbind, columns)), design = design)
}

#' psychonetrics' parameter table keyed as psychnets' (matrix, row, col by variable name)
#'
#' @param parameters The `@parameters` slot of a psychonetrics model.
#' @param vars The variable names, in the model's order.
#' @return A data frame `matrix`, `row`, `col`, `estimate`, `se`, `fixed` (one row per parameter
#'   psychnets reports: the free ones of the saturated model, and the fixed-to-zero ones of a pruned
#'   model).
ref_psychonetrics_parameters <- function(parameters, vars) {
  # The diagonal of an omega (partial correlation) matrix is not a parameter: psychonetrics lists it as fixed.
  kept <- parameters[!(startsWith(parameters$matrix, "omega") & parameters$row == parameters$col), , drop = FALSE]
  data.frame(matrix = kept$matrix, row = vars[kept$row],
             col = ifelse(kept$matrix == "mu", "(mean)", vars[kept$col]),
             estimate = kept$est, se = kept$se, fixed = kept$fixed, stringsAsFactors = FALSE)
}
