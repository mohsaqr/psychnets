# Network Comparison Test (van Borkulo et al. 2022), ported from Nestimate and
# adapted to psychnet: the EBIC-glasso estimator is the clean-room pure-R solver
# and the nearest-correlation projection is base R (no Matrix::nearPD).

# Nearest correlation matrix by one eigenvalue-clamp projection (Higham 2002,
# single step): clamp negative eigenvalues to a small positive floor and
# renormalize to a unit diagonal. Enough to make a barely-non-PD subsample
# correlation matrix usable by the glasso.
#' @noRd
.nearest_pd_cor <- function(S, eig_tol = 1e-8) {
  dn <- dimnames(S)                      # eigen() drops names; restore them so a
  S <- (S + t(S)) / 2                    # projected correlation matrix stays named
  e <- eigen(S, symmetric = TRUE)
  vals <- pmax(e$values, eig_tol)
  S2 <- e$vectors %*% (vals * t(e$vectors))
  d <- sqrt(diag(S2)); d[d < 1e-12] <- 1
  S2 <- S2 / outer(d, d)
  S2 <- (S2 + t(S2)) / 2
  dimnames(S2) <- dn
  S2
}

#' Network Comparison Test
#'
#' Permutation test for whether two groups' Gaussian graphical models differ,
#' on three invariants: global strength (`M`), maximum edge difference (`S`),
#' and per-edge differences (`E`). Networks are EBIC graphical lassos (clean-room
#' pure R). Equivalent in purpose to `NetworkComparisonTest::NCT()`.
#'
#' @param data1,data2 Numeric data frames/matrices with the same columns.
#' @param iter Number of permutations. Default 1000.
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param paired Logical; within-row swapping for paired designs. Default FALSE.
#' @param abs Logical; compare absolute edge weights. Default TRUE.
#' @param weighted Logical; if FALSE, binarize networks first. Default TRUE.
#' @param p_adjust Multiple-comparison adjustment for per-edge p-values
#'   (any [stats::p.adjust] method). Default `"none"`.
#' @return An object of class `psychnet_nct` with `$nw1`, `$nw2`, and `$M`,
#'   `$S`, `$E` (each `observed`, `perm`, `p_value`); `$E` also carries
#'   `edge_names`, a `from`/`to` data frame aligned to the per-edge vector.
#' @examples
#' set.seed(1)
#' a <- matrix(stats::rnorm(150 * 5), 150, 5)
#' b <- matrix(stats::rnorm(150 * 5), 150, 5)
#' colnames(a) <- colnames(b) <- paste0("V", 1:5)
#' fit <- net_compare(a, b, iter = 25)
#' fit
#' @export
net_compare <- function(data1, data2 = NULL, iter = 1000L, gamma = 0.5,
                paired = FALSE, abs = TRUE, weighted = TRUE, p_adjust = "none") {
  # Group object -> compare two of its levels' cross-sectional data. `data2`
  # names the two levels (defaulting to the two levels when there are exactly 2).
  if (inherits(data1, "psychnet_group")) {
    if (!identical(attr(data1, "source"), "data"))
      stop("net_compare() supports group mode for cross-sectional data only.",
           call. = FALSE)
    subs <- attr(data1, "subsets")
    pair <- if (is.null(data2)) names(subs) else as.character(data2)
    if (length(pair) != 2L)
      stop("Specify two group levels to compare, e.g. data2 = c(\"A\", \"B\").",
           call. = FALSE)
    miss <- setdiff(pair, names(subs))
    if (length(miss))
      stop("Group level(s) not found: ", paste(miss, collapse = ", "),
           call. = FALSE)
    cl <- attr(data1, "call")
    return(net_compare(subs[[pair[1L]]], subs[[pair[2L]]], iter = iter,
                       gamma = if (!is.null(cl$gamma)) cl$gamma else gamma,
                       paired = paired, abs = abs, weighted = weighted,
                       p_adjust = p_adjust))
  }
  if (is.null(data2))
    stop("`data2` is required (a second data set, or two group levels when ",
         "`data1` is a psychnet_group).", call. = FALSE)
  # Keep only numeric columns (a stored group subset may carry a stray
  # character/id column); est() runs stats::cor() directly, which would error on
  # a non-numeric matrix. Zero-variance numeric columns are retained -- est()
  # already neutralises their NA correlations.
  .numeric_only <- function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    as.matrix(d[, vapply(d, is.numeric, logical(1)), drop = FALSE])
  }
  data1 <- .numeric_only(data1)
  data2 <- .numeric_only(data2)
  stopifnot(ncol(data1) == ncol(data2), iter >= 1L,
            is.logical(paired), is.logical(abs), is.logical(weighted))
  if (!is.null(colnames(data1)) && !is.null(colnames(data2)) &&
      !identical(colnames(data1), colnames(data2))) {
    stop("`data1` and `data2` must have the same columns in the same order.",
         call. = FALSE)
  }
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)

  iter <- as.integer(iter)
  n1 <- nrow(data1); n2 <- nrow(data2)
  if (paired && n1 != n2) {
    stop("A paired comparison requires the two groups to have equal size.",
         call. = FALSE)
  }
  dataall <- rbind(data1, data2)

  est <- function(x) {
    cx <- suppressWarnings(stats::cor(x))
    # a column that is constant in this (sub)sample has undefined correlations;
    # treat it as unassociated rather than letting NA reach the eigen solver
    cx[is.na(cx)] <- 0
    diag(cx) <- 1
    cor_x <- .nearest_pd_cor(cx)
    ebic_glasso(cor_matrix = cor_x, n = nrow(x), gamma = gamma)$weights
  }
  binarize <- function(m) (m != 0) * 1

  nw1 <- est(data1); nw2 <- est(data2)
  ut <- upper.tri(nw1)
  # Edge labels for the per-edge `E` vector, in the same upper-triangle order as
  # E$observed (column-major over the TRUE cells of `ut`). Lets a caller map each
  # E p-value back to its (from, to) node pair without reconstructing the order.
  labs <- colnames(data1)
  if (is.null(labs)) labs <- paste0("V", seq_len(ncol(data1)))
  ij <- which(ut, arr.ind = TRUE)
  edge_names <- data.frame(from = labs[ij[, 1L]], to = labs[ij[, 2L]],
                           stringsAsFactors = FALSE, row.names = NULL)
  if (!weighted) { nw1 <- binarize(nw1); nw2 <- binarize(nw2) }

  M_obs <- if (abs) base::abs(sum(base::abs(nw1[ut])) - sum(base::abs(nw2[ut])))
           else     base::abs(sum(nw1[ut]) - sum(nw2[ut]))
  diff_real <- base::abs(nw1 - nw2)
  S_obs <- max(diff_real[ut])
  E_obs <- diff_real[ut]
  n_edges <- length(E_obs)

  # Permutation null. A literal for loop is the right tool here: vectorising
  # would materialise iter x n_edges intermediate networks.
  M_perm <- numeric(iter); S_perm <- numeric(iter)
  E_perm <- matrix(0, iter, n_edges)
  for (i in seq_len(iter)) {
    if (paired) {
      s <- sample(c(1L, 2L), n1, replace = TRUE)
      x1p <- rbind(data1[s == 1L, , drop = FALSE], data2[s == 2L, , drop = FALSE])
      x2p <- rbind(data2[s == 1L, , drop = FALSE], data1[s == 2L, , drop = FALSE])
    } else {
      s <- sample(seq_len(n1 + n2), n1, replace = FALSE)
      x1p <- dataall[s, , drop = FALSE]; x2p <- dataall[-s, , drop = FALSE]
    }
    r1 <- est(x1p); r2 <- est(x2p)
    if (!weighted) { r1 <- binarize(r1); r2 <- binarize(r2) }
    M_perm[i] <- if (abs) base::abs(sum(base::abs(r1[ut])) - sum(base::abs(r2[ut])))
                 else     base::abs(sum(r1[ut]) - sum(r2[ut]))
    diff_perm <- base::abs(r1 - r2)
    S_perm[i] <- max(diff_perm[ut])
    E_perm[i, ] <- diff_perm[ut]
  }

  M_pval <- (sum(M_perm >= M_obs) + 1) / (iter + 1)
  S_pval <- (sum(S_perm >= S_obs) + 1) / (iter + 1)
  E_pval <- (colSums(E_perm >= matrix(E_obs, iter, n_edges, byrow = TRUE)) + 1) /
            (iter + 1)
  if (p_adjust != "none") E_pval <- stats::p.adjust(E_pval, method = p_adjust)

  structure(
    list(nw1 = nw1, nw2 = nw2,
         M = list(observed = M_obs, perm = M_perm, p_value = M_pval),
         S = list(observed = S_obs, perm = S_perm, p_value = S_pval),
         E = list(observed = E_obs, perm = E_perm, p_value = E_pval,
                  edge_names = edge_names),
         n_iter = iter, paired = paired,
         params = list(gamma = gamma, abs = abs, weighted = weighted,
                       p_adjust = p_adjust)),
    class = "psychnet_nct"
  )
}

#' Print a Network Comparison Test
#'
#' @param x A `psychnet_nct` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.psychnet_nct <- function(x, ...) {
  cat("Network Comparison Test (", x$n_iter, " permutations)\n", sep = "")
  cat(sprintf("  Global strength (M): observed %.3f, p = %.3f\n",
              x$M$observed, x$M$p_value))
  cat(sprintf("  Network structure (S): observed %.3f, p = %.3f\n",
              x$S$observed, x$S$p_value))
  invisible(x)
}
