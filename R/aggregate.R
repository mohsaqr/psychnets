# Community aggregation: collapse each community of items into a super-node.
# Two families (cf. Nestimate's mcml_pc):
#   * score methods (mean/median/sum/pca/factor/loadings) build one composite
#     COLUMN per community from the raw data and return a reduced data.frame -
#     re-estimate a macro network with psychnet() on the result.
#   * association methods (average/rv/canonical) have no per-row composite; they
#     return the community-by-community macro network directly as a `psychnet`.

.psn_score_methods <- c("mean", "median", "sum", "pca", "factor", "loadings")
.psn_assoc_methods <- c("average", "rv", "canonical")

# Sign-align a composite so it points the same way as its member items.
.psn_align <- function(score, block) {
  r <- suppressWarnings(stats::cor(score, block))
  if (mean(r, na.rm = TRUE) < 0) -score else score
}

#' Aggregate a network's communities into super-nodes
#'
#' Collapses each community of items into a single super-node. `score` methods
#' build a composite column per community from the raw data and return a reduced
#' `data.frame` (re-estimate the macro network with [psychnet()] on it). `assoc`
#' methods summarise each community pair's multivariate association directly and
#' return the macro network as a `psychnet`.
#'
#' @param data A numeric data frame or matrix (rows = observations, columns =
#'   items).
#' @param communities Community membership, one entry per item (column): a vector
#'   aligned to the columns, or a named vector / list keyed by item label.
#' @param method Aggregation method. Score (return reduced data): `"mean"`
#'   (default), `"median"`, `"sum"`, `"pca"` (first principal component),
#'   `"factor"` (1-factor score, falls back to PCA for communities of < 3
#'   items), `"loadings"` (within-community connectivity-weighted mean).
#'   Association (return macro network): `"average"` (mean signed edge weight
#'   between communities), `"rv"` (Escoufier RV coefficient), `"canonical"`
#'   (first canonical correlation).
#' @param estimator Estimator used for the node-level network that `"loadings"`
#'   and `"average"` need (see [psychnet()]). Default `"glasso"`.
#' @param scale Standardise each item before forming score composites. Default
#'   `TRUE`.
#' @param labels Optional item labels (used when `data` has no column names).
#' @param ... Passed to the estimator.
#' @return For a score `method`, a `data.frame` with one column per community
#'   (one row per observation). For an association `method`, a `psychnet` macro
#'   network among communities.
#' @examples
#' net_aggregate(SRL_Claude, communities = c(1, 1, 2, 2, 2))            # reduced data
#' net_aggregate(SRL_Claude, communities = c(1, 1, 2, 2, 2), method = "rv")  # macro net
#' @export
net_aggregate <- function(data, communities, method = "mean",
                          estimator = "glasso", scale = TRUE, labels = NULL,
                          ...) {
  method <- match.arg(method, c(.psn_score_methods, .psn_assoc_methods))
  mat <- .as_numeric_matrix(data)
  p <- ncol(mat)
  labs <- if (!is.null(labels)) labels else colnames(mat)
  if (is.null(labs)) labs <- paste0("V", seq_len(p))
  colnames(mat) <- labs

  # Resolve communities to a character vector aligned to columns.
  if (is.list(communities)) {
    st <- utils::stack(communities)
    communities <- stats::setNames(as.character(st$ind),
                                   as.character(st$values))[labs]
  } else if (!is.null(names(communities))) {
    communities <- communities[labs]
  }
  communities <- as.character(communities)
  if (length(communities) != p)
    stop("`communities` must have one entry per item.", call. = FALSE)
  if (anyNA(communities))
    stop("`communities` could not be aligned to every item.", call. = FALSE)
  ucomm <- sort(unique(communities))
  members <- lapply(ucomm, function(cc) which(communities == cc))
  names(members) <- ucomm

  Z <- if (isTRUE(scale)) scale(mat) else mat

  # ---- score methods: one composite column per community -> reduced data -----
  if (method %in% .psn_score_methods) {
    # loadings needs the node-level weights for within-community strength.
    Wnode <- if (method == "loadings")
      abs(psychnet(mat, method = estimator, labels = labs, ...)$weights) else NULL

    comp <- function(idx) {
      blk <- Z[, idx, drop = FALSE]
      if (length(idx) == 1L) return(as.numeric(blk))
      switch(method,
        mean   = rowMeans(blk),
        median = apply(blk, 1L, stats::median),
        sum    = rowSums(blk),
        pca    = .psn_align(stats::prcomp(blk, rank. = 1)$x[, 1], blk),
        factor = {
          sc <- tryCatch(
            stats::factanal(blk, factors = 1, scores = "regression")$scores[, 1],
            error = function(e) NULL)
          if (is.null(sc)) sc <- stats::prcomp(blk, rank. = 1)$x[, 1]
          .psn_align(sc, blk)
        },
        loadings = {
          w <- vapply(idx, function(j) sum(Wnode[j, setdiff(idx, j)]), numeric(1))
          if (sum(w) == 0) w <- rep(1, length(idx))
          as.numeric(blk %*% (w / sum(w)))
        })
    }
    out <- vapply(members, comp, numeric(nrow(Z)))
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    names(out) <- ucomm
    return(out)
  }

  # ---- association methods: community-by-community macro network -------------
  k <- length(ucomm)
  M <- matrix(0, k, k, dimnames = list(ucomm, ucomm))
  if (method == "average") {
    W <- psychnet(mat, method = estimator, labels = labs, ...)$weights
    for (a in seq_len(k)) for (b in seq_len(k)) if (a != b)
      M[a, b] <- mean(W[members[[a]], members[[b]]])
  } else {
    S <- stats::cov(mat)
    for (a in seq_len(k)) for (b in seq_len(k)) {
      if (a >= b) next
      ia <- members[[a]]; ib <- members[[b]]
      if (method == "rv") {
        Sab <- S[ia, ib, drop = FALSE]
        Saa <- S[ia, ia, drop = FALSE]; Sbb <- S[ib, ib, drop = FALSE]
        val <- sum(Sab^2) / sqrt(sum(Saa^2) * sum(Sbb^2))
      } else {                                    # canonical
        # cancor errors (not just warns) on rank-deficient / wider-than-tall
        # blocks; treat an unresolvable pair as 0 association rather than
        # aborting the whole aggregation.
        val <- tryCatch(
          suppressWarnings(stats::cancor(mat[, ia, drop = FALSE],
                                         mat[, ib, drop = FALSE])$cor[1]),
          error = function(e) NA_real_)
        if (is.na(val)) val <- 0
      }
      M[a, b] <- M[b, a] <- val
    }
  }
  .new_psychnet(M, ucomm, method = paste0("aggregate_", method),
                directed = FALSE, n_obs = nrow(mat))
}
