# Community aggregation: collapse each community of items into a super-node.
# Two families (cf. Nestimate's mcml_pc):
#   * score methods (mean/median/sum/pca/factor/loadings) build one composite
#     COLUMN per community from the raw data and return a reduced data.frame -
#     re-estimate a macro network with psychnet() on the result.
#   * association methods (average/rv/canonical) have no per-row composite; they
#     return the community-by-community macro network directly as a `psychnet`.

.psn_score_methods <- c("mean", "median", "sum", "pca", "factor", "loadings")
.psn_assoc_methods <- c("average", "rv", "canonical")

# Sign-align a composite so it points the same way as its member items. A
# constant composite has no direction (and cor() would warn); it is returned
# as is. Member items are never constant (.as_numeric_matrix drops those).
.psn_align <- function(score, block) {
  if (!(stats::sd(score) > sqrt(.Machine$double.eps))) return(score)
  if (mean(stats::cor(score, block)) < 0) -score else score
}

# Item data for the data-taking framework verbs (net_aggregate, redundancy):
# label the ORIGINAL columns (from `labels`, else the column names, else
# V1..Vp), then reduce to the numeric, non-constant, complete-case matrix the
# estimators use (.as_numeric_matrix). Dropped incomplete rows are reported
# with a classed `psychnet_rows_dropped` warning -- listwise deletion is a
# decision the caller must see, not a silent one. Returns the matrix, the
# original labels and the kept labels, so a caller can tell which original
# columns were dropped.
#' @noRd
.psn_item_data <- function(data, labels = NULL, fn) {
  if (!is.matrix(data) && !is.data.frame(data))
    .psn_stop("`data` must be a data frame or numeric matrix.")
  if (!is.null(labels) && length(labels) != ncol(data))
    .psn_stop(sprintf(
      "`labels` must have one entry per column of `data`: it has %d, `data` has %d.",
      length(labels), ncol(data)))
  .psn_check_labels(labels, ncol(data))
  orig <- if (!is.null(labels)) as.character(labels) else colnames(data)
  if (is.null(orig)) orig <- paste0("V", seq_len(ncol(data)))
  if (anyNA(orig) || anyDuplicated(orig))
    .psn_stop("The columns of `data` need unique, non-missing names (or pass `labels`).")
  colnames(data) <- orig
  mat <- .as_numeric_matrix(data)
  n_drop <- nrow(data) - nrow(mat)
  if (n_drop > 0L)
    .psn_warn(sprintf(
      "%s(): dropped %d of %d rows with missing values (listwise deletion).",
      fn, n_drop, nrow(data)), class = "psychnet_rows_dropped")
  list(mat = mat, labels = orig, kept = colnames(mat))
}

# Community names in output order: a list in its element order, a factor in
# level order, a numeric vector in numeric order, and character labels sorted
# (the order of earlier releases). `resolved` is the character vector
# aligned to the kept columns.
#' @noRd
.psn_community_order <- function(communities, resolved) {
  present <- unique(resolved)
  if (is.list(communities)) {
    cn <- names(communities) %||% as.character(seq_along(communities))
    blank <- is.na(cn) | !nzchar(cn)
    cn[blank] <- as.character(seq_along(communities))[blank]
    return(intersect(cn, present))
  }
  if (is.factor(communities)) return(intersect(levels(communities), present))
  if (is.numeric(communities)) return(present[order(as.numeric(present))])
  sort(present)
}

# Columns whose community is explicitly NA (unnamed vector: NA at that
# position; named vector: an NA value) are left out of the aggregation before
# any preprocessing, so an ID column neither triggers listwise deletion nor
# needs to be numeric-free.
#' @noRd
.psn_left_out_columns <- function(data, communities, labels) {
  p <- ncol(data)
  if (is.list(communities) || !length(communities)) return(rep(FALSE, p))
  if (is.null(names(communities))) {
    # a length mismatch is reported by the community resolver
    if (length(communities) != p) return(rep(FALSE, p))
    return(is.na(communities))
  }
  orig <- if (!is.null(labels)) as.character(labels) else colnames(data)
  if (is.null(orig)) return(rep(FALSE, p))
  orig %in% names(communities)[is.na(communities)]
}

#' Aggregate a network's communities into super-nodes
#'
#' Collapses each community of items into a single super-node. `score` methods
#' build a composite column per community from the raw data and return a reduced
#' `data.frame` (re-estimate the macro network with [psychnet()] on it). `assoc`
#' methods summarise each community pair's multivariate association directly and
#' return the macro network as a `psychnet`.
#'
#' **Data handling.** Non-numeric columns are dropped, and so is any column
#' that is constant once incomplete rows are removed. Rows with a missing value
#' in any numeric column are dropped (listwise deletion), with a warning of
#' class `psychnet_rows_dropped` giving the count. `communities` is matched to
#' the columns of `data` *before* anything is dropped, so a dropped column must
#' not carry a community (give it `NA`, or leave it out of a named vector or
#' list): one that does is an error naming it.
#'
#' @param data A numeric data frame or matrix (rows = observations, columns =
#'   items).
#' @param communities Community membership of the columns of `data`, in one of
#'   three forms: an unnamed vector with one entry per column of `data`, in
#'   column order (`NA` for a column to leave out, e.g. an ID column); a vector
#'   named by item label; or a list of item-label vectors, one element per
#'   community (the element names are the community names). Every retained item
#'   must belong to exactly one community: a label listed twice, a label that is
#'   not a column, or a retained item without a community is an error of class
#'   `psychnet_bad_input`. Output communities are ordered by list order, factor
#'   level, numeric value, or (character labels) alphabetically. A column
#'   whose community is `NA` is left out before any preprocessing.
#' @param method Aggregation method. Score (return reduced data): `"mean"`
#'   (default), `"median"`, `"sum"`, `"pca"` (first principal component),
#'   `"factor"` (1-factor regression score; a community of fewer than 3 items
#'   uses its first principal component, since one factor is not identified
#'   there, and a community where `stats::factanal()` fails also falls back to
#'   its first principal component with a `psychnet_fit_failed` warning naming
#'   it), `"loadings"` (within-community connectivity-weighted mean).
#'   Association (return macro network): `"average"` (mean signed edge weight
#'   between communities), `"rv"` (Escoufier RV coefficient), `"canonical"`
#'   (first canonical correlation; a pair where `stats::cancor()` fails gets
#'   weight 0 with a `psychnet_fit_failed` warning naming it).
#' @param estimator Estimator used for the node-level network that `"loadings"`
#'   and `"average"` need (see [psychnet()]). Default `"glasso"`.
#' @param scale Standardise each item before forming score composites. Default
#'   `TRUE`.
#' @param labels Optional item labels, one per column of `data`, used in place
#'   of its column names (and matched by a named `communities`).
#' @param ... Passed to the estimator.
#' @return For a score `method`, a `data.frame` with one column per community
#'   (one row per retained observation). For an association `method`, a
#'   `psychnet` macro network among communities.
#' @examples
#' net_aggregate(SRL_Claude, communities = c(1, 1, 2, 2, 2))            # reduced data
#' net_aggregate(SRL_Claude, communities = c(1, 1, 2, 2, 2), method = "rv")  # macro net
#' @export
net_aggregate <- function(data, communities, method = "mean",
                          estimator = "glasso", scale = TRUE, labels = NULL,
                          ...) {
  method <- .psn_match_arg(method, c(.psn_score_methods, .psn_assoc_methods))
  .psn_check_flag(scale, "scale")
  left_out <- .psn_left_out_columns(data, communities, labels)
  if (any(left_out)) {
    data <- data[, !left_out, drop = FALSE]
    if (!is.null(labels)) labels <- labels[!left_out]
    communities <- if (is.null(names(communities))) communities[!left_out]
                   else communities[!is.na(communities)]
  }
  prep <- .psn_item_data(data, labels, "net_aggregate")
  mat <- prep$mat
  labs <- prep$kept

  # Resolve communities against the ORIGINAL columns, then keep the retained.
  comm_all <- .psn_resolve_communities(communities, prep$labels, "item")
  dropped <- setdiff(prep$labels, labs)
  carry <- dropped[!is.na(comm_all[dropped])]
  if (length(carry))
    .psn_stop(sprintf(paste0(
      "`communities` assigns a community to column(s) dropped from `data` ",
      "(non-numeric, or constant once incomplete rows are removed): %s. ",
      "Remove them from `data`, or give them community NA."),
      paste(carry, collapse = ", ")),
      class = c("psychnet_dropped_column", "psychnet_bad_input"))
  communities_kept <- unname(comm_all[labs])
  if (anyNA(communities_kept))
    .psn_stop("`communities` gives no community for item(s): ",
              paste(labs[is.na(communities_kept)], collapse = ", "), ".")
  ucomm <- .psn_community_order(communities, communities_kept)
  members <- lapply(ucomm, function(cc) which(communities_kept == cc))
  names(members) <- ucomm

  Z <- if (isTRUE(scale)) scale(mat) else mat
  tol <- sqrt(.Machine$double.eps)

  # ---- score methods: one composite column per community -> reduced data -----
  if (method %in% .psn_score_methods) {
    # loadings needs the node-level weights for within-community strength.
    Wnode <- if (method == "loadings")
      abs(psychnet(mat, method = estimator, labels = labs, ...)$weights) else NULL

    first_pc <- function(blk) stats::prcomp(blk, rank. = 1)$x[, 1]
    comp <- function(idx, cname) {
      blk <- Z[, idx, drop = FALSE]
      if (length(idx) == 1L) return(as.numeric(blk))
      switch(method,
        mean   = rowMeans(blk),
        median = apply(blk, 1L, stats::median),
        sum    = rowSums(blk),
        pca    = .psn_align(first_pc(blk), blk),
        factor = {
          # One factor is not identified on < 3 items: the documented fallback
          # is the first principal component. Otherwise a factanal() failure
          # is surfaced (classed warning naming the community), not hidden.
          sc <- if (length(idx) < 3L) first_pc(blk) else tryCatch(
            stats::factanal(blk, factors = 1, scores = "regression")$scores[, 1],
            error = function(e) {
              .psn_warn(sprintf(paste0(
                "net_aggregate(method = \"factor\"): factanal() failed for ",
                "community '%s' (%s); using its first principal component."),
                cname, conditionMessage(e)), class = "psychnet_fit_failed")
              first_pc(blk)
            })
          .psn_align(sc, blk)
        },
        loadings = {
          w <- vapply(idx, function(j) sum(Wnode[j, setdiff(idx, j)]), numeric(1))
          if (sum(w) < tol) w <- rep(1, length(idx))
          as.numeric(blk %*% (w / sum(w)))
        })
    }
    out <- vapply(ucomm, function(cc) as.numeric(comp(members[[cc]], cc)),
                  numeric(nrow(Z)))
    out <- matrix(out, nrow = nrow(Z))
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    names(out) <- ucomm
    return(out)
  }

  # ---- association methods: community-by-community macro network -------------
  k <- length(ucomm)
  M <- matrix(0, k, k, dimnames = list(ucomm, ucomm))
  directed <- FALSE
  if (method == "average") {
    node_fit <- psychnet(mat, method = estimator, labels = labs, ...)
    W <- node_fit$weights
    directed <- isTRUE(node_fit$directed)
    ab <- which(row(M) != col(M), arr.ind = TRUE)
    M[ab] <- vapply(seq_len(nrow(ab)), function(r)
      mean(W[members[[ab[r, 1L]]], members[[ab[r, 2L]]]]), numeric(1))
  } else {
    S <- stats::cov(mat)
    ab <- which(upper.tri(M), arr.ind = TRUE)
    vals <- vapply(seq_len(nrow(ab)), function(r) {
      ia <- members[[ab[r, 1L]]]; ib <- members[[ab[r, 2L]]]
      if (method == "rv") {
        Sab <- S[ia, ib, drop = FALSE]
        Saa <- S[ia, ia, drop = FALSE]; Sbb <- S[ib, ib, drop = FALSE]
        return(sum(Sab^2) / sqrt(sum(Saa^2) * sum(Sbb^2)))
      }
      # canonical: cancor errors (not just warns) on rank-deficient /
      # wider-than-tall blocks. Such a pair gets weight 0, said out loud.
      tryCatch(
        stats::cancor(mat[, ia, drop = FALSE], mat[, ib, drop = FALSE])$cor[1],
        error = function(e) {
          .psn_warn(sprintf(paste0(
            "net_aggregate(method = \"canonical\"): cancor() failed for ",
            "communities '%s' and '%s' (%s); their edge is set to 0."),
            ucomm[ab[r, 1L]], ucomm[ab[r, 2L]], conditionMessage(e)),
            class = "psychnet_fit_failed")
          0
        })
    }, numeric(1))
    M[ab] <- vals
    M[ab[, c(2L, 1L), drop = FALSE]] <- vals
  }
  .new_psychnet(M, ucomm, method = paste0("aggregate_", method),
                directed = directed, n_obs = nrow(mat))
}
