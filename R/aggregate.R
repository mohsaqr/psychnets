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

#' Aggregate communities of items into a network between communities
#'
#' @description
#' `net_aggregate()` builds a network between communities of items, such as the
#' items of several scales, in which each community is a single node. It offers
#' two routes: the items of each community are combined into a composite score,
#' and a network is then estimated among the composites; or the association
#' between two communities is summarized directly from their items.
#'
#' The score methods are the mean, median or sum of the items, the first
#' principal component, a one-factor score, or a mean in which each item is
#' weighted by its connections to the other items of its community in the
#' item-level network (`"loadings"`). The association methods are the mean
#' signed edge weight between the items of two communities in the item-level
#' network (`"average"`), Escoufier's RV coefficient (`"rv"`), or the first
#' canonical correlation (`"canonical"`).
#'
#' A score method returns a data frame with one composite column per community,
#' ready for [psychnet()]; an association method returns the network between
#' communities directly.
#'
#' The centrality of single items with respect to the communities is computed
#' by [net_bridge()]; pairs of nearly identical items can be found with
#' [redundancy()] before aggregating.
#'
#' @details
#' Non-numeric columns are dropped, and so is any column that is constant once
#' incomplete rows are removed. Rows with a missing value in any numeric
#' column are dropped (listwise deletion). `communities` is matched to the
#' columns of `data` before anything is dropped, so a dropped column must not
#' carry a community: give it `NA`, or leave it out of a named vector or list.
#'
#' A community of one item is that item (standardized when `scale = TRUE`).
#' Principal-component and factor scores are signed to point in the direction
#' of the community's items. A factor score needs at least three items, so a
#' smaller community, or one for which `stats::factanal()` fails, uses its
#' first principal component instead. In
#' `"loadings"`, each item's weight is the sum of its absolute edge weights to
#' the other items of its community (equal weights when these are all zero).
#' A pair of communities for which `stats::cancor()` fails gets weight 0 under
#' `"canonical"`.
#'
#' @param data A numeric data frame or matrix with one row per observation and
#'   one column per item.
#' @param communities The community of each column of `data`, in one of three
#'   forms: an unnamed vector with one entry per column, in column order (`NA`
#'   for a column to leave out, such as an identifier); a vector named by item
#'   label; or a list of item-label vectors, one per community, whose names are
#'   the community names. Every retained item must belong to exactly one
#'   community. Communities are ordered by list order, factor level, numeric
#'   value, or alphabetically for character labels.
#' @param method Aggregation method. Score methods (return a data frame):
#'   `"mean"` (default), `"median"`, `"sum"`, `"pca"`, `"factor"` or
#'   `"loadings"`. Association methods (return a network): `"average"`, `"rv"`
#'   or `"canonical"`.
#' @param estimator Estimator of the item-level network used by `"loadings"`
#'   and `"average"`, as in [psychnet()]. Default `"glasso"`.
#' @param scale If `TRUE` (default), standardize each item before forming
#'   composite scores.
#' @param labels Optional item labels, one per column of `data`, used in place
#'   of its column names (and matched by a named `communities`).
#' @param ... Further arguments for the estimator.
#' @return For a score `method`, a data frame with one row per retained
#'   observation and one column per community, named by community. For an
#'   association `method`, a fitted network ([psychnet-class]) with one node
#'   per community; it is directed only for `"average"` with a directed
#'   estimator.
#' @section Conditions:
#'   `psychnet_bad_input` for an invalid argument, or when `communities` lists
#'   an item twice, names an item that is not a column, or leaves a retained
#'   item without a community; `psychnet_dropped_column` (an error, also
#'   `psychnet_bad_input`) when a dropped column carries a community;
#'   `psychnet_rows_dropped` (warning)
#'   when incomplete rows are dropped; `psychnet_constant_variable` (warning)
#'   when a constant column is dropped; `psychnet_dropped_columns` (a message)
#'   when non-numeric columns are dropped; `psychnet_fit_failed` (warning)
#'   when `stats::factanal()` or `stats::cancor()` fails for a community,
#'   naming it. With `"loadings"` or `"average"`, the estimator raises its own
#'   conditions.
#' @seealso [psychnet()] to estimate a network from the composite scores;
#'   [net_bridge()] for bridge centrality; [redundancy()] for redundant items.
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
