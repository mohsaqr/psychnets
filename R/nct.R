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

#' @noRd
.nct_bad_input <- function(msg) {
  stop(errorCondition(msg, class = c("psychnet_bad_input", "psychnet_error"),
                      call = NULL))
}

#' @noRd
.nct_degenerate <- function(msg) {
  stop(errorCondition(msg, class = c("psychnet_degenerate_sample",
                                     "psychnet_error"), call = NULL))
}

#' @noRd
.nct_validate <- function(iter, gamma, paired, abs, weighted) {
  flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  ok <- c(
    "`iter` must be a single whole number >= 1" =
      is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 1 &&
      base::abs(iter - round(iter)) < sqrt(.Machine$double.eps),
    "`gamma` must be a single non-negative number" =
      is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) && gamma >= 0,
    "`paired` must be TRUE or FALSE"   = flag(paired),
    "`abs` must be TRUE or FALSE"      = flag(abs),
    "`weighted` must be TRUE or FALSE" = flag(weighted)
  )
  if (!all(ok)) .nct_bad_input(paste(names(ok)[!ok], collapse = "; "))
  invisible(TRUE)
}

# Keep the numeric columns of each group and check that the two groups hold the
# same variables in the same order. A position-wise comparison of differently
# ordered columns would test the wrong pairs, so a mismatch is an error, not a
# silent reorder. Unnamed input is taken by position and labelled V1..Vp.
#' @noRd
.nct_columns <- function(data1, data2) {
  numeric_only <- function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    as.matrix(d[, vapply(d, is.numeric, logical(1)), drop = FALSE])
  }
  d1 <- numeric_only(data1)
  d2 <- numeric_only(data2)
  if (ncol(d1) < 2L || ncol(d2) < 2L)
    .nct_bad_input("`data1` and `data2` each need at least 2 numeric variables.")
  # A value is either observed (finite) or missing (NA). An infinite value has
  # no correlation, and would otherwise reach the network as a silently
  # isolated variable.
  infinite <- c(data1 = sum(is.infinite(d1)), data2 = sum(is.infinite(d2)))
  if (any(infinite > 0L))
    .nct_bad_input(sprintf(
      "`data1` and `data2` must hold finite values or NA; found %s.",
      paste(sprintf("%d infinite value(s) in `%s`", infinite[infinite > 0L],
                    names(infinite)[infinite > 0L]), collapse = " and ")))
  n1 <- colnames(d1); n2 <- colnames(d2)
  if (!is.null(n1) && !is.null(n2) && !identical(n1, n2)) {
    listing <- function(v) if (length(v)) paste(v, collapse = ", ") else "none"
    .nct_bad_input(if (setequal(n1, n2) && length(n1) == length(n2))
      "`data1` and `data2` must have the same columns in the same order; put the columns of `data2` in the order of `data1`."
    else sprintf(
      "`data1` and `data2` must have the same columns. Only in data1: %s. Only in data2: %s.",
      listing(setdiff(n1, n2)), listing(setdiff(n2, n1))))
  }
  if (ncol(d1) != ncol(d2))
    .nct_bad_input(sprintf(
      "`data1` has %d numeric variables but `data2` has %d.", ncol(d1), ncol(d2)))
  labels <- n1 %||% n2 %||% paste0("V", seq_len(ncol(d1)))
  colnames(d1) <- colnames(d2) <- labels
  list(data1 = d1, data2 = d2, labels = labels)
}

# Correlation matrix of one (sub)sample, or NULL when a variable pair is never
# co-observed (fewer than 2 shared rows): that correlation is undefined and is
# refused. A variable without variation (over the rows where a pair is observed)
# enters as an isolated node -- its correlations are 0 -- the convention every
# psychnets estimator uses for a constant column (see .psn_refit()). Missing
# values are never filled: with NAs the correlation is pairwise. stats::cor()
# announces a constant variable with a "standard deviation is zero" warning;
# that one warning (matched in the session's own translation) is muffled because
# the case is handled here, and any other warning passes through.
#' @noRd
.nct_cor <- function(x) {
  pairwise <- anyNA(x)
  if (pairwise) {
    co <- crossprod(!is.na(x))
    if (any(co[upper.tri(co)] < 2L)) return(NULL)
  }
  zero_sd <- gettext("the standard deviation is zero", domain = "R-stats")
  S <- withCallingHandlers(
    stats::cor(x, use = if (pairwise) "pairwise.complete.obs" else "everything"),
    warning = function(w) {
      if (identical(conditionMessage(w), zero_sd)) invokeRestart("muffleWarning")
    }
  )
  # Only a variable without variation (over the rows a pair shares) has an
  # undefined correlation here; it is zeroed, the isolated-node convention. Any
  # other non-finite correlation (e.g. overflow of huge values) is an error,
  # never a silently isolated variable.
  undefined <- which(!is.finite(S) & upper.tri(S), arr.ind = TRUE)
  if (nrow(undefined)) {
    flat <- function(v) length(unique(v[!is.na(v)])) < 2L
    explained <- vapply(seq_len(nrow(undefined)), function(k) {
      i <- undefined[k, 1L]; j <- undefined[k, 2L]
      shared <- !is.na(x[, i]) & !is.na(x[, j])
      flat(x[shared, i]) || flat(x[shared, j])
    }, logical(1))
    if (!all(explained)) {
      bad <- undefined[!explained, , drop = FALSE]
      .nct_bad_input(sprintf(
        "The correlation of %s is not finite although both variables vary (are the values too large?).",
        paste(sprintf("%s-%s", colnames(x)[bad[, 1L]], colnames(x)[bad[, 2L]]),
              collapse = ", ")))
    }
    S[undefined] <- 0
    S[undefined[, 2:1, drop = FALSE]] <- 0
  }
  diag(S) <- 1
  S
}

# EBIC-glasso network of one (sub)sample, or NULL when a correlation is
# undefined. Complete data: Pearson correlation -> nearest-PD projection, n =
# rows. Missing data (pairwise): pairwise correlation -> nearest-PD projection,
# n = median pairwise co-observation count, exactly as .cor_input() does for the
# estimators.
#' @noRd
.nct_estimate <- function(x, gamma) {
  S <- .nct_cor(x)
  if (is.null(S)) return(NULL)
  n <- if (anyNA(x)) {
    co <- crossprod(!is.na(x))
    max(round(stats::median(co[upper.tri(co)])), 2L)
  } else nrow(x)
  ebic_glasso(cor_matrix = .nearest_pd_cor(S), n = n, gamma = gamma)$weights
}

# Names of the variables with fewer than two distinct observed values.
#' @noRd
.nct_constant_vars <- function(x) {
  const <- vapply(seq_len(ncol(x)), function(j) {
    v <- x[!is.na(x[, j]), j]
    length(unique(v)) < 2L
  }, logical(1))
  colnames(x)[const]
}

#' @noRd
.nct_group <- function(x, levels, iter, gamma, paired, abs, weighted, p_adjust,
                       na_method, gamma_supplied = FALSE,
                       na_method_supplied = FALSE, paired_supplied = FALSE) {
  if (!identical(attr(x, "source"), "data"))
    .nct_bad_input("net_compare() supports group mode for cross-sectional data only.")
  subs <- attr(x, "subsets")
  pair <- if (is.null(levels)) names(subs) else as.character(levels)
  if (length(pair) != 2L || anyNA(pair) || anyDuplicated(pair))
    .nct_bad_input("Specify two different group levels to compare, e.g. c(1, 4) or c(\"A\", \"B\").")
  miss <- setdiff(pair, names(subs))
  if (length(miss))
    .nct_bad_input(paste0("Group level(s) not found: ", paste(miss, collapse = ", "),
                          ". Available: ", paste(names(subs), collapse = ", "), "."))
  # net_compare() re-estimates each level as an EBIC glasso from (gamma,
  # na_method). The group's own settings are used unless the caller supplies
  # them; say so whenever the networks compared are not the group's own.
  cl <- attr(x, "call")
  method <- .resolve_method(cl$method %||% "glasso")
  dots <- cl$dots %||% list()
  group_gamma <- cl$gamma %||% 0.5          # ebic_glasso()'s default
  group_na <- dots$na_method %||% "pairwise"
  use_gamma <- if (gamma_supplied) gamma else group_gamma
  use_na <- if (na_method_supplied) na_method else group_na
  ignored <- setdiff(names(dots), "na_method")
  if (!isTRUE(all.equal(cl$threshold %||% 0, 0))) ignored <- c(ignored, "threshold")
  overridden <- c(
    if (identical(method, "glasso") && is.numeric(use_gamma) &&
        length(use_gamma) == 1L && !isTRUE(all.equal(use_gamma, group_gamma)))
      sprintf("`gamma = %s` replaces the group's %s", format(use_gamma),
              format(group_gamma)),
    if (!identical(use_na, group_na))
      sprintf("`na_method = \"%s\"` replaces the group's \"%s\"", use_na, group_na))
  if (!identical(method, "glasso") || length(ignored) || length(overridden)) {
    warning(warningCondition(sprintf(
      "net_compare() compares EBIC-glasso networks estimated from `gamma` and `na_method`; %s.",
      paste(c(
        if (!identical(method, "glasso"))
          sprintf("this group's networks were estimated with method = \"%s\"", method),
        if (length(ignored))
          sprintf("the group's %s %s not used",
                  paste(sprintf("`%s`", ignored), collapse = ", "),
                  if (length(ignored) > 1L) "are" else "is"),
        overridden),
        collapse = " and ")),
      class = "psychnet_estimator_mismatch", call = NULL))
  }
  d1 <- subs[[pair[1L]]]
  d2 <- subs[[pair[2L]]]
  # A group built with `within =` knows each row's person: a paired comparison
  # pairs by person, keeping the persons observed in both levels.
  persons <- attr(x, "persons")
  # the same persons at every level: the paired test is the right one, so it
  # is the default for such a group (an explicit `paired = FALSE` still wins)
  if (!paired_supplied && !is.null(persons)) paired <- TRUE
  if (isTRUE(paired) && !is.null(persons)) {
    p1 <- persons[[pair[1L]]]
    p2 <- persons[[pair[2L]]]
    if (anyDuplicated(p1) || anyDuplicated(p2))
      .nct_bad_input(sprintf(paste0(
        "A paired comparison needs one row per person in each level; `%s` ",
        "repeats within level %s."), attr(x, "within"),
        paste(pair[c(anyDuplicated(p1) > 0, anyDuplicated(p2) > 0)], collapse = " and ")))
    # sorted, so the pairs (and hence the permutation draws under a seed) do
    # not depend on the row order of the input
    common <- sort(intersect(p1, p2))
    n_out <- length(p1) + length(p2) - 2L * length(common)
    if (n_out > 0L)
      .psn_warn(sprintf(paste0(
        "%d row(s) of persons observed in only one of levels %s and %s were ",
        "left out of the paired comparison."), n_out, pair[1L], pair[2L]),
        class = "psychnet_rows_dropped")
    d1 <- d1[match(common, p1), , drop = FALSE]
    d2 <- d2[match(common, p2), , drop = FALSE]
  }
  # The group's node labels name the compared variables, as in its networks.
  if (!is.null(cl$labels)) {
    relabel <- function(d) {
      d <- as.data.frame(d, stringsAsFactors = FALSE)
      num <- vapply(d, is.numeric, logical(1))
      if (sum(num) != length(cl$labels))
        .nct_bad_input(sprintf(
          "The group's `labels` name %d variables but its data hold %d numeric variables.",
          length(cl$labels), sum(num)))
      names(d)[num] <- as.character(cl$labels)
      d
    }
    d1 <- relabel(d1)
    d2 <- relabel(d2)
  }
  res <- net_compare(d1, d2, iter = iter, gamma = use_gamma, paired = paired,
                     abs = abs, weighted = weighted, p_adjust = p_adjust,
                     na_method = use_na)
  res[["groups"]] <- pair
  res
}

# Display name of a group: the argument's name when it is a plain variable
# (as t.test() reports "data: x and y"), otherwise "group 1" / "group 2".
#' @noRd
.nct_label <- function(expr, default) {
  if (is.name(expr)) as.character(expr) else default
}

#' Network Comparison Test
#'
#' Permutation test for whether two groups' Gaussian graphical models differ,
#' on three invariants (van Borkulo et al. 2022): network structure (`M`, the
#' largest absolute edge difference), global strength (`S`, the difference in
#' summed absolute edge weights; with `abs = FALSE`, summed signed weights,
#' i.e. global expected influence) and individual edges (`E`). Networks are
#' EBIC graphical lassos (clean-room pure R). Equivalent in purpose to `NetworkComparisonTest::NCT()`, and equal to
#' it in its permutation draws: under the same seed both permute the same rows,
#' so the two differ only by the glasso solver's convergence tolerance.
#'
#' The two groups must hold the same variables in the same order (non-numeric
#' columns such as an id or a group label are dropped first).
#'
#' Missing data follow `na_method`. `"pairwise"` (default) keeps every row,
#' permutes rows as usual, and estimates each (sub)sample from its
#' pairwise-complete correlation matrix (projected to the nearest positive
#' definite correlation matrix), with the median pairwise co-observation count
#' as the sample size -- the same treatment as [ebic_glasso()]. `"listwise"`
#' drops incomplete rows first (in a paired design, a pair is dropped when
#' either member is incomplete). Complete data are unaffected by the choice.
#'
#' A variable without variation in a group, or in a permuted split (a rare
#' binary item can be constant in a split), enters that network as an isolated
#' node -- the convention of every psychnets estimator for a constant column.
#' The statistic is then defined on every split, so the permutation test stays
#' exact. A constant variable in an observed group raises a
#' `psychnet_constant_variable` warning.
#'
#' @param data1,data2 Numeric data frames/matrices with the same variables,
#'   holding finite values or `NA` (an infinite value is a
#'   `psychnet_bad_input` error). When `data1` is a [psychnet_group], `data2`
#'   gives the two levels to compare as they appear in the grouping column,
#'   e.g. `c(1, 4)` for weeks 1 and 4 (default: its two levels when there are
#'   exactly two); the group's node `labels` name the edges, and the group's
#'   `gamma` and `na_method` are used unless you supply `gamma` / `na_method`,
#'   which then win (with a `psychnet_estimator_mismatch` warning when they
#'   differ from the group's).
#' @param iter Number of permutations. Default 1000.
#' @param gamma EBIC hyperparameter. Default 0.5.
#' @param paired Logical; within-row swapping for paired designs (row `i` of
#'   `data1` and of `data2` are the same unit). Default `FALSE`, except for a
#'   group built by `psychnet(..., within = )`: its levels hold the same
#'   persons, so they are compared paired by person unless `paired = FALSE`
#'   (persons seen in only one level are left out, with a
#'   `psychnet_rows_dropped` warning).
#' @param abs Logical; compare absolute edge weights. Default TRUE.
#' @param weighted Logical; if FALSE, binarize networks first. Default TRUE.
#' @param p_adjust Multiple-comparison adjustment for per-edge p-values
#'   (any [stats::p.adjust] method). Default `"none"`.
#' @param na_method `"pairwise"` (default) or `"listwise"` handling of missing
#'   values; see Details.
#' @param x,object A `psychnet_nct` object.
#' @param alpha For the methods: significance level for the `significant`
#'   column, the `*` marks and the edge colouring. Default 0.05.
#' @param n For `print()`: number of edges to list, largest difference first.
#'   Default `Inf` (all).
#' @param type For `plot()`: `"strength"` (default) or `"structure"` draw the
#'   permutation null of S or M with the observed value and p-value marked;
#'   `"edges"` draws the observed per-edge absolute differences, significant
#'   ones coloured and marked `*`.
#' @param row.names,optional Ignored (S3 consistency).
#' @param ... Unused.
#' @return An object of class `psychnet_nct`, a list with `$nw1`, `$nw2` (the
#'   two estimated weighted adjacency matrices); the three tests, each with
#'   `observed`, `perm` and `p_value`: `$structure` (statistic M, the largest
#'   absolute edge difference), `$strength` (statistic S, the absolute
#'   difference in global strength; also `by_group`, each network's global
#'   strength) and `$E` (per-edge absolute differences; also `edge_names`, a
#'   `from`/`to` data frame aligned to the per-edge vector); plus `n_iter`, `n`
#'   (rows used per group), `groups` (display names: the argument names, or
#'   the group levels), `paired` and `params`. The letters follow van
#'   Borkulo et al. (2022) and `NetworkComparisonTest`; `$M` and `$S` (which
#'   before psychnets 0.6.0 held the two global tests the other way round) raise
#'   a `psychnet_renamed` error.
#'
#'   `print()` shows every result: the two global tests, each group's global
#'   strength, the number of differing edges and every edge, largest
#'   difference first. `summary()` returns a data frame with one row per
#'   global test, columns `test` (`"network structure"`; `"global strength"`,
#'   or `"global expected influence"` when `abs = FALSE`), `statistic` (`"M"`,
#'   `"S"`), `observed`, `p_value` and `significant` (`p_value < alpha`).
#'   `as.data.frame()` returns one row per node pair, columns `from`, `to`,
#'   `weight_1`, `weight_2` (the edge in each group's network; 0/1 when
#'   `weighted = FALSE`), `difference` (absolute difference, statistic E),
#'   `p_value` (permutation p, adjusted by `p_adjust`), `p_raw` (the same
#'   permutation p before adjustment; equal to `p_value` when `p_adjust =
#'   "none"`) and `significant` (`p_value < alpha`, i.e. on the adjusted p).
#'   The global M and S tests are never adjusted.
#'   `print()` and `plot()` return `x` invisibly.
#'
#'   Errors of class `psychnet_bad_input` for malformed arguments or groups
#'   whose columns differ (in name or order); of class
#'   `psychnet_degenerate_sample` when (pairwise) a variable pair shares fewer
#'   than two observed rows in a group or a permuted split. Warnings of class
#'   `psychnet_constant_variable` (a constant variable in an observed group) and
#'   `psychnet_estimator_mismatch` (a group object built with an estimator or
#'   settings other than EBIC glasso, or compared with a `gamma` / `na_method`
#'   other than the group's own).
#' @references van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski,
#'   J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022).
#'   Comparing network structures on three aspects: A permutation test.
#'   *Psychological Methods*. \doi{10.1037/met0000476}
#' @examples
#' set.seed(1)
#' # group a has a stronger chain than group b
#' a <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.5^abs(outer(1:4, 1:4, "-")))
#' b <- matrix(stats::rnorm(100 * 4), 100, 4) %*% chol(0.2^abs(outer(1:4, 1:4, "-")))
#' colnames(a) <- colnames(b) <- paste0("V", 1:4)
#' # iter = 10 keeps the example fast (smallest p = 1/11); use iter >= 1000.
#' fit <- net_compare(a, b, iter = 10)
#' fit
#' print(fit, n = 3)          # the three largest edge differences
#' summary(fit)
#' as.data.frame(fit)
#' plot(fit)                  # global strength permutation null
#' plot(fit, type = "edges")  # per-edge differences
#' @export
net_compare <- function(data1, data2 = NULL, iter = 1000L, gamma = 0.5,
                        paired = FALSE, abs = TRUE, weighted = TRUE,
                        p_adjust = "none", na_method = c("pairwise", "listwise")) {
  na_method_supplied <- !missing(na_method)
  na_method <- .psn_match_arg(na_method)
  groups <- c(.nct_label(substitute(data1), "group 1"),
              .nct_label(substitute(data2), "group 2"))
  if (inherits(data1, "psychnet_group")) {
    return(.nct_group(data1, data2, iter = iter, gamma = gamma, paired = paired,
                      abs = abs, weighted = weighted, p_adjust = p_adjust,
                      na_method = na_method, gamma_supplied = !missing(gamma),
                      na_method_supplied = na_method_supplied,
                      paired_supplied = !missing(paired)))
  }
  if (is.null(data2))
    .nct_bad_input(paste0("`data2` is required (a second data set, or two ",
                          "group levels when `data1` is a psychnet_group)."))
  .nct_validate(iter, gamma, paired, abs, weighted)
  p_adjust <- .psn_match_arg(p_adjust, stats::p.adjust.methods)
  cols <- .nct_columns(data1, data2)
  data1 <- cols$data1
  data2 <- cols$data2
  labs <- cols$labels

  if (paired && nrow(data1) != nrow(data2))
    .nct_bad_input("A paired comparison requires the two groups to have equal size.")
  if (na_method == "listwise") {
    keep1 <- stats::complete.cases(data1)
    keep2 <- stats::complete.cases(data2)
    if (paired) keep1 <- keep2 <- keep1 & keep2
    data1 <- data1[keep1, , drop = FALSE]
    data2 <- data2[keep2, , drop = FALSE]
  }
  iter <- as.integer(round(iter))
  n1 <- nrow(data1); n2 <- nrow(data2)
  if (n1 < 3L || n2 < 3L)
    .nct_bad_input(sprintf("Each group needs at least 3 rows; got %d and %d.",
                           n1, n2))
  dataall <- rbind(data1, data2)

  undefined <- function(where) .nct_degenerate(sprintf(
    "A variable pair shares fewer than 2 observed rows in %s, so its correlation is undefined. Use na_method = \"listwise\" or drop the sparsely observed variable.",
    where))
  est <- function(x, where) {
    nw <- .nct_estimate(x, gamma)
    if (is.null(nw)) undefined(where)
    nw
  }
  binarize <- function(m) (m != 0) * 1

  constant <- list(data1 = .nct_constant_vars(data1),
                   data2 = .nct_constant_vars(data2))
  constant <- constant[lengths(constant) > 0L]
  if (length(constant))
    warning(warningCondition(sprintf(
      "%s; a constant variable enters that network as an isolated node.",
      paste(sprintf("constant in `%s`: %s", names(constant),
                    vapply(constant, paste, character(1), collapse = ", ")),
            collapse = "; ")),
      class = "psychnet_constant_variable", call = NULL))

  nw1 <- est(data1, "`data1`"); nw2 <- est(data2, "`data2`")
  ut <- upper.tri(nw1)
  # Edge labels for the per-edge `E` vector, in the same upper-triangle order as
  # E$observed (column-major over the TRUE cells of `ut`).
  ij <- which(ut, arr.ind = TRUE)
  edge_names <- data.frame(from = labs[ij[, 1L]], to = labs[ij[, 2L]],
                           stringsAsFactors = FALSE, row.names = NULL)
  if (!weighted) { nw1 <- binarize(nw1); nw2 <- binarize(nw2) }

  # Notation of van Borkulo et al. (2022) and NetworkComparisonTest: M is the
  # network-structure invariance statistic (largest absolute edge difference),
  # S the global-strength invariance statistic (difference in summed edge
  # weights; summed absolute weights when `abs`, i.e. global strength, else
  # signed weights, i.e. global expected influence), E the per-edge statistic.
  global <- function(m) if (abs) sum(base::abs(m[ut])) else sum(m[ut])
  strength_by_group <- c(global(nw1), global(nw2))
  str_obs <- base::abs(strength_by_group[1L] - strength_by_group[2L])
  diff_real <- base::abs(nw1 - nw2)
  max_obs <- max(diff_real[ut])
  E_obs <- diff_real[ut]
  n_edges <- length(E_obs)

  # Permutation null. A literal for loop is the right tool here: vectorising
  # would materialise iter x n_edges intermediate networks.
  str_perm <- numeric(iter); max_perm <- numeric(iter)
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
    where <- sprintf("permuted split %d", i)
    r1 <- est(x1p, where); r2 <- est(x2p, where)
    if (!weighted) { r1 <- binarize(r1); r2 <- binarize(r2) }
    str_perm[i] <- base::abs(global(r1) - global(r2))
    diff_perm <- base::abs(r1 - r2)
    max_perm[i] <- max(diff_perm[ut])
    E_perm[i, ] <- diff_perm[ut]
  }

  str_pval <- (sum(str_perm >= str_obs) + 1) / (iter + 1)
  max_pval <- (sum(max_perm >= max_obs) + 1) / (iter + 1)
  E_pval <- (colSums(E_perm >= matrix(E_obs, iter, n_edges, byrow = TRUE)) + 1) /
            (iter + 1)
  # Keep the unadjusted per-edge p-values next to the adjusted ones, so a
  # corrected analysis still reports both (M and S are single tests and are
  # never adjusted, as in NetworkComparisonTest).
  E_pval_raw <- E_pval
  if (p_adjust != "none") E_pval <- stats::p.adjust(E_pval, method = p_adjust)

  structure(
    list(nw1 = nw1, nw2 = nw2,
         structure = list(observed = max_obs, perm = max_perm, p_value = max_pval),
         strength = list(observed = str_obs, perm = str_perm, p_value = str_pval,
                         by_group = strength_by_group),
         E = list(observed = E_obs, perm = E_perm, p_value = E_pval,
                  p_raw = E_pval_raw, edge_names = edge_names),
         groups = groups, n_iter = iter, n = c(n1, n2), paired = paired,
         params = list(gamma = gamma, abs = abs, weighted = weighted,
                       p_adjust = p_adjust, na_method = na_method)),
    class = "psychnet_nct"
  )
}

#' @describeIn net_compare The M and S tests as a table.
#' @export
summary.psychnet_nct <- function(object, alpha = 0.05, ...) {
  stopifnot("`alpha` must be a single number in (0, 1)" =
              is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)
  p <- c(object$structure$p_value, object$strength$p_value)
  data.frame(
    test = c("network structure",
             if (isFALSE(object$params$abs)) "global expected influence"
             else "global strength"),
    statistic = c("M", "S"),
    observed = c(object$structure$observed, object$strength$observed),
    p_value = p, significant = p < alpha,
    stringsAsFactors = FALSE)
}

#' @describeIn net_compare The per-edge tests, one row per node pair.
#' @export
as.data.frame.psychnet_nct <- function(x, row.names = NULL, optional = FALSE,
                                       alpha = 0.05, ...) {
  stopifnot("`alpha` must be a single number in (0, 1)" =
              is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1)
  ut <- upper.tri(x$nw1)
  data.frame(from = x$E$edge_names$from, to = x$E$edge_names$to,
             weight_1 = x$nw1[ut], weight_2 = x$nw2[ut],
             difference = x$E$observed, p_value = x$E$p_value,
             p_raw = x$E$p_raw %||% x$E$p_value,
             significant = x$E$p_value < alpha,
             stringsAsFactors = FALSE, row.names = NULL)
}

#' @describeIn net_compare Every result, edges largest difference first.
#' @export
print.psychnet_nct <- function(x, n = Inf, alpha = 0.05, ...) {
  cat(.nct_display(x, n = n, alpha = alpha), sep = "\n")
  invisible(x)
}

# The printed report, as lines. Fixed-width columns; numbers to 3 decimals;
# p-values below 0.001 shown as "<0.001"; "*" marks p < alpha.
#' @noRd
.nct_display <- function(x, n = Inf, alpha = 0.05) {
  g <- x$groups %||% c("group 1", "group 2")
  gl <- summary(x, alpha = alpha)
  ed <- as.data.frame(x, alpha = alpha)
  num <- function(v) formatC(v, format = "f", digits = 3)
  pval <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
  star <- function(p) ifelse(p < alpha, "*", "")
  pad <- function(v, w, right = TRUE) formatC(v, width = if (right) w else -w)
  adj <- x$params$p_adjust %||% "none"

  header <- c(
    sprintf("Network Comparison Test: %s vs %s", g[1L], g[2L]),
    sprintf("%d permutations, %s groups, n = %d and %d",
            x$n_iter, if (isTRUE(x$paired)) "paired" else "independent",
            x$n[1L], x$n[2L]),
    "")
  n_sig <- sum(ed$significant)
  labels <- c("M  network structure", paste("S ", gl$test[2L]), "E  edges differing")
  wl <- max(nchar(labels))
  test_lines <- c(
    paste0("  ", pad("Test", wl, FALSE), "  ", pad("Statistic", 9), "  ", pad("p", 6)),
    paste0("  ", pad(labels[1:2], wl, FALSE), "  ", pad(num(gl$observed), 9), "  ",
           pad(pval(gl$p_value), 6), " ", star(gl$p_value)),
    sprintf("       %s %s vs %s %s", g[1L], num(x$strength$by_group[1L]),
            g[2L], num(x$strength$by_group[2L])),
    paste0("  ", pad(labels[3L], wl, FALSE), "  ",
           pad(sprintf("%d of %d", n_sig, nrow(ed)), 9)))

  k <- min(nrow(ed), n)
  edge_lines <- character(0)
  if (k > 0) {
    top <- ed[order(-ed$difference, ed$from, ed$to)[seq_len(k)], ]
    edge <- paste(top$from, "-", top$to)
    gn <- substr(g, 1L, 12L)
    we <- max(nchar(c("Edge", edge)))
    wg <- pmax(nchar(gn), 6L)
    edge_lines <- c(
      "",
      if (k < nrow(ed)) sprintf("  Largest edge differences (%d of %d)", k, nrow(ed))
      else "  Edges, largest difference first",
      paste0("  ", pad("Edge", we, FALSE), "  ", pad(gn[1L], wg[1L]), "  ",
             pad(gn[2L], wg[2L]), "  ", pad("|diff|", 6), "  ", pad("p", 6)),
      paste0("  ", pad(edge, we, FALSE), "  ", pad(num(top$weight_1), wg[1L]),
             "  ", pad(num(top$weight_2), wg[2L]), "  ",
             pad(num(top$difference), 6), "  ", pad(pval(top$p_value), 6), " ",
             star(top$p_value)))
  }
  footer <- c("", sprintf("  * p < %g%s", alpha,
                          if (identical(adj, "none")) ""
                          else sprintf("; edge p-values %s-adjusted", adj)))
  c(header, test_lines, edge_lines, footer)
}

#' @noRd
#' @export
`$.psychnet_nct` <- function(x, name) {
  if (name %in% c("M", "S")) {
    stop(errorCondition(sprintf(paste0(
      "`$%s` was removed in psychnets 0.6.0: it held the %s test, which ",
      "van Borkulo et al. (2022) and NetworkComparisonTest call %s. Use ",
      "summary(x), or `$structure` (M) / `$strength` (S)."),
      name, if (name == "M") "global-strength" else "network-structure",
      if (name == "M") "S" else "M"),
      class = c("psychnet_renamed", "psychnet_error"), call = NULL))
  }
  .subset2(x, name)
}
