# Dichotomize numeric/ordinal columns to 0/1 -- the standard preprocessing step
# before fitting an Ising network to Likert data.

#' Split numeric items into binary 0/1 variables
#'
#' @description
#' `dichotomize()` splits numeric items into binary 0/1 variables, each
#' column independently at its median (default), at its mean, or by rank.
#' Ising network models ([ising_fit()], [ising_sampler()]) describe binary
#' variables, so Likert-type or continuous items must be split before such a
#' model can be fitted; the resulting network then describes associations
#' between scoring high and scoring low on the items rather than between
#' their full response scales.
#'
#' Missing values stay missing. The result is an integer 0/1 matrix with the
#' dimensions of the input.
#'
#' Splitting discards the information in the response scale. To keep it,
#' estimate a Gaussian graphical model, for example [ebic_glasso()] with
#' `cor_method = "auto"` (polychoric correlations for ordinal items), or a
#' mixed graphical model with [mgm_fit()].
#'
#' @param data A numeric matrix or data frame with one row per observation.
#' @param method Split rule, applied to each column:
#'   \describe{
#'     \item{`"median"`}{(default) `1` when the value is at or above the
#'       column median.}
#'     \item{`"mean"`}{`1` when the value is above the column mean.}
#'     \item{`"rank"`}{`1` when the value's average rank is above the
#'       midpoint. Tied values stay together, so the split can be
#'       unbalanced but never depends on row order.}
#'   }
#'   On an item with many responses at its floor, the median can equal the
#'   minimum, so `"median"` puts every observation in the upper group and the
#'   output column is constant (an Ising fit would drop it); `"mean"` or
#'   `"rank"` splits such items.
#' @return An integer matrix of 0/1 values with the dimensions and dimnames
#'   of `as.matrix(data)`; missing values stay `NA`.
#' @section Conditions:
#'   `psychnet_bad_input` when `data` is not numeric or has no columns, a
#'   column is constant or entirely missing, or `method` is not one of the
#'   rules; `psychnet_constant_variable` (warning) when a split puts every
#'   observation of a column in one group.
#' @seealso [ising_fit()] and [ising_sampler()], which need binary data.
#' @examples
#' b <- dichotomize(SRL_GPT, method = "median")
#' table(b)                          # values are 0/1 only
#' @export
dichotomize <- function(data, method = c("median", "mean", "rank")) {
  method <- .psn_match_arg(method)
  X <- as.matrix(data)
  if (!is.numeric(X) || ncol(X) < 1L)
    .psn_stop("`data` must be a numeric matrix or data frame with at least one column.")
  nms <- colnames(X); if (is.null(nms)) nms <- paste0("V", seq_len(ncol(X)))
  n_distinct <- function(x) length(unique(x[!is.na(x)]))
  constant <- vapply(seq_len(ncol(X)), function(j) n_distinct(X[, j]) < 2L,
                     logical(1))
  if (any(constant)) {
    .psn_stop("Cannot dichotomize constant or all-missing column(s): ",
         paste(nms[constant], collapse = ", "), ".", call. = FALSE)
  }
  split_col <- switch(method,
    median = function(x) as.integer(x >= stats::median(x, na.rm = TRUE)),
    mean   = function(x) as.integer(x > mean(x, na.rm = TRUE)),
    rank   = function(x) {
      n_obs <- sum(!is.na(x))
      as.integer(rank(x, ties.method = "average", na.last = "keep") >
                   (n_obs + 1) / 2)
    })
  b <- vapply(seq_len(ncol(X)), function(j) split_col(X[, j]), integer(nrow(X)))
  dim(b) <- dim(X)
  dimnames(b) <- dimnames(X)
  one_group <- vapply(seq_len(ncol(b)), function(j) n_distinct(b[, j]) < 2L,
                      logical(1))
  if (any(one_group)) {
    .psn_warn(sprintf(paste0(
      "method = \"%s\" put every observation of %d column(s) in one group, so ",
      "the output is constant: %s. A constant column carries no information ",
      "for an Ising network; try method = \"mean\" or \"rank\"."),
      method, sum(one_group), paste(nms[one_group], collapse = ", ")),
      class = "psychnet_constant_variable")
  }
  b
}
