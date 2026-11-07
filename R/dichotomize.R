# Dichotomize numeric/ordinal columns to 0/1 -- the standard preprocessing step
# before fitting an Ising network to Likert data.

#' Dichotomize numeric columns to 0/1
#'
#' Splits each column of a numeric matrix or data frame into a binary 0/1
#' variable. This is the usual preprocessing step before fitting an Ising
#' network ([ising_fit()], [ising_sampler()]) to Likert or other
#' ordinal/continuous data, which require binary input.
#'
#' @param data Numeric matrix or data frame (rows = observations).
#' @param method Split rule, applied independently to each column:
#'   \describe{
#'     \item{`"median"`}{(default) `1` if the value is `>=` the column median.}
#'     \item{`"mean"`}{`1` if the value is `>` the column mean.}
#'     \item{`"rank"`}{`1` for observations whose average rank is above the
#'       midpoint. Tied values always remain together, so the split may be
#'       unbalanced rather than depending on row order.}
#'   }
#'   On a floor- or ceiling-heavy item the median can equal the column minimum,
#'   so `"median"` (`>=`) puts every observation in the upper group and the
#'   output column is constant. Such columns are reported by a
#'   `psychnet_constant_variable` warning (an Ising fit would drop them);
#'   `method = "mean"` or `"rank"` splits them.
#' @return An integer matrix of `0`/`1` values with the same dimensions and
#'   dimnames as `data`. Missing values stay `NA` in the output (they are not
#'   assigned to either group). A non-numeric `data`, an empty one, or a
#'   constant / all-missing input column is an error of class
#'   `psychnet_bad_input`.
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
