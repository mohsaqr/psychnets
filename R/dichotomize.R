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
#'     \item{`"rank"`}{`1` for the upper half of the column by rank, giving a
#'       balanced (~50/50) split that is robust to ties (useful for coarse
#'       Likert items where a median split is badly unbalanced).}
#'   }
#' @return An integer matrix of `0`/`1` values with the same dimensions and
#'   dimnames as `data`.
#' @examples
#' b <- dichotomize(SRL_GPT, method = "median")
#' table(b)                          # values are 0/1 only
#' @export
dichotomize <- function(data, method = c("median", "mean", "rank")) {
  method <- match.arg(method)
  X <- as.matrix(data)
  stopifnot(is.numeric(X), ncol(X) >= 1L)
  split_col <- switch(method,
    median = function(x) as.integer(x >= stats::median(x)),
    mean   = function(x) as.integer(x > mean(x)),
    rank   = function(x) as.integer(rank(x, ties.method = "first") > length(x) / 2))
  b <- vapply(seq_len(ncol(X)), function(j) split_col(X[, j]), integer(nrow(X)))
  dimnames(b) <- dimnames(X)
  b
}
