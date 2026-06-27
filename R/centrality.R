# Node centrality for psychometric networks (pure base R). Reports the two
# measures recommended for psychometric networks -- strength and expected
# influence (Robinaugh, Millner & McNally 2016) -- which, unlike closeness and
# betweenness, are well-defined on signed, weighted graphs.

#' Node centrality
#'
#' @param x A [psychnet] object or a weighted adjacency matrix.
#' @param ... Unused.
#' @return A tidy `data.frame`, one row per node, with columns `node`,
#'   `strength` (sum of absolute edge weights) and `expected_influence` (sum of
#'   signed edge weights).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' centrality(ebic_glasso(cor_matrix = S, n = 250))
#' @export
centrality <- function(x, ...) {
  if (inherits(x, "psychnet")) {
    g <- x$weights
    labs <- x$nodes$label
  } else {
    if (!is.matrix(x) && !is.data.frame(x)) {
      stop("`x` must be a psychnet object or a square weighted matrix.",
           call. = FALSE)
    }
    g <- as.matrix(x)
    if (!is.numeric(g) || nrow(g) != ncol(g)) {
      stop("`x` must be a square numeric weighted adjacency matrix.",
           call. = FALSE)
    }
    labs <- colnames(g)
    if (is.null(labs)) labs <- paste0("V", seq_len(ncol(g)))
  }
  diag(g) <- 0
  data.frame(
    node               = labs,
    strength           = rowSums(abs(g)),
    expected_influence = rowSums(g),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
