# Node predictability (Haslbeck & Waldorp 2018), clean-room base R. How well
# each node is predicted by its neighbours in the fitted network: R-squared for
# continuous (Gaussian) nodes, and classification accuracy (CC) plus the
# normalized accuracy (nCC, accuracy above the marginal baseline) for binary
# nodes. Equivalent in purpose to mgm::predictability().
#
# For a Gaussian graphical model the value is closed-form from the precision
# matrix and needs no raw data: the residual variance of node j given the rest
# is 1 / Theta_jj, so R^2_j = 1 - 1 / (Theta_jj * S_jj). For the nodewise models
# (ising, IsingSampler, mgm) predictability is computed from the data using the
# stored nodewise coefficients.

#' Node predictability
#'
#' Reports how well each node is predicted by the others in a fitted network.
#' For Gaussian graphical models this is the closed-form variance explained
#' (R-squared) from the precision matrix and needs no data. For the nodewise
#' models ([ising_fit()], [ising_sampler()], [mgm_fit()]) it requires the data
#' and reports R-squared for Gaussian nodes and classification accuracy (`CC`)
#' plus normalized accuracy (`nCC`) for binary nodes.
#'
#' @param x A [psychnet] object.
#' @param data The data the network was estimated from; required for the
#'   nodewise models (ising / IsingSampler / mgm), ignored for the GGMs.
#' @param ... Unused.
#' @return A tidy `data.frame`, one row per node, with columns `node`, `type`
#'   (`"gaussian"` or `"binary"`), `metric` (`"R2"` or `"nCC"`),
#'   `predictability`, and `accuracy` (classification accuracy for binary nodes,
#'   `NA` for Gaussian).
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' predictability(ebic_glasso(cor_matrix = S, n = 250))
#' @export
predictability <- function(x, data = NULL, ...) {
  stopifnot(inherits(x, "psychnet"))
  labs <- x$nodes$label
  p <- nrow(x$nodes)

  # --- nodewise models: compute from data ------------------------------------
  if (!is.null(x$nodewise)) {
    if (is.null(data)) {
      stop(sprintf("`data` is required for predictability of a '%s' network.",
                   x$method), call. = FALSE)
    }
    mat <- .as_numeric_matrix(data)
    if (!is.null(colnames(mat))) {
      # named data must contain every node, selected by name (not by position)
      if (!all(labs %in% colnames(mat))) {
        stop("`data` is missing columns for some network nodes.", call. = FALSE)
      }
      mat <- mat[, labs, drop = FALSE]
    }
    if (ncol(mat) != p) {
      stop("`data` columns do not match the network nodes.", call. = FALSE)
    }
    nw <- x$nodewise
    Xstd <- sweep(sweep(mat, 2L, nw$center, "-"), 2L, nw$scale, "/")
    rows <- lapply(seq_len(p), function(i) {
      eta <- nw$intercept[i] + as.numeric(Xstd %*% nw$beta_std[i, ])
      y <- mat[, i]
      if (nw$families[i] == "binomial") {
        cls  <- as.numeric(1 / (1 + exp(-eta)) > 0.5)
        cc   <- mean(cls == y)
        marg <- max(mean(y), 1 - mean(y))
        ncc  <- if (marg < 1) (cc - marg) / (1 - marg) else 0
        data.frame(node = labs[i], type = "binary", metric = "nCC",
                   predictability = ncc, accuracy = cc,
                   stringsAsFactors = FALSE)
      } else {
        # mgm scales the gaussian response, so eta is on the standardized scale;
        # put y there too (R-squared is scale-invariant, so this is a no-op for
        # estimators that did not scale the response).
        ys <- if (!is.null(nw$resp_scale)) {
          (y - nw$resp_center[i]) / nw$resp_scale[i]
        } else y
        r2 <- 1 - sum((ys - eta)^2) / sum((ys - mean(ys))^2)
        data.frame(node = labs[i], type = "gaussian", metric = "R2",
                   predictability = r2, accuracy = NA_real_,
                   stringsAsFactors = FALSE)
      }
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  # --- Gaussian graphical model: closed form from the precision --------------
  if (is.null(x$precision)) {
    stop(sprintf("predictability() is not defined for a '%s' network.",
                 x$method), call. = FALSE)
  }
  theta <- x$precision
  s_diag <- if (!is.null(x$cor_matrix)) diag(x$cor_matrix) else rep(1, p)
  r2 <- 1 - 1 / (diag(theta) * s_diag)
  data.frame(node = labs, type = "gaussian", metric = "R2",
             predictability = as.numeric(r2), accuracy = NA_real_,
             row.names = NULL, stringsAsFactors = FALSE)
}
