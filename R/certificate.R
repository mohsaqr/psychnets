# A single tidy accessor for every estimator's correctness certificate, so the
# user never reaches into `$kkt` or calls a method-specific certificate function
# by hand. Each estimator already stores its residual; this verb reads it and
# reports it uniformly.

# What a method's certificate measures: "none" (cor/pcor: no optimization),
# "structural" (TMFG/relimp: a defining identity, not an optimum) or "kkt"
# (every regularized / constrained estimator). Shared with print.psychnet so
# the printed label and certificate()'s `kind` can never disagree.
#' @noRd
.certificate_kind <- function(method) {
  switch(method,
    cor = "none", pcor = "none",
    tmfg = "structural", relimp = "structural",
    "kkt")
}

#' Check that a fitted network solves its estimator
#'
#' @description
#' `certificate()` measures how far a fitted network is from the solution of
#' its estimator. For a regularized or constrained estimator such as the
#' graphical lasso, this is the largest violation of the optimality
#' (Karush–Kuhn–Tucker) conditions; a value near zero means the network is
#' the optimum, up to numerical precision.
#'
#' Such an estimator defines the network as the solution of an optimization
#' problem; the graphical lasso, for example, maximizes a penalized
#' likelihood. Any solution must satisfy the problem's optimality conditions,
#' so whether a returned network is that solution can be checked without a
#' second solver; the amount by which the network violates the conditions is
#' its certificate. The information-filtering (TMFG) and relative-importance
#' networks are not optima of an objective; for them the certificate checks
#' the identity that defines the graph. A network changed after fitting, by a
#' threshold or by pruning, is no longer the solution that was checked, and
#' its certificate is reported as missing.
#'
#' `certificate()` returns the value together with the kind of certificate
#' and whether it falls below a tolerance.
#'
#' Correlation and partial-correlation networks involve no optimization and
#' have no certificate. To check a precision matrix or regression
#' coefficients that are not stored in a fitted network, call
#' [glasso_kkt()], [ggm_support_kkt()] or [glm_lasso_kkt()] directly.
#'
#' @details
#' The optimality certificates are computed by [glasso_kkt()] (graphical lasso
#' networks), [ggm_support_kkt()] (unregularized networks restricted to a
#' graph) and [glm_lasso_kkt()] (nodewise Ising and mixed graphical models).
#' The structural certificates are [tmfg_certificate()], the number of
#' violated TMFG properties (0 to 4), and [lmg_certificate()], the largest
#' difference between the summed importance shares of a node's predictors and
#' that node's explained variance.
#'
#' @param x A fitted network, an object of class `psychnet` ([psychnet-class]).
#' @param tol Tolerance at or below which the network is reported as
#'   certified: a single finite number >= 0. Default `1e-6`.
#' @return A data frame with one row and columns `method` (the estimator),
#'   `certificate` (the violation; smaller is better), `kind` (`"kkt"` for an
#'   optimality certificate, `"structural"` for TMFG and relative importance,
#'   `"none"` for correlation and partial-correlation networks) and
#'   `certified` (`TRUE` when `certificate` is at or below `tol`). For a
#'   network changed after fitting, or a fit that did not converge,
#'   `certificate` and `certified` are `NA`. For correlation and
#'   partial-correlation networks `certificate` is `NA` and `certified` is
#'   `TRUE`, because no optimization claim is made.
#' @section Conditions:
#'   `psychnet_bad_input` when `x` is not a `psychnet` network or `tol` is
#'   not a single finite number >= 0.
#' @seealso [glasso_kkt()], [ggm_support_kkt()], [glm_lasso_kkt()],
#'   [tmfg_certificate()] and [lmg_certificate()] for the individual checks;
#'   [psychnet()] for the estimators.
#' @examples
#' S <- 0.4^abs(outer(1:6, 1:6, "-"))
#' certificate(ebic_glasso(cor_matrix = S, n = 250))
#' certificate(tmfg_network(cor_matrix = S))
#' @export
certificate <- function(x, tol = 1e-6) {
  if (!inherits(x, "psychnet"))
    .psn_stop("`x` must be a psychnet network object.")
  .psn_check_number(tol, "tol", lower = 0)
  kind <- .certificate_kind(x$method)
  value <- switch(x$method,
    cor = NA_real_, pcor = NA_real_,
    tmfg = tmfg_certificate(x),
    if (is.null(x$kkt)) NA_real_ else x$kkt)
  data.frame(
    method      = x$method,
    certificate = value,
    kind        = kind,
    certified   = if (kind == "none") TRUE else if (is.na(value)) NA else value <= tol,
    row.names   = NULL, stringsAsFactors = FALSE
  )
}
