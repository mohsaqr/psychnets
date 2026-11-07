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

#' Correctness certificate of a fitted network
#'
#' Regularized or constrained estimators in `psychnets` self-certify: each
#' fitted network carries a residual that reports how far it sits from the
#' optimum of its own objective (a KKT / stationarity residual), or -- for the
#' structural methods (TMFG, relative importance) -- whether the graph
#' satisfies the identity that defines it. This verb returns that certificate
#' as a tidy one-row `data.frame`, so correctness is read the same way for
#' every method.
#'
#' The residual is near machine zero for a correctly solved problem. A graph
#' altered after fitting (post-estimation `threshold`, Wald pruning in
#' [ising_sampler()]) or an external-engine fit without an available residual
#' reports `NA` rather than being presented as certified. `cor` and `pcor` have
#' no optimization to certify and report `NA` with kind `"none"`.
#'
#' @param x A [psychnet] object.
#' @param tol Tolerance below which the fit is flagged `certified = TRUE`: a
#'   single finite number >= 0. Default `1e-6`. Anything else is an error of
#'   class `psychnet_bad_input`.
#' @return A one-row `data.frame` with columns `method`, `certificate` (the
#'   residual; smaller is better), `kind` (`"kkt"` for the optimization
#'   certificates, `"structural"` for TMFG/relimp, `"none"` for cor/pcor), and
#'   `certified` (logical: residual at or below `tol`; `NA` when no applicable
#'   certificate is available). Correlation networks use kind `"none"` and are
#'   reported as `TRUE` because no optimization claim is made.
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
