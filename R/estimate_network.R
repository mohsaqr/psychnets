# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...).

#' Estimate a psychometric network
#'
#' Single entry point that routes to the requested estimator and returns a
#' common [psychnet] object, so callers can swap estimators without rewiring
#' downstream code. Mirrors `bootnet::estimateNetwork(data, default = ...)`.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator: `"EBICglasso"`, `"cor"`, `"pcor"`, `"ising"`,
#'   `"mgm"`, `"huge"`, `"ggmModSelect"`, `"TMFG"`, `"LoGo"`, `"relimp"`, or
#'   `"IsingSampler"`. Common aliases are accepted (e.g. `"glasso"` ->
#'   `"EBICglasso"`, `"IsingFit"` -> `"ising"`).
#' @param threshold Absolute-weight threshold below which edges are zeroed.
#' @param gamma EBIC hyperparameter for the regularized methods.
#' @param labels Optional node labels.
#' @param ... Passed to the underlying estimator (e.g. `npn=` for `"huge"`,
#'   `alpha=` for the correlation / `"IsingSampler"` methods).
#' @return A `psychnet` object.
#' @examples
#' x <- matrix(stats::rnorm(200 * 5), 200, 5)
#' estimate_network(x, method = "EBICglasso")
#' estimate_network(x, method = "pcor")
#' @export
estimate_network <- function(data,
                             method = c("EBICglasso", "cor", "pcor",
                                        "ising", "mgm", "huge", "ggmModSelect",
                                        "TMFG", "LoGo", "relimp",
                                        "IsingSampler"),
                             threshold = 0, gamma = 0.5, labels = NULL, ...) {
  method <- .resolve_method(method)
  switch(
    method,
    cor          = cor_network(data, threshold = threshold, labels = labels, ...),
    pcor         = pcor_network(data, threshold = threshold, labels = labels, ...),
    EBICglasso   = ebic_glasso(data, gamma = gamma, threshold = threshold,
                               labels = labels, ...),
    ising        = ising_fit(data, gamma = gamma, labels = labels, ...),
    mgm          = mgm_fit(data, gamma = gamma, labels = labels, ...),
    huge         = huge_network(data, gamma = gamma, threshold = threshold,
                                labels = labels, ...),
    ggmModSelect = ggm_modselect(data, gamma = gamma, threshold = threshold,
                                 labels = labels, ...),
    TMFG         = tmfg_network(data, labels = labels, ...),
    LoGo         = logo_network(data, threshold = threshold, labels = labels, ...),
    relimp       = relimp_network(data, labels = labels, ...),
    IsingSampler = ising_sampler(data, labels = labels, ...),
    stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
}

# Resolve method name + aliases.
#' @noRd
.resolve_method <- function(method) {
  if (length(method) > 1L) method <- method[1L]
  aliases <- c(glasso = "EBICglasso", ebicglasso = "EBICglasso",
               EBICglasso = "EBICglasso", isingfit = "ising", IsingFit = "ising",
               ising = "ising", cor = "cor", correlation = "cor",
               pcor = "pcor", partial = "pcor", mgm = "mgm",
               huge = "huge", nonparanormal = "huge", npn = "huge",
               ggmmodselect = "ggmModSelect", ggmModSelect = "ggmModSelect",
               modselect = "ggmModSelect", stepwise = "ggmModSelect",
               tmfg = "TMFG", TMFG = "TMFG",
               logo = "LoGo", LoGo = "LoGo",
               relimp = "relimp", relativeimportance = "relimp",
               isingsampler = "IsingSampler", IsingSampler = "IsingSampler")
  key <- aliases[method]
  if (is.na(key)) key <- aliases[tolower(method)]
  if (is.na(key)) stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  unname(key)
}
