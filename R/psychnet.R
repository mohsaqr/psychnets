# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...).

#' Estimate a psychometric network
#'
#' The package's main entry point: routes to the requested estimator and returns
#' a common `psychnet` object, so callers can swap estimators without rewiring
#' downstream code. Mirrors `bootnet::estimateNetwork(data, default = ...)`.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator: `"EBICglasso"`, `"cor"`, `"pcor"`, `"ising"`,
#'   `"mgm"`, `"huge"`, `"ggmModSelect"`, `"TMFG"`, `"LoGo"`, `"relimp"`, or
#'   `"IsingSampler"`. Common aliases are accepted (e.g. `"glasso"` ->
#'   `"EBICglasso"`, `"IsingFit"` -> `"ising"`).
#' @param threshold Absolute-weight threshold below which edges are zeroed
#'   (forwarded only to the methods that take it: `cor`, `pcor`, `EBICglasso`,
#'   `huge`, `ggmModSelect`, `LoGo`).
#' @param gamma EBIC hyperparameter. `NULL` (default) keeps each method's own
#'   default (0.5 for the Gaussian graphical models, 0.25 for `ising`/`mgm`);
#'   set it to override. Forwarded only to the regularized methods.
#' @param labels Optional node labels.
#' @param ... Passed to the underlying estimator (e.g. `cor_method=` for the
#'   correlation-based methods, `npn=` for `"huge"`, `rule=` for the Ising
#'   methods, `alpha=` for the correlation / `"IsingSampler"` methods).
#' @return A `psychnet` object.
#' @examples
#' x <- matrix(stats::rnorm(200 * 5), 200, 5)
#' psychnet(x, method = "EBICglasso")
#' psychnet(x, method = "pcor", cor_method = "spearman")
#' @export
psychnet <- function(data,
                             method = c("EBICglasso", "cor", "pcor",
                                        "ising", "mgm", "huge", "ggmModSelect",
                                        "TMFG", "LoGo", "relimp",
                                        "IsingSampler"),
                             threshold = 0, gamma = NULL, labels = NULL, ...) {
  method <- .resolve_method(method)
  fn <- switch(
    method,
    cor = cor_network, pcor = pcor_network, EBICglasso = ebic_glasso,
    ising = ising_fit, mgm = mgm_fit, huge = huge_network,
    ggmModSelect = ggm_modselect, TMFG = tmfg_network, LoGo = logo_network,
    relimp = relimp_network, IsingSampler = ising_sampler,
    stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
  # Forward each shared argument only to the methods that accept it, so the
  # callee's own defaults stand otherwise (notably gamma: ising/mgm = 0.25).
  args <- c(list(data = data, labels = labels), list(...))
  if (method %in% c("cor", "pcor", "EBICglasso", "huge", "ggmModSelect", "LoGo"))
    args$threshold <- threshold
  if (!is.null(gamma) &&
      method %in% c("EBICglasso", "huge", "ggmModSelect", "ising", "mgm"))
    args$gamma <- gamma
  do.call(fn, args)
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
