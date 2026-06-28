# Unified front door, mirroring bootnet::estimateNetwork(data, default = ...).

#' Estimate a psychometric network
#'
#' The package's main entry point: routes to the requested estimator and returns
#' a common `psychnet` object, so callers can swap estimators without rewiring
#' downstream code.
#'
#' `method` speaks the package's own short vocabulary -- `"glasso"`, `"ggm"`,
#' `"tmfg"`, `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`, `"huge"`,
#' `"mgm"`, `"cor"`, `"pcor"`. For interoperability it **also** accepts the
#' `qgraph`/`bootnet` spellings (`"EBICglasso"`, `"ggmModSelect"`, `"TMFG"`,
#' `"LoGo"`, `"IsingFit"`, `"IsingSampler"`), which resolve to the same
#' estimators. Whichever you pass in, the stored `$method` is the short name.
#'
#' @param data Numeric data frame or matrix (rows = observations).
#' @param method Estimator. One of `"glasso"` (default), `"ggm"`, `"tmfg"`,
#'   `"logo"`, `"relimp"`, `"ising"`, `"ising_sampler"`, `"huge"`, `"mgm"`,
#'   `"cor"`, `"pcor"`. The `qgraph`/`bootnet` names are accepted as aliases
#'   (e.g. `"EBICglasso"` -> `"glasso"`, `"ggmModSelect"` -> `"ggm"`,
#'   `"LoGo"` -> `"logo"`).
#' @param threshold Absolute-weight threshold below which edges are zeroed
#'   (forwarded only to the methods that take it: `cor`, `pcor`, `glasso`,
#'   `huge`, `ggm`, `logo`).
#' @param gamma EBIC hyperparameter. `NULL` (default) keeps each method's own
#'   default (0.5 for the regularized Gaussian graphical models, 0 for `ggm`,
#'   0.25 for `ising`/`mgm`); set it to override. Forwarded only to the
#'   regularized methods.
#' @param labels Optional node labels.
#' @param ... Passed to the underlying estimator (e.g. `cor_method=` for the
#'   correlation-based methods, `npn=` for `"huge"`, `rule=` for the Ising
#'   methods, `alpha=` for the correlation / `"ising_sampler"` methods).
#' @return A `psychnet` object.
#' @examples
#' x <- matrix(stats::rnorm(200 * 5), 200, 5)
#' psychnet(x, method = "glasso")
#' psychnet(x, method = "pcor", cor_method = "spearman")
#' psychnet(x, method = "EBICglasso")  # qgraph alias, same result
#' @export
psychnet <- function(data,
                     method = c("glasso", "cor", "pcor", "ising", "mgm",
                                "huge", "ggm", "tmfg", "logo", "relimp",
                                "ising_sampler"),
                     threshold = 0, gamma = NULL, labels = NULL, ...) {
  method <- .resolve_method(method)
  fn <- switch(
    method,
    cor = cor_network, pcor = pcor_network, glasso = ebic_glasso,
    ising = ising_fit, mgm = mgm_fit, huge = huge_network,
    ggm = ggm_modselect, tmfg = tmfg_network, logo = logo_network,
    relimp = relimp_network, ising_sampler = ising_sampler,
    stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  )
  # Forward each shared argument only to the methods that accept it, so the
  # callee's own defaults stand otherwise (notably gamma: ising/mgm = 0.25).
  args <- c(list(data = data, labels = labels), list(...))
  if (method %in% c("cor", "pcor", "glasso", "huge", "ggm", "logo"))
    args$threshold <- threshold
  if (!is.null(gamma) &&
      method %in% c("glasso", "huge", "ggm", "ising", "mgm"))
    args$gamma <- gamma
  do.call(fn, args)
}

# Resolve a method name to its canonical short form. Accepts the package's own
# short names and the qgraph/bootnet spellings (dual vocabulary); matching is
# case- and separator-insensitive, so "EBICglasso", "ebic_glasso" and "glasso"
# all land on "glasso".
#' @noRd
.resolve_method <- function(method) {
  if (length(method) > 1L) method <- method[1L]
  key <- gsub("[^a-z]", "", tolower(method))         # strip case/underscores
  canon <- c(
    cor = "cor", correlation = "cor",
    pcor = "pcor", partial = "pcor", partialcor = "pcor",
    glasso = "glasso", ebicglasso = "glasso", ebic = "glasso",
    huge = "huge", npn = "huge", nonparanormal = "huge",
    ggm = "ggm", ggms = "ggm", ggmmodselect = "ggm", modselect = "ggm",
    stepwise = "ggm",
    tmfg = "tmfg",
    logo = "logo",
    relimp = "relimp", relativeimportance = "relimp", lmg = "relimp",
    ising = "ising", isingfit = "ising",
    isingsampler = "ising_sampler",
    mgm = "mgm"
  )
  out <- canon[key]
  if (is.na(out)) stop(sprintf("Unknown method '%s'.", method), call. = FALSE)
  unname(out)
}
