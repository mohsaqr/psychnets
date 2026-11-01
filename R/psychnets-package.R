#' psychnets: Clean-Room Base-R Psychometric Network Estimation
#'
#' Estimates psychometric network models in pure base R, with no compiled
#' dependencies and dependency-free correctness certificates. See [psychnet()]
#' for the unified entry point to the cross-sectional estimators
#' ([cor_network()], [pcor_network()], [ebic_glasso()], [huge_network()],
#' [ggm_modselect()], [tmfg_network()], [logo_network()], [relimp_network()],
#' [ising_fit()], [ising_sampler()], [mgm_fit()]) and [panel_gvar()] for the
#' panel graphical VAR. Fitted networks are read with `as.data.frame()`,
#' [net_centralities()], [net_predict()] and [certificate()], and their accuracy
#' and stability are assessed with [net_boot()], [net_stability()] and
#' [net_compare()].
#'
#' @section Certification:
#' Regularized estimators are graded against their own objective rather than
#' an external solver. For the Gaussian graphical model, [glasso_kkt()] gives
#' the stationarity (KKT) residual: zero certifies the unique global optimum of
#' the convex objective. Every fitted network stores its residual, and
#' [certificate()] reports it as a one-row table for any estimator.
#'
#' @keywords internal
#' @importFrom stats cor cov2cor sd complete.cases qnorm pnorm
"_PACKAGE"
