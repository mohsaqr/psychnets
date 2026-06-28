#' psychnets: Clean-Room Base-R Psychometric Network Estimation
#'
#' Estimates cross-sectional psychometric network models in pure base R, with
#' no compiled dependencies and dependency-free correctness certificates. See
#' [psychnet()] for the unified entry point and [ebic_glasso()],
#' [cor_network()], [pcor_network()], [ising_fit()], and [mgm_fit()] for the
#' individual estimators.
#'
#' @section Certification:
#' Regularized estimators are graded against their own convex objective rather
#' than an external solver. For the Gaussian graphical model, [glasso_kkt()]
#' returns the stationarity (KKT) residual: zero certifies the unique global
#' optimum. Every [ebic_glasso()] result carries this value in `$kkt`.
#'
#' @keywords internal
#' @importFrom stats cor cov2cor sd complete.cases qnorm pnorm
"_PACKAGE"
