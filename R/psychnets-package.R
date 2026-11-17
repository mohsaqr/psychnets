#' psychnets: Psychometric Network Estimation
#'
#' @description
#' `psychnets` estimates psychometric networks and assesses their accuracy,
#' stability and differences. A psychometric network represents
#' psychological variables, such as items, symptoms or behaviours, as nodes
#' and their associations as edges; it shows which variables are directly
#' related, which are central, and how the structure differs between groups.
#'
#' The cross-sectional estimators are reached through [psychnet()] or called
#' directly: Gaussian graphical models ([ebic_glasso()], [ggm_modselect()],
#' [huge_network()]), information-filtering networks ([tmfg_network()],
#' [logo_network()]), relative-importance networks ([relimp_network()]),
#' correlation and partial-correlation networks ([cor_network()],
#' [pcor_network()]), Ising models for binary data ([ising_fit()],
#' [ising_sampler()]) and mixed graphical models ([mgm_fit()]). Panel data from
#' a few waves are modelled with the panel graphical VAR ([panel_gvar()]).
#' Fitted networks are read with `as.data.frame()`, [net_centralities()],
#' [net_predict()] and [certificate()]; their accuracy, stability and group
#' differences are assessed with [net_boot()], [net_stability()] and
#' [net_compare()].
#'
#' @section Certificates:
#' Each regularized or constrained estimator stores how far its network is
#' from satisfying the optimality conditions of its objective, so that the
#' fitted network can be checked without a second solver. [certificate()]
#' reports this value for any fitted network.
#'
#' @keywords internal
#' @importFrom stats cor cov2cor sd complete.cases qnorm pnorm
"_PACKAGE"
