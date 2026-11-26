# psychnets 0.8.9

* Article titles name the package ("... with the psychnets R package").

# psychnets 0.8.8

* The articles were revised: fuller introductions, the model before the data,
  and references to current methodological work. Their code is unchanged.
* The Ising and mixed graphical model articles now quote the results of the
  current version (they predated 0.8.6).
* `SRL_GPT` cites its source correctly: Vogelsmeier, Oliveira, Misiejuk,
  López-Pernas and Saqr (2025).

# psychnets 0.8.7

* `net_centralities()` documents that the closeness of an isolated node is 0
  (igraph gives `NaN`). No number changed.

# psychnets 0.8.6

* `tmfg_certificate()` tests planarity for any graph (left-right planarity
  test). Before, a non-planar graph that was not chordal could be counted
  planar. The verdict for a TMFG is unchanged.
* `ising_fit()` and `mgm_fit()` with the native solver select each node's
  penalty on glmnet's default path, as `IsingFit::IsingFit()` and `mgm::mgm()`
  do. `lambda_min_ratio` now defaults to glmnet's rule. Edge weights can change
  (up to 0.0165 on six bfi items) and now agree with IsingFit and mgm. Pass
  `lambda_min_ratio = 0.01` for the old grid.

# psychnets 0.8.5

* `as.data.frame()` of grouped and multilevel results passes `...` to each
  network, so `include_zero = TRUE` returns every node pair. It was ignored.

# psychnets 0.8.3

* A slow but finite nodewise logistic fit is no longer reported as separated
  in `ising_sampler()`.

# psychnets 0.8.2

## New features

* `panel_gvar()` estimates the panel graphical VAR model (Epskamp, 2020) from
  few-wave panel data: temporal, contemporaneous and between-person networks.
  Edges are selected with `select = c("none", "prune", "ebic")`. New dataset
  `srl_panel` and a vignette.
* `psychnet(within = )` estimates within-person networks.
* `net_compare()` accepts two networks of one grouped `psychnet()` result,
  gains `na_method`, and reports the unadjusted p-value (`p_raw`).
* `net_stability()` and `net_casedrop_reliability()` gain `cor_method`.
* New `as.data.frame()` methods and `what =` options for stability, bootstrap,
  moderated MGM and panel GVAR results; `summary()` and `plot()` for
  multilevel networks.

## Changes to results

* `net_compare()` stores the network-structure test as `$structure` (M) and
  the global-strength test as `$strength` (S), as in van Borkulo et al.
  (2022). `$M` and `$S` raise an error.
* `net_stability()` computes the CS-coefficient as `bootnet::corStability()`;
  `net_casedrop_reliability()` defaults to `cor_method = "pearson"`.
* `ebic_glasso()` refits the lambdas near the EBIC minimum before choosing one,
  and returns absent edges as exact zeros.
* `mgm_fit()` solves gaussian nodes to optimality (weights move by up to 1e-7).
* `mgm_fit()`, `ising_fit()` and `ising_sampler()` rescale `weights` to a
  maximum of 1, as `mgm::mgm()` does.
* `net_bridge(normalize = TRUE)` divides the 2-step bridge expected influence
  by the number of nodes, as `networktools::bridge()`.
* `net_split_reliability()` gives a split with an empty half a correlation
  of 0.
* `condition()` reads `value` in the moderator's own coding.
* `ising_sampler()` warns under (quasi-)complete separation and does not
  certify such a fit.

## Other

* Errors and warnings are classed conditions (`psychnet_bad_input`, ...).
* Ordered factors are kept as ordinal items; `ising_fit()` codes logical and
  two-level factor columns 0/1.
* Plots use the Okabe-Ito palette and show significance by fill and line type
  as well as colour.

# psychnets 0.5.2

* Faster examples: the resampling examples use `method = "pcor"`. No
  estimator changed.

# psychnets 0.5.1

## Renamed

* `casedrop_reliability()` is now `net_casedrop_reliability()` and
  `network_reliability()` is now `net_split_reliability()`. The old names,
  also exported by Nestimate with a different first argument, are removed.
  Arguments, results and numbers are unchanged.

## New

* `ebic_glasso()` gains `lambda_path =` (a fixed penalty grid, for resampling),
  `penalize_diagonal =` (with a matching `glasso_kkt()` argument) and
  `refit =` (`"unregularized"` refits the selected graph without a penalty).
* `mgm_fit(native = FALSE)` models categorical nodes with more than two
  levels, as `mgm::mgm()` does. A numeric column with 10 or fewer integer
  values is treated as categorical, with a warning
  (`psychnet_type_promotion`).
* Any two-level coding (`{1, 2}`, `{"no", "yes"}`) gives the same fit as
  `{0, 1}`.

## Fixes

* The base nodewise solver returns exact zeros for inactive coefficients, so
  the selected edges no longer depend on the platform. One weak edge of the
  mixed-data mgm network goes.
* Resampling verbs pass each draw to the estimator as a direct call would
  (factor columns and incomplete rows kept), and keep the observed node set.
* `cor_auto()` rejects pairs with fewer than two joint observations;
  correlation p-values use pair-specific sample sizes.
* Bootstrap difference tests split ties, so identical distributions give
  `p = 1`.
* Corrected node-betweenness normalisation; relative-importance aggregation
  keeps direction; rank dichotomisation no longer splits ties by row order.
* Unavailable certificates report `certified = NA`.
* The glmnet paths honour `nlambda` and `lambda_min_ratio`.
* `net_predict()` refuses an mgm network with a categorical node of more than
  two levels.
* Node order follows the input columns and is never sorted.
* Estimator numerics are stable within a minor version: a change to a selected
  edge set, or to a weight beyond 1e-10, comes with a minor version and is
  listed here.

# psychnets 0.4.3

* First CRAN release.
* The moderated mixed graphical model (`mgm_fit(moderators = )`) runs on
  psychnets' own base-R solver; every estimator is base R by default, with
  `glmnet` and `glasso` as optional engines (`native = FALSE`).
* "Network reliability" and "Visualizing networks" are website articles at
  <https://pak.dynasite.org/psychnets>.
