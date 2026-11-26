# Changelog

## psychnets 0.8.9

- Article titles name the package (“… with the psychnets R package”).

## psychnets 0.8.8

- The articles were revised: fuller introductions, the model before the
  data, and references to current methodological work. Their code is
  unchanged.
- The Ising and mixed graphical model articles now quote the results of
  the current version (they predated 0.8.6).
- `SRL_GPT` cites its source correctly: Vogelsmeier, Oliveira, Misiejuk,
  López-Pernas and Saqr (2025).

## psychnets 0.8.7

- [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
  documents that the closeness of an isolated node is 0 (igraph gives
  `NaN`). No number changed.

## psychnets 0.8.6

- [`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md)
  tests planarity for any graph (left-right planarity test). Before, a
  non-planar graph that was not chordal could be counted planar. The
  verdict for a TMFG is unchanged.
- [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
  and
  [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  with the native solver select each node’s penalty on glmnet’s default
  path, as
  [`IsingFit::IsingFit()`](https://rdrr.io/pkg/IsingFit/man/isingfit.html)
  and [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) do.
  `lambda_min_ratio` now defaults to glmnet’s rule. Edge weights can
  change (up to 0.0165 on six bfi items) and now agree with IsingFit and
  mgm. Pass `lambda_min_ratio = 0.01` for the old grid.

## psychnets 0.8.5

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of
  grouped and multilevel results passes `...` to each network, so
  `include_zero = TRUE` returns every node pair. It was ignored.

## psychnets 0.8.3

- A slow but finite nodewise logistic fit is no longer reported as
  separated in
  [`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md).

## psychnets 0.8.2

CRAN release: 2026-11-18

### New features

- [`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md)
  estimates the panel graphical VAR model (Epskamp, 2020) from few-wave
  panel data: temporal, contemporaneous and between-person networks.
  Edges are selected with `select = c("none", "prune", "ebic")`. New
  dataset `srl_panel` and a vignette.
- `psychnet(within = )` estimates within-person networks.
- [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  accepts two networks of one grouped
  [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  result, gains `na_method`, and reports the unadjusted p-value
  (`p_raw`).
- [`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  and
  [`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  gain `cor_method`.
- New [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods and `what =` options for stability, bootstrap, moderated MGM
  and panel GVAR results;
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  multilevel networks.

### Changes to results

- [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  stores the network-structure test as `$structure` (M) and the
  global-strength test as `$strength` (S), as in van Borkulo et al.
  (2022). `$M` and `$S` raise an error.
- [`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  computes the CS-coefficient as `bootnet::corStability()`;
  [`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  defaults to `cor_method = "pearson"`.
- [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
  refits the lambdas near the EBIC minimum before choosing one, and
  returns absent edges as exact zeros.
- [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  solves gaussian nodes to optimality (weights move by up to 1e-7).
- [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md),
  [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
  and
  [`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
  rescale `weights` to a maximum of 1, as
  [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) does.
- `net_bridge(normalize = TRUE)` divides the 2-step bridge expected
  influence by the number of nodes, as
  [`networktools::bridge()`](https://rdrr.io/pkg/networktools/man/bridge.html).
- [`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
  gives a split with an empty half a correlation of 0.
- [`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md)
  reads `value` in the moderator’s own coding.
- [`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
  warns under (quasi-)complete separation and does not certify such a
  fit.

### Other

- Errors and warnings are classed conditions (`psychnet_bad_input`, …).
- Ordered factors are kept as ordinal items;
  [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
  codes logical and two-level factor columns 0/1.
- Plots use the Okabe-Ito palette and show significance by fill and line
  type as well as colour.

## psychnets 0.5.2

CRAN release: 2026-08-23

- Faster examples: the resampling examples use `method = "pcor"`. No
  estimator changed.

## psychnets 0.5.1

### Renamed

- `casedrop_reliability()` is now
  [`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  and `network_reliability()` is now
  [`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md).
  The old names, also exported by Nestimate with a different first
  argument, are removed. Arguments, results and numbers are unchanged.

### New

- [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
  gains `lambda_path =` (a fixed penalty grid, for resampling),
  `penalize_diagonal =` (with a matching
  [`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
  argument) and `refit =` (`"unregularized"` refits the selected graph
  without a penalty).
- `mgm_fit(native = FALSE)` models categorical nodes with more than two
  levels, as [`mgm::mgm()`](https://rdrr.io/pkg/mgm/man/mgm.html) does.
  A numeric column with 10 or fewer integer values is treated as
  categorical, with a warning (`psychnet_type_promotion`).
- Any two-level coding (`{1, 2}`, `{"no", "yes"}`) gives the same fit as
  `{0, 1}`.

### Fixes

- The base nodewise solver returns exact zeros for inactive
  coefficients, so the selected edges no longer depend on the platform.
  One weak edge of the mixed-data mgm network goes.
- Resampling verbs pass each draw to the estimator as a direct call
  would (factor columns and incomplete rows kept), and keep the observed
  node set.
- [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)
  rejects pairs with fewer than two joint observations; correlation
  p-values use pair-specific sample sizes.
- Bootstrap difference tests split ties, so identical distributions give
  `p = 1`.
- Corrected node-betweenness normalisation; relative-importance
  aggregation keeps direction; rank dichotomisation no longer splits
  ties by row order.
- Unavailable certificates report `certified = NA`.
- The glmnet paths honour `nlambda` and `lambda_min_ratio`.
- [`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
  refuses an mgm network with a categorical node of more than two
  levels.
- Node order follows the input columns and is never sorted.
- Estimator numerics are stable within a minor version: a change to a
  selected edge set, or to a weight beyond 1e-10, comes with a minor
  version and is listed here.

## psychnets 0.4.3

CRAN release: 2026-07-30

- First CRAN release.
- The moderated mixed graphical model (`mgm_fit(moderators = )`) runs on
  psychnets’ own base-R solver; every estimator is base R by default,
  with `glmnet` and `glasso` as optional engines (`native = FALSE`).
- “Network reliability” and “Visualizing networks” are website articles
  at <https://pak.dynasite.org/psychnets>.
