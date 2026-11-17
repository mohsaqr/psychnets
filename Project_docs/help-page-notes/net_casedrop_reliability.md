# Implementation notes removed from `?net_casedrop_reliability` (2026-11-17)

Still true of R/reliability.R; kept for maintainers.

- **Provenance and name**: adapted from Nestimate. Renamed from
  `casedrop_reliability()` (2026-08-07) because Nestimate exports a function of
  that name with a different first argument; the old name was removed, not
  aliased.
- **Metrics** (`.psn_edge_metrics()`): `mean_abs_dev`, `median_abs_dev`,
  `correlation` (`.psn_cor0()`, `constant = "zero"`), `max_abs_dev`, on the
  upper-triangle edge vector (every off-diagonal cell when directed,
  `.psn_edge_vector()`).
- **Summaries** (`.psn_summarise_draws()`): over finite values only; `sd` needs
  two values; `NA`, never `NaN`, when too few are defined.
- **Warnings**: a constant full-sample edge vector reuses the class
  `psychnet_constant_centrality` (the same class as `net_stability()`).
  `psychnet_undefined_statistic` counts estimated subsets with a non-finite
  correlation.
- **Stored attributes**: `cs`, `correlations` (iter x proportions),
  `prop_above`, `threshold`, `certainty`, `iter`, `cor_method`, `method`,
  `n_cases`, `n_edges`, `n_success`, `n_failed`. There is no accessor for the
  CS-coefficient other than `print()`.
- Subsets and CS rule as in `net_stability()` (see `net_stability.md`); at
  least one subset must succeed.
