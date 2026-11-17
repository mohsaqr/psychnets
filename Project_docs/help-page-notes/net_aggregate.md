# Implementation notes removed from `?net_aggregate` (2026-11-17)

Still true of R/aggregate.R; kept for maintainers.

- Inspired by Nestimate's mgm / `mcml_pc` aggregation. The PCA composite
  matches `networktools::net_reduce()`; the canonical method matches
  `stats::cancor()`.
- Columns whose community is `NA` are removed before any preprocessing
  (`.psn_left_out_columns()`), so an identifier neither triggers listwise
  deletion nor needs to be numeric.
- `scale` affects only the score composites. `"rv"` and `"canonical"` use the
  unscaled data (`stats::cov(mat)` for RV).
- `.psn_align()` flips a PCA/factor score whose mean correlation with the
  community's items is negative; a constant composite is left as is.
- `"loadings"` fits the item-level network once with
  `psychnet(mat, method = estimator, ...)`; weights below
  `sqrt(.Machine$double.eps)` in total fall back to equal weights.
- `"average"` averages all cross-community entries of the item-level weight
  matrix, zeros included; it is directed when that network is.
- Association networks are built with `.new_psychnet()`, method
  `"aggregate_<method>"`, `n_obs` = retained rows. `cancor()` errors on
  rank-deficient or wider-than-tall blocks; that pair gets weight 0.
