# Implementation notes removed from `?net_split_reliability` (2026-11-17)

Still true of R/reliability.R; kept for maintainers.

- **Provenance and name**: adapted from Nestimate. Renamed from
  `network_reliability()` (2026-08-07) because Nestimate exports a function of
  that name with a different signature; the old name was removed, not aliased.
- **Halves**: `n_half = max(2, round(n * split))`; the second half is the
  complement. Halves are drawn in the parent, one per iteration.
- **Full-sample fit**: run only when `labels` is `NULL`, to read the node set
  both halves are re-expanded to (a factor column is a node for mgm, dropped by
  glasso). With `labels` supplied there is no full-sample fit.
- **Correlation rule**: `.psn_cor0(constant = "na")`: `NA` when both halves are
  constant, 0 when exactly one is.
- **Percentiles**: `lower`/`upper` are the 2.5% and 97.5% `quantile()` values
  over the defined splits. At least two splits must succeed.
- **Stored attributes**: `iterations` (data frame of per-split metrics),
  `iter`, `split`, `cor_method`, `method`, `n`, `n_success`, `n_failed`.
