# Implementation notes removed from `?net_compare` (2026-11-17)

Still true of R/nct.R; kept for maintainers.

- **Provenance**: ported from Nestimate. Each network is `ebic_glasso()` on the
  package's own solver; the nearest-correlation projection is base R
  (`.nearest_pd_cor()`: one eigenvalue clamp at `1e-8` and renormalization to a
  unit diagonal, Higham 2002 single step), not `Matrix::nearPD`.
- **Parity with NetworkComparisonTest**: under the same seed both permute the
  same rows (independent: `sample(n1 + n2, n1)`; paired:
  `sample(c(1, 2), n1, replace = TRUE)` per pair), so p-values agree and
  estimates differ only by the glasso solver's convergence tolerance. Fixtures:
  `nct-reference.rds` (parity), `nct-baseline.rds` (frozen complete-data
  output). Calibration and tie audits live in
  `local_testing_and_equivalence/nct/`.
- **p-values**: `(sum(perm >= observed) + 1) / (iter + 1)` for M, S and each
  edge; only the per-edge p-values are adjusted.
- **Letters (version history)**: before 2026-10-12 (psychnets 0.6.0) the
  fields `$M` and `$S` held the two global tests the other way round from van
  Borkulo et al. (2022). The fields are now `structure` (M) and `strength` (S,
  with `by_group`); `$.psychnet_nct` raises `psychnet_renamed` for `$M`/`$S`
  with a message naming the replacement. (The old help described this swap in
  the return value.)
- **Stored fields**: `nw1`, `nw2`, `structure` and `strength` (`observed`,
  `perm`, `p_value`), `strength$by_group`, `E` (`observed`, `perm` iter x
  edges, `p_value`, `p_raw`, `edge_names`), `groups` (argument names, or the
  levels for a group comparison), `n_iter`, `n`, `paired`, `params`.
- **Missing data**: `.nct_cor()` uses `pairwise.complete.obs` when any value is
  missing; returns `NULL` (-> `psychnet_degenerate_sample`) when a pair shares
  fewer than two rows. The sample size is the median pairwise co-observation
  count, at least 2. A missing value is never zero-filled (it once was, and one
  `NA` erased a variable).
- **Constant variables**: an undefined correlation caused by a variable without
  variation over the shared rows is set to 0 (isolated node); the
  "standard deviation is zero" warning of `stats::cor()` is muffled by its
  translated message. Any other non-finite correlation is a
  `psychnet_bad_input` error.
- **Group mode** (`.nct_group()`, `.nct_members()`): uses the group's `gamma`
  (default 0.5) and `na_method` unless supplied; warns
  `psychnet_estimator_mismatch` for non-glasso groups, other dots, a nonzero
  threshold, or overridden settings. Paired-by-person comparisons sort the
  common persons so the pairs do not depend on row order.
- **abs**: the old help said "compare absolute edge weights"; in the code
  `abs` changes only S (sum of absolute vs signed weights); M and E always use
  `abs(nw1 - nw2)`. The new page says so.
- **CLAUDE.md drift**: CLAUDE.md mentions `as.data.frame(x, what = "global")`;
  `as.data.frame.psychnet_nct()` has no `what` argument. The global tests are
  returned by `summary()`, and the page documents that.
