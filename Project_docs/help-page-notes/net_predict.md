# Implementation notes removed from `?net_predict` (2026-11-17)

Still true of R/predictability.R; kept for maintainers.

- Same quantities as the predictability reported by the `mgm` package (R2,
  CC, nCC); implemented in base R without it. (The source comment names
  `mgm::net_predict()`, which is not an mgm function.)
- GGM closed form uses `S_jj = diag(solve(theta))` when the diagonal was
  penalized (`penalize_diagonal`), else the diagonal of the stored
  `cor_matrix`, else 1. The table is not clamped; `node_predictability()` and
  the predictability stored on the node table clamp to \[0, 1\].
- Nodewise models: coefficients `beta_std` on the standardized predictors
  (`center`, `scale`); binary prediction thresholds the logistic probability
  at 0.5; `nCC = (CC - m) / (1 - m)` with `m = max(p, 1 - p)`, 0 when the node
  is constant. Gaussian R2 is computed on the response scale the estimator
  used (`resp_center`, `resp_scale` for mgm).
- Category coding: logical columns are coded by their training levels; numeric
  0/1 for a node trained on labels is accepted as the estimator's codes.
- A non-`psychnet` `x` fails an unclassed `stopifnot()` (not
  `psychnet_bad_input`).
- An mgm network with a multi-level categorical node is refused (no
  one-coefficient-per-variable design).
