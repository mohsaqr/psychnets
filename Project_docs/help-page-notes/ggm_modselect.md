# Implementation notes removed from `?ggm_modselect` (2026-11-17)

Still true of R/ggmmodselect.R; kept for maintainers.

- **Old framing**: "Equivalent in purpose to `qgraph::ggmModSelect()`, but
  pure base R and self-certified via ggm_support_kkt()".
- **Candidates**: each lambda on the path is solved by `.glasso_solve()` at
  tolerance `1e-4`; the support is `|wi| > 1e-6`, symmetrized, deduplicated
  by edge pattern. Each candidate is refit by `.ggm_fit_support()`
  (block-coordinate descent with unpenalized least squares on each node's
  neighbours, `tol_outer = 1e-8`, up to `1e4` sweeps).
- **EBIC dimension** is the number of edges in the support, not the count of
  nonzero MLE entries (a structurally present edge can estimate to zero).
- **Stepwise**: first-improvement over all pairs in upper-triangle order,
  accepting a flip when EBIC drops by more than `1e-10`; at most 200 passes.
- **Failures**: candidate/stepwise fits that error are counted
  (`.ggm_failure_log()`), reported once as `psychnet_fit_failed`; the final
  support is refit once more so `precision` always equals the fit on
  `support`. The certificate warning threshold is `1e-6`.
- **Stored fields** not named on the new page: `fit_kkt`,
  `certificate_target`.
