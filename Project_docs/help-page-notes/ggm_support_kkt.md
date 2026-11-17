# Implementation notes removed from `?ggm_support_kkt` (2026-11-17)

Still true of R/glasso.R; kept for maintainers.

- **Old framing**: "Constrained Gaussian-MRF (graph-restricted MLE)
  stationarity residual"; "maximum-likelihood / maximum-entropy conditions";
  "certifies the constrained optimum with no reference solver".
- **active_tol is not used**: it is validated (`.psn_check_number`, >= 0) but
  the off-graph term is `max(abs(theta[!support]))` with no tolerance. The old
  help said it is the "magnitude above which an off-support entry counts as a
  nonzero violation"; the new page documents the code (accepted, no effect).
- **Feasibility checks**: as `glasso_kkt()` (asymmetry > 1e-8, eigenvalue
  <= 0, failed `solve()` all return `Inf`). `support` is symmetrized with
  `support | t(support)`.
