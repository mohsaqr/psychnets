# Implementation notes removed from `?ising_sampler` (2026-11-17)

Still true of R/ising_unreg.R and R/lasso_glm.R; kept for maintainers.

- **Solver**: each node is `.glm_lasso_fit(..., lambda = 0)` on standardized
  predictors (the lambda = 0 limit of the penalized IRLS kernel); the
  certificate is `glm_lasso_kkt(..., lambda = 0)`, i.e. the maximum-likelihood
  score residual. The old help called this "self-certified".
- **Wald p-values** (`.logit_wald_p()`): from the Fisher information of the
  intercept-plus-slopes design (weighted when `weights` are given), with the
  working weights floored at 1e-10. A singular information falls back to the
  base-R pseudo-inverse `MASS_ginv()` with a `psychnet_singular_information`
  warning. The p-values are scale-invariant, so they are computed on the
  standardized fit.
- **Pruning**: `p.adjust()` is applied over all off-diagonal entries of the
  asymmetric p-value matrix (p(p-1) tests), then entries with adjusted p >
  alpha are zeroed in both the raw and standardized coefficient matrices
  before the AND/OR rule. No refit. Thresholds are recomputed from the
  pruned raw coefficients: `b0_std - B %*% center`. The stored `p_values`
  are the unadjusted ones.
- **Certificate bookkeeping**: `kkt = NA` when pruning zeroed a nonzero
  coefficient or any node did not converge; `fit_kkt` keeps the unpruned
  worst residual; `certificate_target` is `"nonconverged_fit"`,
  `"pre_pruning_fit"` or `"returned_network"` (in that order of precedence).
- **Separation**: detected by a drift test on a continuation of the IRLS
  fit (`separation` flag in `.glm_lasso_fit()`); it also sets
  `converged = FALSE`. Under separation the score tends to zero while the
  coefficients diverge, so a small residual would certify nothing; hence the
  NA certificate. When separation occurs, only `psychnet_separation` is
  raised (not also `psychnet_no_converge`).
