# Implementation notes removed from `?glm_lasso_kkt` (2026-11-17)

Still true of R/lasso_glm.R; kept for maintainers.

- **Old framing**: "Dependency-free correctness certificate for a nodewise
  lasso, analogous to glasso_kkt()".
- **Contradiction fixed**: the old help gave only the coefficient
  conditions and said `X` must be standardized. The code also checks the
  intercept condition (`|sum(w (y - mu))| / n`) and accepts any `X`; the
  conditions refer to the scale of `X` (the fits penalize standardized
  predictors). Gradients are normalized by the row count `n`, matching
  `.wls_lasso()`, also when weights are given.
- The fits rescale observation weights to a maximum of 1
  (`.check_weights()`); `glm_lasso_kkt()` uses the weights as given.
- `ising_sampler()` calls it with `lambda = 0` and family binomial.
