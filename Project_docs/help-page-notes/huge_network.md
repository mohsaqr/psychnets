# Implementation notes removed from `?huge_network` (2026-11-17)

Still true of R/huge.R and R/glasso.R; kept for maintainers.

- **Old framing**: "refits to the certified optimum. Equivalent in purpose to
  `huge::huge()` (nonparanormal) / `bootnet`'s `"huge"` default, but pure base
  R and self-certified via glasso_kkt() on the transformed correlation."
- **Pipeline**: `.npn_scores()` ranks each column over its observed values
  (with no missing data this is exactly `qnorm(rank / (n + 1))`); the skeptic
  is `2 sin(pi/6 * spearman)` with unit diagonal. The matrix is projected
  with `.nearest_pd_cor()` when its smallest eigenvalue is below `1e-8`, and
  always for pairwise data with missing values. `n` is the median pairwise
  co-observation count with missing data.
- **Selection** reuses `.select_ebic()` with its defaults (refit = TRUE,
  penalize_diagonal = FALSE); there is no `lambda_path`, `refit` or
  `penalize_diagonal` argument. A matrix with no association (max
  |S_ij| <= 1e-12) returns `.empty_glasso()`.
- **Stored fields** not named on the new page: `fit_kkt`,
  `certificate_target`.
