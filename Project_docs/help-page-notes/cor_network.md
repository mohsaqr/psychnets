# Implementation notes removed from `?cor_network` (2026-11-17)

Still true of R/correlation.R; kept for maintainers.

- **Old framing**: "Marginal (zero-order) association network ... Equivalent
  to `bootnet`'s `"cor"` default."
- **Tests**: `.cor_pvalues()` uses `t = r sqrt(df / max(1 - r^2, 1e-12))`;
  with `data` the per-pair co-observation matrix `n_pair` is used as `n`.
  Kendall p-values come from `stats::cor.test(method = "kendall",
  exact = FALSE)` per pair. Adjustment is `p.adjust()` over the upper
  triangle, mirrored. Stored `p_values` are unadjusted.
- **n**: with `cor_matrix` and `alpha`, `n` must be a whole number >= 3
  (`.psn_check_n(n, min = 3)`).
