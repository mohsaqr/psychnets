# Implementation notes removed from `?net_stability` (2026-11-17)

Still true of R/stability.R and R/reliability.R; kept for maintainers.

- **Provenance**: adapted from Nestimate; the CS rule follows
  `bootnet::corStability()` as of bootnet 1.7.1 (`.psn_cs_coefficient()`).
- **Correlation rule** (`.psn_cor0()`, bootnet's `cor0()`): `NA` when either
  side has fewer than two observed values or a non-finite value; 0 when either
  side is constant; otherwise `stats::cor()`. Constancy is judged by
  `.psn_is_constant()`: range <= `100 * .Machine$double.eps * max(abs(v))`,
  i.e. relative to the vector's own magnitude.
- **Subsets**: `.psn_keep_sizes()` keeps `max(2, round(n * (1 - drop_prop)))`
  cases; subsets are drawn in the parent, proportion outer and iteration
  inner (`.psn_casedrop_indices()`), so a seeded call reproduces the former
  nested loop. Each draw runs under `.psn_try_draw()`; at least one subset
  must succeed (`min_success = 1`).
- **Stored fields**: `cs`, `table`, `drop_prop`, `threshold`, `certainty`,
  `cor_method`, `correlations` (iter x proportions per measure), `iter`,
  `method`, `n_success`, `n_failed`. `as.data.frame(what = "correlations")`
  errors for objects made by older versions without `correlations`; `print()`
  falls back to `"spearman"` when an old object lacks `cor_method`.
- **CLAUDE.md drift**: the framework table in CLAUDE.md describes the
  correlation as Spearman; the code's default is `cor_method = "pearson"`
  (as bootnet), and the page documents the code.
