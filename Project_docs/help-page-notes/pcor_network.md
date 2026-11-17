# Implementation notes removed from `?pcor_network` (2026-11-17)

Still true of R/correlation.R; kept for maintainers.

- **Old framing**: "Conditional (full-order) association network ...
  Equivalent to `bootnet`'s `"pcor"` default."
- **Singular input** (`.pd_solve()`): if `solve()` fails with a "singular"
  error, or `rcond(S) < 1e-12`, a `psychnet_singular` warning is raised; on
  failure the nearest PD matrix (`.nearest_pd_cor()`) is inverted.
- **Test**: `k = p - 2` variables partialled out, so `df = n - p`; a
  pairwise fit with unequal `n_pair` refuses `alpha`.
