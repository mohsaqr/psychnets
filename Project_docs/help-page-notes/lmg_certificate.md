# Implementation notes removed from `?lmg_certificate` (2026-11-17)

Still true of R/relimp.R; kept for maintainers.

- Reads `raw_importance` when present (normalized fits), else `weights`, and
  calls `.lmg_residual(weights, r2, cor_matrix)`, the same function that
  fills the stored `kkt` at fit time.
