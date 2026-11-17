# Implementation notes removed from `?glasso_kkt` (2026-11-17)

Still true of R/glasso.R; kept for maintainers.

- **Old framing**: "A dependency-free correctness certificate ... certifies
  correctness independently of any reference solver." The new page says what
  the conditions are and what a near-zero value shows.
- **Feasibility checks**: `Inf` when `max|theta - t(theta)| > 1e-8`, when the
  smallest eigenvalue is `<= 0`, or when `solve(theta)` fails.
- **Value**: `max(diag_v, v_a, v_i)` over the upper triangle, with
  `v_i = max(pmax(|W_ij - S_ij| - rho, 0))` for inactive entries.
