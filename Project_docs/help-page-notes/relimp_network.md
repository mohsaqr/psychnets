# Implementation notes removed from `?relimp_network` (2026-11-17)

Still true of R/relimp.R; kept for maintainers.

- **Old framing**: "Equivalent in purpose to
  `relaimpo::calc.relimp(type = "lmg")` applied nodewise / `bootnet`'s
  `"relimp"` default, pure base R and self-certified via lmg_certificate()".
- **Computation**: subset R^2 = `S[j,A] solve(S[A,A]) S[A,j]` over all
  `2^(p-1) - 1` nonempty predictor subsets (a pseudo-inverse `MASS_ginv()`
  when `solve()` fails); Shapley weights via `lgamma` (finite up to m = 20);
  literal loops, the documented exception to the no-loop rule.
- **normalized = TRUE**: columns whose share sum is below
  `sqrt(.Machine$double.eps) * max(1, max(colSums))` are left at zero rather
  than divided.
- **cor_matrix errors** from `.check_cor_matrix()` are re-raised with the
  prefix "relimp_network() requires a symmetric positive-semidefinite
  matrix:".
- The old help wrote the identity as `colSums($weights) == $r2`.
