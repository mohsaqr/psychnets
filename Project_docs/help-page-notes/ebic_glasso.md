# Implementation notes removed from `?ebic_glasso` (2026-11-17)

Still true of R/glasso.R; kept for maintainers.

- **Old framing**: "pure base R and self-certified", "refits the chosen
  penalty to machine precision so the fitted precision matrix is the
  certified global optimum", "Equivalent in purpose to
  `qgraph::EBICglasso()` / `bootnet`'s `"EBICglasso"` default". The new page
  names `qgraph::EBICglasso()` once, at the end of the description.
- **Selection** (`.select_ebic()`): scan the lambda path at `scan_tol = 1e-4`
  with warm starts; then, unless `refit = FALSE` or `native = FALSE`, refit
  every lambda whose scan EBIC is within two edge penalties
  (`2 * (log n + 4 gamma log p)`) of the scan minimum at `refit_tol = 1e-8`
  (inner tolerance `1e-10`), rescore, and keep the lowest exact EBIC; repeat
  if the minimum moves to an unconfirmed row. Edges are counted at
  `|theta_ij| > 1e-10`; returned entries at or below `1e-10` are set to exact
  zeros. A lambda whose fit errors or is not positive definite scores
  `EBIC = Inf` (visible in `ebic_path`).
- **refit = FALSE**: the old help called it "a bit-compatibility switch, not
  a speed one ... a consumer holding frozen baselines against an unrefitted
  path scan needs the scanned matrix to reproduce them", and said its
  certificate reflects "the loose scan tolerance (~1e-4) rather than the
  refit's ~1e-11". On `0.4^|i-j|` (p = 6, n = 250) the stored certificate
  was 2.2e-6 (refit = FALSE) vs 1.1e-10 (default), measured 2026-11-17.
- **native = FALSE**: delegates each fixed-penalty solve to
  `glasso::glasso()` (`thr = tol`), symmetrizes `wi`; the confirm step is
  skipped and the final refit stays at `scan_tol`, so the result is what the
  glasso package returns. The old help claimed "byte-identical
  `glasso`/`qgraph` output". On the same example the certificate was 1.5e-6,
  i.e. above `certificate()`'s default `tol = 1e-6`.
- **penalize_diagonal = TRUE** holds the working covariance diagonal at
  `diag(S) + lambda` (instead of `diag(S)`), and the certificate uses
  `W_ii = S_ii + rho`.
- **lambda_path** is echoed back unchanged (stored as `lambda_path`); the old
  help explained that recomputing the path per resample lets the grid drift
  with each sample's largest absolute correlation.
- **Stored fields** not named on the new page: `na_method`, `native`,
  `penalize_diagonal`, `refit`, `fit_kkt` (the certificate before
  thresholding), `certificate_target` (`"returned_network"` or
  `"pre_threshold_fit"`).
