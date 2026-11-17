# Implementation notes removed from `?certificate` (2026-11-17)

Still true of R/certificate.R, R/result.R and the estimators; kept for
maintainers.

- **Value read**: `certificate()` reads the stored `kkt` field of the network
  (`NA` when the field is absent), except for TMFG, where it recomputes
  `tmfg_certificate(x)` from the stored `adjacency`. `kind` comes from
  `.certificate_kind()`, which `print.psychnet()` shares, so the printed label
  and `kind` cannot disagree.
- **Kinds**: `"none"` for `cor`/`pcor`; `"structural"` for `tmfg`/`relimp`;
  `"kkt"` for every other method (glasso, huge, ggm, logo, ising,
  ising_sampler, mgm).
- **certified**: `TRUE` for kind `"none"` (no optimization claim; certificate
  `NA`); `NA` when the value is `NA`; otherwise `value <= tol`.
- **Missing certificates**: the estimators store `kkt = NA` after a positive
  `threshold` (ebic_glasso, huge_network, ggm_modselect, logo_network, mgm
  LW/HW threshold that removed a coefficient), after Wald pruning in
  `ising_sampler()`, and for a non-converged `ising_sampler()` fit. The
  residual of the fit before the post-processing is kept as `fit_kkt`, and
  `certificate_target` records which (`"returned_network"`,
  `"pre_threshold_fit"`, `"pre_pruning_fit"`, `"nonconverged_fit"`).
- **External engines**: the old help said an "external-engine fit without an
  available residual reports NA". `ebic_glasso(native = FALSE)` stores the
  residual of the glasso package's solution (its looser tolerance, ~1e-4); the
  moderated mgm on glmnet stores `NA`, but that object is not a `psychnet` and
  `certificate()` refuses it. No reachable `psychnet` currently has an
  engine-missing certificate, so the new page does not mention it.
- **Old wording**: "Regularized or constrained estimators in psychnets
  self-certify ... a KKT / stationarity residual ... near machine zero for a
  correctly solved problem." Typical values on the default engines are of the order
  1e-10 (e.g. 1.1e-10 for the `?certificate` example glasso fit).
- **Contradiction fixed**: the old description said cor/pcor "report NA with
  kind none" and the old @return said they are "reported as TRUE". The code
  gives `certificate = NA`, `certified = TRUE`; the new page states both.
