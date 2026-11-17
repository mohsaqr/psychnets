# Implementation notes removed from `?panel_gvar` (2026-11-17)

Still true of R/panelgvar.R and R/panelgvar_ebic.R; kept for maintainers.

- **Model** (psychonetrics `panelvar` with GGM within/between): stacked
  (wave x variable) covariance is block Toeplitz, `Sigma[t, s] = B^(t-s)
  Sigma0 + Sigma_between` (t >= s), `vec Sigma0 = (I - B %x% B)^-1 vec
  Sigma_zeta`, with `Sigma_zeta` and `Sigma_between` each parameterized as
  `D (I - Omega)^-1 D`; one mean vector shared by all waves. Fit function
  `tr(SK) + r'Kr - log det K` on the ML covariance (denominator n), the
  same as psychonetrics' objective. Parameter order is psychonetrics':
  `mu`, `omega_zeta_within`, `delta_zeta_within`, `beta` (column-major,
  `beta[i, j]` = effect of j at t-1 on i at t), `omega_zeta_between`,
  `delta_zeta_between`.
- **Optimizer**: `nlminb` with the reverse-mode analytic gradient
  (`.pgv_gradient_full`: backpropagation through the lag recursion plus an
  adjoint Lyapunov solve), `rel.tol = 1e-14`, delta lower bound 1e-8, then up
  to 25 Fisher-scoring polish steps with step halving and a pseudo-inverse
  step (so a boundary fit cannot crash it). Domain guard in `.pgv_state()`.
- **Starts** (`.pgv_best_fit`): the saturated fit uses three starts
  (`moments`, `simple`, fitted `baseline`); a pruning refit uses `previous`
  and `baseline`; the EBIC final solve uses the path solution plus the
  saturated fit and the moments start; an EBIC refit uses `selected` and
  `baseline`. Agreement: objective within 1e-6 relative of the best, judged
  on the standardized scale. "no start reached a feasible optimum" raises
  `psychnet_no_converge`.
- **Scaling**: optimization, inference, the certificate, step-up MIs,
  boundary checks, fit indices and the variance split all run on variables
  standardized by one pooled SD per variable (`sctx`); estimates and SEs map
  back with `.pgv_scale_factor()` (mu, delta x s; `beta[i, j]` x s_i/s_j;
  omegas invariant). The data-scale objective differs by the constant
  `2 * waves * sum(log(s))`.
- **Certificate**: ML fits store max |score| over the free parameters by
  forward mode (`.pgv_jacobian` + `.pgv_score`), independent of the
  optimizer's reverse-mode gradient (the two agree to ~1e-15). A delta on its
  lower bound counts as an active constraint (score >= 0 there). EBIC fits
  without refit store the lasso KKT residual `.pgv_pen_kkt`. Warning
  threshold for an interior solution: 1e-6.
- **Inference**: SEs = `sqrt(diag(solve(n * I)))` with I the per-observation
  expected information (psychonetrics' `getVCOV()` for non-Wishart ML).
  Unidentified coordinates (`.pgv_identified_se`, eigenspace on the
  correlation scale of the information) get `NA` SEs in the final fit and
  infinite SEs while pruning (so they are pruned), not a pseudo-inverse.
- **EBIC path** (`.pgv_ebic_path`, psychonetrics PML): penalty on the data
  scale = weighted lasso on the standardized scale with weights `to_data`;
  lambda_max = max |gradient| / w at the fitted baseline; FISTA with
  restart and domain-safe backtracking, then an active-set Newton polish,
  repeated if the zero set moves (default 60 FISTA iterations per round),
  stopping on a boundary ridge. Path fits use warm starts; the selected
  lambda is re-solved at `tol = 1e-10` (3000 FISTA iterations from
  the path start, 400 from the others) and reselected until the EBIC minimum
  is a re-solved (`final`) row. `beta_min = 1e-5` on the data scale; the
  certificate is recomputed at the thresholded point. k = free unpenalized
  parameters + penalized |theta| >= 1e-5; V = waves * p.
- **Pruning / step-up**: one pass = psychonetrics `prune()`; step-up =
  `stepup()` non-greedy: MI = n g^2 / (4 V) with V the Schur complement of
  the per-observation information, EPC start, LR test at `alpha`, BIC
  criterion = LR > log n. Default-on because on ESM weekly means one-pass
  pruning removed 22 parameters that were jointly significant (p = 5e-8).
- **Parity**: psychonetrics' default optimizer stops at a gradient of about
  1e-3, so estimates differ by up to ~1e-3; equivalence is asserted as the
  same pruned support and psychnets' F <= psychonetrics' F under
  psychonetrics' own fit function. When the saturated fit is on a boundary
  the two can start the greedy search from different points and select
  different adequate models.
- **Baseline floor**: the baseline model (no temporal or network structure)
  is nested in every fit; a fit ending above it errors
  (`psychnet_no_converge`), except an unrefitted EBIC fit.
- **Old help wording**: the old description said the EBIC lambda is
  "re-solved from three starts"; that is true (path solution, saturated,
  moments). It said p-values reported in tables are unadjusted; still true.
