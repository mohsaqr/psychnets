# Implementation notes removed from `?mgm_fit` (2026-11-17)

Still true of R/mgm.R, R/lasso_glm.R, R/lasso_glmnet.R and
R/mgm_moderated.R; kept for maintainers.

- **Provenance wording removed**: "pure base R", "dependency-free",
  "self-certified", nodewise KKT about 1e-9.
- **native = FALSE parity**: the old help said the edge magnitudes
  byte-match `abs(mgm::mgm()$pairwise$wadj)` "to ~1e-6". The glmnet path
  (`.mgm_fit_glmnet()`) scales gaussian columns once up front (as
  `mgm::mgm` does; per-node scaling diverges), expands categorical nodes
  with `model.matrix()`, fits gaussian lasso / multinomial lasso (a binary
  node is a 2-class multinomial), EBIC with `npar` = number of DUMMY
  columns, LW/HW thresholding, AND/OR on magnitudes. Edge sign is recovered
  from a gaussian endpoint where one exists; a binary-binary sign is
  undefined there and stored positive. Certificate: `.self_lambda_kkt()` /
  `.self_lambda_kkt_multinomial()` (glmnet's tolerance; 2.6e-6 on the help
  example).
- **Base path scaling**: continuous columns standardized with the sample
  SD (frequency-weighted, divisor `sum(w) - 1`, with `weights`); the
  composite raw -> standardized-predictor transform is stored for
  `net_predict()` (`nodewise$center/scale`, `resp_center/resp_scale`).
- **Binary response scale**: the logistic coefficient is halved
  (`beta_graph / 2`) to put it on mgm's two-class multinomial scale before
  thresholding and symmetrization.
- **LW threshold** (base path, unmoderated): `tau_i = ||B[i, -i]||_2 *
  sqrt(log(p - 1) / n)`; HW: `sqrt(log(p - 1) / n)`. The moderated path
  (`.mmg_threshold()`) uses an extra factor `sqrt(2)` in LW and `npar` =
  number of design columns (main effects plus interactions).
- **Certificate bookkeeping**: `kkt = NA` and `certificate_target =
  "pre_threshold_fit"` when the threshold zeroed any selected coefficient;
  `fit_kkt` keeps the residual of the penalized fits. Same rule on both
  engines.
- **Constant columns**: dropped up front by `.usable_columns()`, silently
  (no `psychnet_constant_variable` warning, unlike `ising_fit()`). The old
  help said "dropped, as in every other estimator".
- **Listwise deletion** happens before type detection, so levels seen only
  in incomplete rows do not count.
- **Unreachable check**: the "declared 'c' but not coded 0/1" error on the
  base path cannot fire, because categorical columns are recoded to 0-based
  level codes before the check.
- **Multi-level categorical**: base kernel is binomial only and errors by
  name; `net_predict()` refuses a network with a k-level node (k-1 design
  columns, no one-coefficient-per-variable matrix).
