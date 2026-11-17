# Implementation notes removed from `?ising_fit` (2026-11-17)

Still true of R/ising.R, R/lasso_glm.R and R/lasso_glmnet.R; kept for
maintainers.

- **Provenance wording removed**: the old help described the default solver
  as "pure base R", "dependency-free" and "self-certified", with a nodewise
  KKT of about 1e-9 (the example gives 1.45e-09). The new page says only that
  each regression is checked against its optimality conditions.
- **Default engine** (`native = TRUE`): `.nodewise_ebic()` (penalized IRLS
  with inner coordinate descent, `.wls_lasso()`), per-node lambda path of
  `nlambda` values down to `lambda_min_ratio`, EBIC selection. Predictors are
  standardized with (weighted) centers/scales from `.standardize(X,
  weights)`; raw-scale slopes are stored in the edge matrix and the
  standardized ones in `nodewise$beta_std`.
- **Thresholds**: raw-scale node thresholds are the standardized intercept
  de-centered by the raw slopes, `tau_i = b0_std_i - sum_j beta_raw_ij *
  center_j`, so they are coherent with the raw-scale edge weights.
- **Certificate**: `kkt` = max over nodes of the nodewise KKT residual. The
  AND/OR combination is not treated as a post-fit change: the certificate
  stays defined.
- **native = FALSE** (`.ising_fit_glmnet()`): glmnet per node on raw 0/1
  predictors with the IsingFit EBIC path, AND/OR symmetrization, thresholds =
  glmnet `a0`. The old help claimed `weights`/`thresholds` byte-match
  `IsingFit::IsingFit()` "to ~1e-16". The certificate is
  `.self_lambda_kkt()`: the effective penalty is recovered as the mean active
  |gradient| on glmnet's standardized scale (max |gradient| for an empty
  model), so it reflects glmnet's tolerance (3.3e-6 on the `SRL_GPT`
  example, above `certificate()`'s default tolerance 1e-6). No convergence
  warning is raised on this path. `nodewise$center/scale` are 0/1 there and
  `beta_std` holds the raw coefficients.
- **lambda_min_ratio with native = FALSE**: when the argument is missing,
  `NULL` is passed so glmnet's own default applies (IsingFit parity); an
  explicit value is honoured.
- **weights**: rescaled to max 1 (`mgm::mgm()` convention) after `min_sum`
  filtering, so a removed row cannot change the scale; the EBIC uses
  `log(sum(weights))`, which is why a fixed scale is needed.
- **Pairwise NA fill** (`.na_prep_nodewise()`): a 0/1 column is filled with
  `as.numeric(mean(obs) >= 0.5)`, i.e. the modal value with ties to 1.
- **Factor coding**: a two-level factor's first level is 0; the training
  coding (`category_levels`) is kept so `net_predict()` encodes new data by
  label.
