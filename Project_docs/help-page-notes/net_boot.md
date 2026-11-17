# Implementation notes removed from `?net_boot` (2026-11-17)

Still true of R/bootstrap.R, R/reliability.R, R/correlation.R and R/plot.R;
kept for maintainers.

- **Provenance**: adapted from Nestimate's `boot_glasso` into a general
  resample-and-re-estimate loop over any estimator routed by `psychnet()`
  (`.psn_refit_network()`).
- **cores** (`.resolve_cores()`): `NULL` -> `floor(2/3 * detectCores())`,
  at least 1; an explicit value is capped at the detected count; capped at 2
  when `_R_CHECK_LIMIT_CORES_` is set (CRAN); forced to 1 on Windows, where
  `parallel::mclapply()` cannot fork. `detectCores()` returning `NA` counts as 1.
- **Core-count invariance**: all resample indices are drawn in the parent
  (`sample.int(n, n, replace = TRUE)`, one per resample) before fitting, and
  each draw runs under `.psn_try_draw()` inside its worker, so a failing draw
  costs only itself (with forked workers an uncaught error would poison every
  draw prescheduled to the same worker). The result is therefore identical for
  any number of cores.
- **Failed draws**: recorded as `NA` rows in the stored draw matrices;
  `.psn_check_draws()` warns once (`psychnet_resample_failed`, first error
  message) and errors (`psychnet_resample_failed` + `psychnet_estimation_failed`)
  when fewer than two draws succeed. Intervals and means use `na.rm = TRUE`
  (default `quantile()` type 7). The column-drop notices are muted inside a
  draw (`.psn_with_quiet_columns()`).
- **engine**: translated by `.resample_engine()` to the estimators' `native`
  switch (`native = identical(engine, "base")`); allowed values are
  `base`/`glasso` for glasso and huge, `base`/`glmnet` for ising and mgm, and
  nothing for the other methods.
- **prop_nonzero**: share of non-failed draws with `abs(weight) > 1e-12`.
- **Stored draws**: fields `edge_boot` (n_boot x edges), `centrality_boot`
  (list by measure), legacy aliases `str_boot` / `ei_boot` for strength and
  expected influence, `predictability_boot`, plus `edge_labels`,
  `node_labels`, `lambda_path`, `lambda_selected`, `n_success`, `n_failed`.
- **diff_test**: stores `edge_diff_p` (`NULL` past 500 edges,
  `.psn_diff_pmat(max_items = 500)`) and `centrality_diff_p` (list by measure),
  p-values adjusted over the unique pairs by `p_adjust`. No accessor returns
  them. (The old help said these matrices are "drawn by `plot(x, type =
  "edge_diff")` / `"centrality_diff"`"; the code's plot method recomputes the
  tests with `difference_test()` from the stored draws, without `p_adjust`,
  whether or not `diff_test = TRUE`. The new page no longer ties them to the
  plot.)
- **Predictability**: `.psn_predictability()` uses `diag(solve(K))` as the
  variance when the diagonal was penalized, else the diagonal of the stored
  correlation matrix, else 1; values clamped to \[0, 1\]. A node absent from a
  draw has `NA`.
- **mgm types**: `.resample_input()` decides node types once on the full data
  (`.detect_types()`) so a borderline column (e.g. 11 distinct integers in
  full, <= 10 in a draw) cannot change type between draws.
- **Group mode**: `.group_data_apply()` rebuilds `estimator_args` from the
  group's stored call (dots, gamma, threshold where the method takes it) and
  returns class `c("psychnet_bootstrap_group", "psychnet_result_group")`.
