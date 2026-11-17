# Implementation notes removed from `?condition` (2026-11-17)

Still true of R/mgm_moderated.R; kept for maintainers.

- **Mirrors `mgm::condition()`**: applies the AND-rule pre-filter
  (`pw_alive`, `int_alive` computed at fit time with the fit's rule; with
  `rule = "OR"` everything is alive), absorbs the moderator value into the
  main-effect coefficients (`.mmg_cond_indicator()`), and re-aggregates the
  pairwise edges as the mean absolute value over dummy columns (and over
  classes for a multinomial response), then AND/OR over the two sides.
- **Moderator scale**: a categorical moderator's value is matched against
  its original level labels (`moderator_levels`) and mapped to the internal
  0-based level code; a continuous moderator is mapped to
  `(value - mean) / sd` (mgm's `scale = TRUE`). `mgm::condition()` takes the
  standardized value: `condition(fit, v)` here equals
  `mgm::condition(ref, list(k = (v - mean) / sd))`. Model-matrix columns are
  named `V<k>.<code>`.
- **Interaction design**: each non-moderator node is regressed on all main
  effects plus `predictor * moderator` terms; the moderator node on all main
  effects plus every pairwise product of the other variables. Interaction
  columns whose variables are all gaussian are standardized after the model
  matrix is built (mgm main, lines 216-226).
- **Base engine**: a binary node is fitted by the binomial lasso and its
  coefficient split into glmnet's symmetric 2-class form `(-b/2, +b/2)`; the
  EBIC on the base path selects on the same objective as the glmnet path
  (see the comment above `.mmg_fit_node_base()`).
- **Certificate**: the moderated fit stores `kkt` / `fit_kkt` /
  `certificate_target` (base engine only; `NA` with glmnet, and `NA` after a
  threshold that removed a coefficient), but it is not a `psychnet`, so
  `certificate()` refuses it with `psychnet_bad_input`; there is currently
  no accessor for it. The network returned by `condition()` (method
  `"mgm_moderated"`) stores no `kkt`, so `certificate()` reports `NA`.
