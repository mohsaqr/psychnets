# psychnets 0.5.0

## Pre-release correctness hardening

* `cor_auto()` now rejects variable pairs with fewer than two joint
  observations instead of optimizing an empty polychoric likelihood.
  Marginal-correlation p-values use pair-specific effective sample sizes;
  full-order partial-correlation inference rejects heterogeneous pairwise
  sample sizes.
* Bootstrap, stability, case-drop, and split-half refits retain the observed
  node set when a resample makes a node constant. Failed-fit counts are exposed,
  invalid small-sample splits are rejected, and `estimator_args` carries
  estimator options whose names collide with diagnostic options. The legacy
  `net_boot(engine=)` argument is translated to the public `native` switch.
* Bootstrap difference tests split exact ties between both tails, so identical
  draw distributions return `p = 1`.
* Corrected undirected node-betweenness normalization, preserved directedness
  in relative-importance aggregation, and prevented rank dichotomization from
  splitting tied values by row order.
* Certificates that are unavailable, including post-threshold GGM graphs and
  conditioned moderated networks, now report `certified = NA`. GGM objects
  retain the pre-threshold residual in `$fit_kkt`.
* The glmnet Ising/MGM paths and moderated MGM now honor `nlambda` and
  `lambda_min_ratio`; moderated MGM honors both `LW` and `HW`. MGM consistently
  stores `$levels`, and prediction accepts original two-level factor columns.
* GGM predictability now uses the fitted marginal variance when the precision
  diagonal is penalized.

Unblocks Nestimate's delegation of its cross-sectional psychometric estimators
to psychnets (see `AGENT-NOTE-NESTIMATE-DELEGATION.md`, gaps A and B). Default
*numeric* output of every estimator is unchanged; all new behaviour is opt-in.
Two functions are **renamed** — see below.

## BREAKING: two reliability verbs renamed

* `casedrop_reliability()` is now **`net_casedrop_reliability()`**.
* `network_reliability()` is now **`net_split_reliability()`**.

The old names are **removed, not deprecated.** Both were also exported by
Nestimate with a *different first argument* (a fitted network there, raw data
here), so `library()` load order silently decided which one a user got and the
failure surfaced as an argument-matching error inside their own code. Nestimate
is the older package and now depends on psychnets, so psychnets is the side that
renames. An erroring alias would not have fixed anything — psychnets would still
own the name and still shadow Nestimate's verb.

Neither the arguments, the return value, nor the numbers changed; only the name.
Update calls with a search and replace. The new names also bring both verbs into
the `net_*` style every other framework verb already uses.

`net_edge_betweenness()` is **not** renamed: both packages give it the same
contract, so it is genuinely one verb rather than a collision.

## `ebic_glasso()` can now serve a resampling caller

* New `lambda_path =` argument: scan a caller-supplied penalty grid, in the
  order given, instead of recomputing one from the input correlation matrix.
  Resampling callers compute the path once from the original matrix and reuse it
  for every resample, so selected penalties stay comparable across draws rather
  than drifting with each sample's largest absolute correlation.

* New `penalize_diagonal =` argument (default `FALSE`, unchanged behaviour).
  `glasso_kkt()` gained a matching argument: with a penalized diagonal the
  optimality condition becomes `W_ii = S_ii + rho`, so grading such a fit with
  the default would report a spurious violation of exactly `rho`. Verified
  against the `glasso` Fortran package to ~1e-11.

* New `refit =` argument. `TRUE` (default) keeps today's behaviour of re-solving
  the selected penalty to machine precision. `FALSE` returns the path fit exactly
  as scanned -- a bit-compatibility switch for consumers holding frozen baselines
  against an unrefitted path scan, not a performance knob. `"unregularized"`
  keeps the selected support and refits the unpenalized graph-restricted MLE on
  it, certified by `ggm_support_kkt()` instead of `glasso_kkt()`. `$support`
  is stored on the object in that case only.

## `mgm_fit()` supports multi-level categorical nodes

* `mgm_fit(native = FALSE)` now models a `factor`/`character` column, or a
  numeric column with more than two levels, as a single k-level categorical node.
  Previously these errored, and the suggested workaround (one-hot encoding) was
  not equivalent: it turns one k-level node into k separate nodes, a different
  model. Matches `mgm::mgm()` to ~1e-15 with exact structure on 3- and 4-level
  factors.

* Node type detection now follows `mgm::mgm()`'s own rule: a `factor`/`character`
  column is categorical, and so is a numeric column with 10 or fewer distinct
  integer values. A numeric column promoted by that integer rule is **warned
  about by name** (condition class `psychnet_type_promotion`), because it
  switches the node from gaussian to categorical -- a different model, not a
  different encoding. Likert and count items hit this rule routinely. No
  previously-working call changes its result: every input newly classified as
  categorical raised an error before, and constant / all-NA columns are still
  dropped as before. A `character` column with more than 10 distinct values is
  refused as an identifier (use `factor()` to insist), and a column declared
  `"c"` via `types =` with more than 10 distinct values is refused as a
  mislabelled continuous variable rather than routed to a multinomial.

* `na_method` applies identically on both engines: the glmnet path is fed the
  same imputed (or listwise-reduced) rows as the base path, and a missing
  multi-level code is imputed with the modal level rather than a mean.

* The certificate of a multi-level node is the multinomial (softmax)
  stationarity residual. Grading each class as an independent logistic
  regression reported ~0.2 for a fit glmnet had converged; the true residual is
  ~1e-5. Related: on the glmnet paths (`mgm_fit` / `ising_fit`,
  `native = FALSE`) a node whose EBIC-selected model is *empty* now certifies as
  optimal (residual ~0) instead of reporting its whole gradient as a violation
  -- the empty solution is exact at glmnet's own `lambda_max`. Only `$kkt`
  changes; weights are untouched.

* The resampling verbs (`net_boot()`, `net_stability()`,
  `net_casedrop_reliability()`, `net_split_reliability()`) now hand each draw
  to the estimator **exactly as the full data would be**: a data.frame keeps
  its factor columns and its incomplete rows, a matrix is coerced without
  dropping rows. So resampling `mgm_fit(native = FALSE)` keeps a factor node
  (previously it was silently dropped), and `na_method = "pairwise"` -- every
  estimator's default -- resamples the rows a direct call would use instead of
  a listwise-reduced subset (previously the verbs listwise-deleted up front,
  so the bootstrapped network was not the network `psychnet()` returns on the
  same data). Node labels are read off the observed fit. For `method = "mgm"`
  node types are decided once on the full data and pinned via `types=` for
  every draw, so a borderline column cannot flip between gaussian and
  categorical from draw to draw, and the promotion warning is raised once per
  verb call. Complete numeric data is byte-identical to before.

* `mgm_fit()`: constant / all-NA columns are dropped before anything is indexed
  against the columns, and `moderators` (a positional index) and a `labels`
  vector given for all original columns now follow that drop; a moderator that
  is itself constant is refused by name. `na_method = "listwise"` is applied
  first, so a node left with one level after deletion is dropped like any
  other constant column on both engines. An explicit `factor` with more than
  10 levels is accepted (that is the documented escape hatch); only `character`
  columns and `types = "c"` declarations on non-factor columns are subject to
  the 10-level guard.

* A 2-level column is recoded through its factor levels, so any binary coding
  (`{1,2}`, `{"no","yes"}`) now works and gives the same fit as `{0,1}`.

* The base kernel (`native = TRUE`) is a binomial solver and still refuses a
  node with more than two levels, now naming the column and pointing at
  `native = FALSE`.

* `net_predict()` refuses an mgm network containing a multi-level categorical
  node: a k-level predictor occupies k-1 design columns, so there is no
  one-coefficient-per-variable matrix to predict from. It errors rather than
  computing from a truncated model.

## Documentation and stability

* Node order is now a documented guarantee, not observed behaviour: `$weights`,
  `$nodes`, and `$precision` follow the input column order and are never sorted.
  Stated in `?ebic_glasso`, `?cor_network`, `?pcor_network`, `?mgm_fit` and
  pinned by a test using non-alphabetical names.

* A **numerical-stability policy** is now documented: estimator numerics are
  stable within a minor version; any change to a selected support, or to a weight
  beyond 1e-10, ships in a minor bump and is listed here. Enforced by
  `tests/testthat/test-golden.R` against frozen fixtures at 1e-12, so downstream
  consumers can pin `psychnets (>= x.y)` against the promise.

* `net_crosswalk()` gained rows for `penalize.diagonal`, `refit`, and
  `lambda_path`, which it previously listed as reference-only.

# psychnets 0.4.3

* First CRAN release.

* `mgm_fit(..., moderators = )` (the moderated mixed graphical model, read with
  `condition()`) now runs on psychnets' own base-R penalized-IRLS kernel and is
  KKT-certified. It previously required the compiled `glmnet` package and
  silently ignored `native =`. Every estimator in the package is now pure base R
  by default; `glmnet` and `glasso` remain optional reference engines selected
  with `native = FALSE`.

* The base moderated kernel covers gaussian and binary nodes (the documented
  `mgm_fit()` scope) and errors explicitly on a categorical node with more than
  two levels, which needs `native = FALSE`.

* The two resampling-heavy long-form guides, "Network reliability" and
  "Visualizing networks", are now website articles rather than installed
  vignettes, so they keep their full iteration counts while the package check
  stays inside CRAN's time budget. Read them at
  <https://pak.dynasite.org/psychnets>. The eight remaining vignettes are
  unchanged and still ship with the package.
