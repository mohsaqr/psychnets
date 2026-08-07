# psychnets 0.5.0

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
  it, certified by `ggm_support_kkt()` instead of `glasso_kkt()`.

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
  about by name**, because it switches the node from gaussian to categorical --
  a different model, not a different encoding. Likert and count items hit this
  rule routinely. No previously-working call changes its result: every input
  newly classified as categorical raised an error before.

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
