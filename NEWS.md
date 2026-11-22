# psychnets 0.8.6

* `tmfg_certificate()` tests planarity for any graph, by the left-right
  planarity test (de Fraysseix & Rosenstiehl 1985; Brandes 2009), as
  NetworkX's `check_planarity()` does. The earlier stacked-triangulation check
  was valid only for chordal graphs, so a non-planar graph that is also not
  chordal was counted planar: K3,3 now violates 3 properties (was 2), as does
  the bfi TMFG with an edge added. The verdict (0 = a TMFG) is unchanged.
* `ising_fit()` and `mgm_fit()` with the native solver select each node's
  penalty on glmnet's default path, as `IsingFit::IsingFit()` and
  `mgm::mgm()` do: `lambda_min_ratio` now defaults to `NULL`, glmnet's rule
  (1e-4, or 0.01 with fewer observations than predictors), and the path stops
  early by glmnet's rule. The old grid ended at 0.01 of the largest penalty
  and could miss the EBIC choice: on six bfi items the Ising edge weights
  moved by up to 0.0165 and now agree with IsingFit to 1.6e-5; on a mixed
  Gaussian/binary model an edge of 0.0497 that mgm finds was missing, and the
  magnitudes now agree with mgm run with glmnet converged to 6e-10. Pass
  `lambda_min_ratio = 0.01` for a grid that ends at 0.01.

# psychnets 0.8.5

* `as.data.frame()` of a grouped (`psychnet_group`) or multilevel
  (`psychnet_multilevel`) result now passes `...` on to each network, so
  `include_zero = TRUE` returns every node pair of every group or level. It
  was ignored: zero edges were dropped silently (57 rows where 2 x 45 = 90
  pairs exist on the bundled data by gender).

# psychnets 0.8.3

* A nodewise logistic fit that approaches a finite optimum slowly (nearly
  collinear predictors) is no longer reported as separated. 0.8.2 could stop
  such a fit early with a separation warning in `ising_sampler()`.

# psychnets 0.8.2

## New features

* `panel_gvar()` estimates the panel graphical VAR model (Epskamp, 2020) from
  few-wave panel data: a temporal, a contemporaneous and a between-person
  network, by maximum likelihood. Edges are selected with
  `select = c("none", "prune", "ebic")`; `standardize` follows psychonetrics.
  New dataset `srl_panel` and a vignette.
* `psychnet(within = )` estimates within-person networks: each variable is
  centered and scaled within persons, as `mlVAR::mlVAR(scaleWithin = TRUE)`.
  With `group`, `net_compare()` then pairs the occasions by person.
* `net_compare()` accepts two networks of one grouped `psychnet()` result,
  gains `na_method`, and reports the unadjusted p-value (`p_raw`) next to the
  adjusted one.
* `net_stability()` and `net_casedrop_reliability()` gain `cor_method`.
* New `as.data.frame()` methods and `what =` options for stability,
  bootstrap, moderated MGM and panel GVAR results; `summary()` and `plot()`
  for multilevel networks.

## Changes to results

* `net_compare()` stores the network-structure test as `$structure` (M) and
  the global-strength test as `$strength` (S), with the letters of van
  Borkulo et al. (2022); `$M` and `$S` raise a `psychnet_renamed` error.
* `net_stability()` computes the CS-coefficient as
  `bootnet::corStability()`; `net_casedrop_reliability()` now defaults to
  `cor_method = "pearson"`.
* `ebic_glasso()` refits the lambdas near the EBIC minimum before choosing
  one, and returns absent edges as exact zeros.
* `mgm_fit()` solves gaussian nodes to optimality (edge weights move by up to
  about 1e-7).
* `mgm_fit()`, `ising_fit()` and `ising_sampler()` rescale `weights` to a
  maximum of 1, as `mgm::mgm()` does.
* `net_bridge(normalize = TRUE)` divides the 2-step bridge expected influence
  by the number of nodes, as `networktools::bridge()`.
* `net_split_reliability()` gives a split with one empty half a correlation
  of 0.
* `condition()` reads `value` in the moderator's own coding.
* `ising_sampler()` warns under (quasi-)complete separation, does not certify
  such a fit, and stops it early.

## Input and conditions

* Errors and warnings are classed conditions (`psychnet_error`,
  `psychnet_bad_input`, ...), and arguments are validated.
* Ordered factors are kept as ordinal items; `ising_fit()` codes logical and
  two-level factor columns 0/1.

## Plots

* All plots use the Okabe-Ito palette and show significance by fill and line
  type as well as colour.

# psychnets 0.5.2

## CRAN incoming-check runtime

* Resubmission addressing the automated incoming pre-test NOTE (12 examples
  over 10 s and an overall checktime over 10 min on the Windows builder). No
  estimator numerics changed.
* Resampling examples (`net_boot()`, `difference_test()`, `net_stability()`,
  `net_casedrop_reliability()`, `net_split_reliability()`, and their plot
  methods) now demonstrate the workflow with the closed-form
  `method = "pcor"` estimator and note that the `"glasso"` default with full
  iteration counts is what a real analysis should use. `net_compare()`
  examples use fewer permutations on smaller data.
* The cograph visualization vignette's bootstrap chunk uses `method = "pcor"`.
* Four heavyweight numerical test files are skipped on CRAN (they run on every
  GitHub Actions matrix job, R-hub platform, and local `devtools::test()`);
  the frozen 1e-12 golden estimator baselines still run on CRAN.

# psychnets 0.5.1

## Pre-release correctness hardening

* The base nodewise kernel (`ising_fit()`, `mgm_fit()`, moderated MGM with
  `native = TRUE`) now returns exact zeros for inactive lasso coordinates.
  Coordinate descent could leave a BLAS-rounding crumb (~1e-17) on a
  coordinate whose gradient sits exactly at the penalty boundary; whether the
  crumb was `0.0` or `4e-17` differed between BLAS implementations, and the
  AND rule then counted it as an edge on some platforms but not others. The
  EBIC degrees of freedom already ignored such coefficients (`|beta| <=
  1e-10`); the returned support now agrees. One golden-fixture edge changes:
  the base-kernel mixed-data mgm network loses a spurious weak edge, matching
  the glmnet reference path, which always reported that edge as absent.
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
