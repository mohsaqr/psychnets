# psychnets 0.7.1

* Tests: the frozen network-comparison baseline is compared to a tolerance of
  1e-10 instead of exact identity. It was frozen on macOS, and other
  platforms' floating-point arithmetic differs in the last bits (~1e-15), which
  failed the check on Linux and Windows. No code or returned number changed.

# psychnets 0.7.0

This release follows a full function-by-function review of the package. Most
changes are input validation and conditions that were silent before; the ones
marked **(numbers change)** alter returned values, as the numerical-stability
policy requires for a minor version.

## Changes to returned numbers

* **(numbers change)** `net_stability()` computes the CS-coefficient with the
  rule of `bootnet::corStability()`, the source of the 0.25/0.5 guidance:
  Pearson correlations (`cor0`: 0 when one side is constant) and strict `>`
  for both `threshold` and `certainty`. New `cor_method = c("pearson",
  "spearman")`; Spearman was the only option before. Undefined correlations
  still count as failing (bootnet drops them). `net_casedrop_reliability()`
  uses the same rule and its `cor_method` default changes from `"spearman"` to
  `"pearson"`. `threshold` and `certainty` must lie in (0, 1).
* **(numbers change)** `mgm_fit()`: gaussian nodes are solved to optimality
  (coordinate descent on the Gram matrix plus an exact active-set solve)
  instead of stopping after 100 sweeps. On collinear items the residual falls
  from about 5e-5 (uncertified) to about 1e-15, and the fit is faster. Default
  edge weights move by up to about 1e-7 (3.4e-8 on the golden data, no edge
  added or removed); the golden fixture was regenerated. A fit that still
  misses its optimum warns `psychnet_no_converge`.
* **(numbers change for weighted fits)** `mgm_fit()`, `ising_fit()` and
  `ising_sampler()` rescale `weights` to a maximum of 1, as `mgm::mgm()`
  does, so only relative weights matter (weights summing to 1 used to return
  an empty network). With `min_sum`, the rescaling uses the retained rows only.
* **(numbers change)** `net_bridge(normalize = TRUE)`: the 2-step bridge
  expected influence is normalized by the number of nodes, as in
  `networktools::bridge()`, not the number of communities.
* `condition()` reads `value` in the moderator's own coding: a level of a
  categorical moderator (1/2, `"treat"`, ...) or the raw scale of a continuous
  one. It used an internal 0-based code or a z value, and a value that is not
  a level silently returned the reference group; that is now an error.
* `panel_gvar()`: character waves that are all numbers are ordered
  numerically ("10" used to sort before "2", fitting the waves out of time
  order); other character waves warn `psychnet_wave_order`. The Fisher-scoring
  polish halves its step instead of stopping, removing spurious
  `psychnet_no_converge` warnings.
* `net_boot()` and the case-dropping verbs catch a failing resample inside its
  own draw, so results are identical for any `cores` (one failing draw used to
  void a whole forked chunk).

## Certificates

* `mgm_fit()`: when the LW/HW threshold removes a coefficient, `$kkt` is `NA`
  with the pre-threshold residual in `$fit_kkt`, as for every other
  post-thresholded estimator (it reported the pre-threshold residual as
  certified).
* `ising_sampler()` warns `psychnet_separation` under (quasi-)complete
  separation and no longer certifies such a fit. Separation is detected by a
  drift test (the objective still rising along the direction the fit moves),
  so an extreme but finite fit is still certified.
* A binary node left constant by `min_sum` filtering or zero weights is
  degenerate (no finite intercept): `psychnet_degenerate_sample` warning and no
  certificate, instead of a certified fit with a huge threshold.
* `lmg_certificate()` also checks each node's full-model R^2 against its
  closed form (the efficiency identity alone holds for any set function).
* `logo_network()` refuses ill-conditioned clique blocks (`psychnet_singular`)
  instead of returning a network that fails its own certificate.
* `ggm_modselect()` reports failed candidate fits (`psychnet_fit_failed`) and
  uncertified results (`psychnet_no_converge`); `ggm_support_kkt()`
  symmetrizes a triangular support; `print()` labels the TMFG and relimp
  certificates as structural.

## Input handling

* Every input-contract error is a classed condition
  (`psychnet_bad_input`, subclass of `psychnet_error`), and `threshold`,
  `alpha`, `n`, `labels`, `gamma`, `nlambda`, `lambda_min_ratio`, `ci`,
  `cores`, `iter`, `drop_prop`, `moderators` and the certificate helpers'
  arguments are validated. `threshold = TRUE` (qgraph's logical form) and
  other non-numeric thresholds used to return an empty network.
* Data frames: ordered factors are kept as ordinal items (so
  `cor_method = "auto"` treats them like `qgraph::cor_auto()`); other
  non-numeric columns are dropped with a `psychnet_dropped_columns` message and
  constant columns with a `psychnet_constant_variable` warning. `ising_fit()`
  codes logical and two-level factor columns 0/1 instead of dropping them, and
  `net_predict()` encodes new data with the same training coding.
* Significance thresholding with too few observations (df <= 0) errors
  (`psychnet_degenerate_sample`) instead of deleting every edge;
  `pcor_network()` warns `psychnet_singular` for a singular correlation matrix.
* `psychnet()` warns `psychnet_ignored_argument` when `threshold` or `gamma`
  is given to a method that does not use it, and refuses several methods.
* Resampling verbs report the estimator's own error (`psychnet_estimation_failed`)
  instead of "estimation failed", and a failed draw's message
  (`psychnet_resample_failed`).
* `net_compare()`: the edge table gains `p_raw`, the unadjusted permutation
  p-value, next to `p_value` (adjusted by `p_adjust`; `significant` uses it).
  Previously the unadjusted values were discarded once an adjustment was
  requested. The global M and S tests are never adjusted.
* `net_compare()` refuses infinite values, lets a `gamma`/`na_method` you
  supply override a group's (warning `psychnet_estimator_mismatch`), and
  labels group-mode edges with the group's `labels`.
* Group mode follows factor level order, refuses empty-string levels and
  levels with fewer than 3 rows, and warns for small levels and for rows with
  a missing group. Event data refuse missing actors or sessions.
* `net_centralities(centrality_fn =)` aligns a named result by node label;
  `net_edge_betweenness()` treats an asymmetric matrix as directed;
  `net_bridge()`/`net_aggregate()` refuse duplicated or unknown community
  labels; `net_aggregate()`/`redundancy()` report listwise-dropped rows.
  `net_aggregate()` leaves out a column whose community is `NA` (e.g. an ID)
  before any preprocessing, and orders numeric communities numerically
  (`"2"` before `"10"`); character communities stay alphabetical.
* `relimp_network()` gains an ignored `n` argument (as `tmfg_network()`).

## New verbs

* `scale_within()` centers each variable within persons and, by default,
  divides by the person's own SD (the within-person standardization of
  `mlVAR::mlVAR(scaleWithin = TRUE)`, which it reproduces exactly), as a data
  step before any analysis -- e.g. a paired `net_compare()` of two occasions'
  within-person networks. It warns (`psychnet_mirror_deviations`) when no
  person has more than two observations, where the occasions' deviations are
  mirror images and their networks identical by construction, and sets values
  without within-person variation to `NA` (`psychnet_constant_within`).
* `psychnet(within = "person")` applies the same scaling inside the front door
  (over all rows, before `group` splits them), drops the person column from
  the nodes, and records each row's person, so
  `net_compare(nets, c(1, 4))` compares two occasions paired by person, by
  default (independent of row order; persons seen in one level only are left
  out with a `psychnet_rows_dropped` warning; `paired = FALSE` overrides).
  Other groups keep the independent-groups default.

## New accessors

* `as.data.frame()` for `net_stability()` results (`what = c("table", "cs",
  "correlations")`) and for moderated MGM fits; `as.data.frame()` on a
  bootstrap gains `what = c("edges", "centrality", "predictability",
  "thresholded")`; on a panel GVAR, `what = "path"` (EBIC lambda path) and
  `what = "covariance"`. `psychnet_multilevel` gains `summary()` and `plot()`.
* **(numbers change)** `net_split_reliability()`: a split in which one half
  network is empty and the other is not has edge correlation 0 (bootnet's
  `cor0()`); such splits used to be dropped, inflating the mean correlation.
  Only splits where both halves are constant are undefined.
* `net_split_reliability()` gains `n_valid` and reports NA rather than NaN;
  stability and case-drop tables gain `actual_drop` and `n_valid`.

## Plots

* All plots use the Okabe-Ito palette. Significance is shown by filled vs open
  points and solid vs dashed segments, series by line type and point shape as
  well as colour, box-matrix cells print the sign, and the bridge plot styles
  nodes by community. Every plot restores `par()`.

## Documentation

* S3 methods are documented on the page of the verb that creates the object
  (`?net_boot` covers its `print`, `plot` and `as.data.frame` methods); the
  `psychnet` object has its own page, `?"psychnet-class"`. The manual shrinks
  from 69 to 48 topics.
* `cor_auto()` documents that its polyserial is the Olsson ad-hoc estimator;
  `net_edge_betweenness()` returns path counts; `net_centralities()`
  betweenness and closeness are normalized (as igraph, unlike
  `qgraph::centrality_auto()`); the crosswalk maps qgraph's
  `EBICglasso(refit = TRUE)` to `refit = "unregularized"`.

# psychnets 0.6.2

## Documentation

* Vignettes rechecked against current output. The Ising and mixed graphical
  model vignettes described networks from an earlier release (7 edges and a
  self-regulation to test-anxiety edge); their prose now matches the current
  fits (6 edges, test anxiety isolated under the AND rule; the OR rule adds the
  negative edge). Each vignette shows the edge table with `as.data.frame()`,
  which `summary()` returns only invisibly, and lists its method references.
* Examples no longer reach into result objects with `$`: certificates are read
  with `certificate()`, stability with `print()`, and the certificate helpers
  (`glasso_kkt()`, `ggm_support_kkt()`, `glm_lasso_kkt()`) are shown on exact
  closed-form optima. `panel_gvar()` examples use the bundled `srl_panel`; Ising
  examples use `dichotomize(SRL_GPT)`. Every print method has an example.
* Method references added to the estimator, centrality, predictability and
  resampling help pages; DESCRIPTION references now give authors and years.

# psychnets 0.6.1

## Complete function-review fixes

* Case-dropping stability now counts undefined correlations and failed fits
  against the requested certainty, rather than inflating stability by omitting
  those draws. Two-node edge-reliability diagnostics return available deviation
  metrics with an unavailable correlation instead of failing.
* Resampling preserves source-column identities when custom labels and ignored
  nonnumeric columns are present. Grouped diagnostics forward explicit estimator
  options, composite grouping keys cannot pool distinct tuples, and groups of
  within/between models can be printed and summarized.
* Event-feature missingness follows the requested handling instead of deleting
  an otherwise varying feature. Within-person centering uses available values
  per variable. Generic independent-row resampling and NCT now refuse repeated-
  person event groups: reconstructing these as cross-sectional samples had lost
  their effective sample size and could change even the observed graph.
* Panel models refuse missing person identifiers, support empty EBIC networks
  in all accessors and plots, return CFI = 1 when both estimated noncentralities
  are zero, and match edge significance without ambiguous label concatenation.
  Plot abbreviations override defaults, custom line styles no longer carry a
  misleading significance key, and `nlambda` must be a finite whole integer.
* Kendall marginal-network significance uses a tie-corrected Kendall test.
  Kendall significance without raw data, or for full-order partial
  correlations, is refused rather than applying Pearson p-values. Polyserial
  estimates are invariant to the numeric spacing of ordinal codes.
* MGM automatic predictability uses raw-scale Gaussian observations. Public
  prediction preserves original source columns independently of display labels,
  applies the training category encoding to numeric and factor inputs alike,
  and refuses unseen categories.
* Wald-pruned IsingSampler fits retain their pre-pruning score as `fit_kkt` and
  report the altered network's certificate as unavailable. Covariance PSD
  validation runs in correlation units, including for relative-importance fits.
* Two-node external Ising/MGM fits work through an inert design-padding column;
  moderation requires at least three nodes. TMFG certificates now also validate
  planar stacked-triangulation structure.
* Small-world null rewiring randomizes undirected edge orientations and retains
  rejected proposals, removing node-order restrictions and accepted-move bias.
  NCT plots identify signed global expected influence when `abs = FALSE`.

## Review fixes

* Panel EBIC path diagnostics now recompute the objective and KKT residual
  after small coefficients are rounded to zero. Infeasible rounded candidates
  are excluded from selection instead of retaining a stale certificate.
* `ebic_glasso(refit = "unregularized")` now uses the same confirmed lambda
  and support as the default penalized fit before the unpenalized refit.
* Glasso's returned `ebic_path` includes the rescored EBIC values, so its
  minimum agrees with the selected lambda and reported EBIC.
* Panel plots accept overrides such as `node_fill`, `layout`, `title`, and
  edge colours through `...` without duplicate-argument errors.

## New estimator: panel graphical VAR

* `panel_gvar()` estimates the panel graphical vector-autoregressive model
  (Epskamp 2020) from few-wave, many-person panel data: a directed temporal
  network (within-person lag-1 effects), a contemporaneous network (partial
  correlations of the within-person innovations), and a between-person network
  (partial correlations of the stable person means). The within-person process
  is a stationary VAR(1) around a random intercept with stationary means -- the
  model of `psychonetrics::panelgvar()` -- fitted by maximum likelihood in base
  R (nlminb with an analytic reverse-mode gradient, then a Fisher-scoring
  polish). Long data (`vars`, `id`, `wave`) or wide data (`vars` as a
  variables-by-waves design matrix) are accepted.
* Robust on raw-scale, small-sample data: the optimizer works on standardized
  variables (an exact reparameterization -- the model is scale-equivariant)
  from three deterministic starts (moment-based, neutral, and the fitted
  baseline), keeps the best, and refuses (`psychnet_no_converge` error) to
  return a fit worse than its nested baseline. On weekly means of an
  experience-sampling panel (41 persons, 0-100 scales) this reaches a lower
  objective than psychonetrics on every design compared, where a single
  unscaled start had diverged.
* `plot()` draws the three networks side by side on one shared circle layout,
  with the temporal autoregressive effects as self-loops, edge weight labels,
  solid (significant) vs dashed (not) edges, blue positive / vermillion
  negative (Okabe-Ito), and "no edges" on an empty network.
* `as.data.frame(fit, what = "variance")`: per-variable split of the implied
  variance into the between-person and stationary within-person parts (an
  intraclass correlation, `share_between`).
* New dataset `srl_panel`: a synthetic panel of 100 students x 4 weeks on five
  self-regulated-learning indicators (0-100), simulated once from a
  panel-GVAR whose structure was taken from real experience-sampling data
  (the anonymized `esm_srl` data of the idiographic package). New vignette
  "Panel graphical VAR: temporal, contemporaneous, and between-person
  networks" works through it.
* `standardize = c("none", "z", "z_per_wave")`, as psychonetrics'
  `panelvar(standardize =)`: `"z_per_wave"` z-scores every variable within
  every wave, removing drift in wave means and variances (matches
  psychonetrics' chi-square to three decimals).
* Penalized EBIC selection, `select = "ebic"` -- the panel counterpart of
  `ebic_glasso()` and psychonetrics' `panelgvar(estimator = "PML")`: lasso on
  every temporal (autoregressions included), contemporaneous and between
  parameter, 50-value lambda path from lambda_max with warm starts, EBIC
  (`gamma = 0.5`), optional unpenalized `refit`. Accelerated proximal gradient
  that never leaves the model's domain, then an active-set Newton polish;
  certified by the lasso KKT conditions (~1e-10); the selected lambda is
  re-solved from three starts and their agreement reported. Matches
  psychonetrics' PML on every simulated design compared (same lambda, same
  zero pattern, estimates within 3e-3, lower penalized objective); on weekly
  ESM means psychonetrics' PML returns an invalid model (non-positive-definite
  contemporaneous block) while `select = "ebic"` returns a certified one.
  Lasso estimates carry no p-values (`significant` is `NA`; the plot draws
  selected edges solid).
* Robustness (from an adversarial review): the EBIC choice is reconsidered
  after the final certified solve (a final solve can change a lambda's EBIC;
  the path table marks fully solved rows in `final`); only variances on the
  optimizer's actual bound are treated as active constraints (a near-zero but
  interior variance no longer has its score suppressed in the certificate);
  singular-information SEs use the pseudo-inverse on the identifiable
  eigenspace (dropping unidentified coordinates had understated the SEs of
  parameters conditioned on their identifiable combinations); all post-fit
  computation runs on the standardized scale, so data whose variables differ
  in units by orders of magnitude (or are tiny) fit identically instead of
  failing a covariance check or crashing in a data-scale solve.
* Edge selection is one argument, `select = c("none", "prune", "ebic")`
  (replacing the earlier logical `prune`).
* Default output is the full model: `prune = FALSE` by default (the
  psychonetrics default, where `prune()` is a separate step). Every edge is
  reported with its SE, p-value, Wald CI and a new `significant` column
  (p < `alpha`); `print()` and `summary()` count significant edges per network,
  and `plot()` draws significant edges solid and non-significant ones dashed
  (sign by colour plus the signed weight label, with a line-type key).
  `select = "prune"` runs prune + step-up.
* Step-up search after pruning (`stepup = TRUE`, the default when pruning;
  `criterion = c("bic", "none")`), equal to `psychonetrics::stepup()`: re-add
  the pruned edge with the largest modification index while it is significant,
  passes a chi-square difference test, and lowers BIC. One-pass pruning alone
  can remove parameters that substitute for each other (an autoregression and
  correlated cross-lagged effects): on weekly ESM means it removed 22
  parameters that were jointly significant (chi-square 76.9, p = 5e-8), and
  step-up restores the two that matter. Modification indices match
  psychonetrics' (score test, `n g^2 / (4 V)` with the Schur-complement
  information); the selected supports match psychonetrics' `prune()` and
  `prune() |> stepup()` on every well-identified design compared.
  `stepup = FALSE` reproduces `prune()` alone.
* Variances estimated at zero (Heywood cases) are treated as active bound
  constraints: left out of the information (their SE is NA, the others keep
  valid SEs) and certified by the bound KKT condition. A fit whose information
  is singular is pruned with infinite standard errors for the parameters it
  does not identify (rather than a pseudo-inverse, which would understate
  them). Boundary warnings now name the part on (or heading for) the boundary,
  e.g. "between-person variance of V2 at zero".
* `as.data.frame(fit, what = "diagnostics")`: one row with persons, waves, the
  ML score certificate, boundary flag, SE availability, and how many of the
  three starts reached the best optimum (fewer than three = competing local
  optima).
* Pruning follows `psychonetrics::prune()`: Wald p-values from the expected
  Fisher information, `alpha = 0.01`, optional `adjust` (any
  `p.adjust` method) and `recursive`; `prune = FALSE` returns the saturated fit.
* Self-certified: the ML score on the free parameters is recomputed by a
  forward-mode Jacobian (independent of the optimizer's reverse-mode gradient)
  and reported by `print()`; it is ~1e-13 on the validation designs.
* Validated against psychonetrics 0.17.8 on three designs (3 x 4, 5 x 5, and a
  BH-adjusted n = 150 case): identical pruned support in every case, and
  psychonetrics' own fit function rates the psychnets estimate as the deeper
  optimum (its default optimizer stops with a gradient near 1e-3, which moves
  estimates by up to ~1e-3). Fit indices (chi-square, CFI, TLI, RMSEA, AIC,
  BIC) agree to that tolerance.
* Read the fit with `as.data.frame(fit)` (edges with SE, p, Wald CI),
  `as.data.frame(fit, what = "parameters")`, `as.data.frame(fit, what =
  "fit")`, or `as.data.frame(fit, what = "diagnostics")`. The object is a `netobject_group`, so `plot()` /
  `cograph::splot()` draws the three networks as a grid.
* Classed conditions: `psychnet_bad_input` errors (fewer than three waves, too
  few persons, `na_method = "pairwise"`, ...), and `psychnet_boundary`,
  `psychnet_singular_information`, `psychnet_listwise_drop`,
  `psychnet_no_converge` warnings. A near-singular contemporaneous or between
  covariance (partial correlations near +/-1, common with three waves) is
  flagged rather than reported silently.
* Missing data: listwise only. Full-information ML (psychonetrics' default
  when data have missing values) is not yet implemented, and a pairwise
  covariance is refused because it is not a likelihood for this model.
* `net_crosswalk("panelgvar")` maps the psychonetrics arguments.
* No existing estimator's numbers changed; the new golden baseline lives in a
  separate fixture (`golden-panelgvar.rds`).

## `ebic_glasso()`: absent edges are exact zeros

* An edge the EBIC selection counted as absent could be returned in
  `$weights` as a floating-point residue of about 1e-17 instead of 0 (the
  `$edges` table already dropped it). This happened when the selected lambda
  sat at max|S_ij| (the empty graph) after the positive-definite projection,
  and in a substantial share of sparse fits generally (10 of 40 random
  n = 60, p = 12 null draws in the new regression test). Anything that compares
  weights with zero read the residue as an edge; in `net_compare()` two empty
  networks gave M = S = 5.6e-17 and a false p = 0.02 instead of 1 (found by
  the NetworkComparisonTest parity grid). Precision entries with
  |theta_ij| <= 1e-10 -- the same rule the EBIC uses to count edges -- are now
  returned as exact zeros, so the network has exactly the support the EBIC
  scored. Every edge above 1e-10 is unchanged; the golden baselines pass at
  1e-12. Because the returned support can change (residue -> 0), this ships in
  the 0.6.0 minor release per the numerical-stability policy.

## `ebic_glasso()`: the selected lambda is the exact EBIC minimum

* The lambda path is scanned at a loose tolerance (`1e-4`) and only the chosen
  lambda was refit tightly. The scan can leave an entry that is zero at the
  optimum just above the 1e-10 edge-count threshold; one miscounted edge
  raises that lambda's EBIC by exactly log(n) + 4 gamma log(p) and can hide the
  true minimum. The NetworkComparisonTest parity audit found this in 4 of 29
  divergent fits (EBIC up to 9.6 units worse than the optimum). Every lambda
  whose scan EBIC is within two edge penalties of the best is now refit at the
  refit tolerance and rescored, and the lowest exact EBIC wins. On the audited
  fits psychnets now reaches the lower exact EBIC in 29 of 29 (NCT/qgraph in
  4, because qgraph counts a ~1e-17 residue as an edge at lambda = max|S_ij|
  and so never selects an empty graph). Default engine only; `engine =
  "glasso"` is unchanged. About 40% more time per fit; golden baselines,
  equivalence and every estimator suite pass unchanged.

## Network comparison test (`net_compare()`) fixes

* **BREAKING -- the M and S letters were the wrong way round.** In van
  Borkulo et al. (2022) and `NetworkComparisonTest`, M is the
  network-structure test (largest absolute edge difference) and S the
  global-strength test; psychnets printed and stored them swapped, so a write-up
  copied from `print()` reported the wrong statistic. The tests are now stored
  as `$structure` (M) and `$strength` (S); reading the old `$M` / `$S` raises a
  `psychnet_renamed` error that says which field to use, rather than silently
  returning the other test. Printed and tabled labels use the NCT letters.
  `$strength$by_group` (shown by `print()`) holds each network's global
  strength, NCT's `glstrinv.sep`; with `abs = FALSE` the test
  is labelled global expected influence, as in NCT.

On complete data every returned number is unchanged (frozen in
`tests/testthat/fixtures/nct-baseline.rds`), and under a shared seed the M, S
and per-edge p-values equal `NetworkComparisonTest::NCT()` 2.2.2 exactly.

* **Missing values no longer erase a variable.** A single `NA` made every
  correlation of that variable `NA`, which was then replaced by 0, so the
  variable lost all its edges in the observed networks and in every
  permutation, with no warning. New `na_method = c("pairwise", "listwise")`
  (default `"pairwise"`, as for the estimators): pairwise correlations with
  the median co-observation count as the sample size; listwise drops
  incomplete rows (whole pairs in a paired design).
* A constant variable now enters the network as an isolated node through an
  explicit rule (the convention of every psychnets estimator for a constant
  column) instead of a blanket `NA` -> 0 replacement under
  `suppressWarnings()`; a constant variable in an observed group raises a
  `psychnet_constant_variable` warning. A variable that becomes constant in a
  permuted split (a rare binary item) is handled the same way, so the
  statistic is defined on every split and the test stays exact.
* With pairwise missing data, a variable pair sharing fewer than two observed
  rows (in a group or a permuted split) has an undefined correlation and is
  refused with `psychnet_degenerate_sample` instead of being set to 0.
* Columns that differ in name or order raise `psychnet_bad_input` (unchanged
  rule, now a classed condition with the offending names).
* Group mode (`net_compare(psychnet(d, group = "g"))`) forwards the group's
  `gamma` and `na_method`, and warns (`psychnet_estimator_mismatch`) when the
  group was built with a non-glasso `method`, a `threshold`, or other
  estimator settings: the networks compared are then EBIC-glasso networks, not
  the group's own.
* Argument validation: `iter` must be a single whole number (2.7 was
  truncated, `"10"` was accepted); malformed arguments raise
  `psychnet_bad_input`.
* Printing the result shows everything: the two groups by name (the argument
  names, or the group levels in group mode), M and S with each group's global
  strength, the number of differing edges, and every edge (largest difference
  first) with both groups' weights, rounded p-values and `*` for p < 0.05
  (`print(fit, n = 10)` shortens the edge list). As tables, `summary()` gives
  the M and S tests and `as.data.frame()` one row per edge (`from`, `to`,
  `weight_1`, `weight_2`, `difference`, `p_value`, `significant`).
* `plot(type = "edges")` marks significant edges with `*` as well as colour.

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
