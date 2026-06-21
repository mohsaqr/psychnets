# psychnet: dependency-free, self-certifying psychometric network estimation in base R

## Abstract

`psychnet` estimates the cross-sectional network models used in psychometrics.
These are marginal and partial correlation networks, the EBIC-regularized
Gaussian graphical model and its nonparanormal and stepwise variants, two
information-filtering graphs, a relative-importance network, and the Ising and
mixed graphical models. All are reimplemented from first principles in base R.
The package's only import is the `stats` package shipped with R. It calls none
of the reference packages it reproduces, and none of the compiled numerical
kernels those packages depend on. Each regularized fit returns a correctness
certificate: the stationarity residual of the convex objective it solves, which
is the graphical-lasso objective for the Gaussian models and each nodewise
penalized regression for the Ising and mixed models. The residual is near zero
when the fit is at the optimum, and it is computed without reference to any
external solver. We describe the construction, report agreement with the
reference packages across one hundred simulated datasets, examine the cases where
`psychnet` returns a sparser and arguably more accurate result than `qgraph`, and
state the cost, which is run time. On the graphical lasso `psychnet` is about two
orders of magnitude slower than the Fortran solver used by `qgraph`. On the Ising
model it is about twenty times slower than `IsingFit`. On the mixed graphical
model it is slightly faster than `mgm`. For the problem sizes common in
psychometrics, tens of nodes, the run time is small in absolute terms. We argue
that self-certification, the absence of a compiled dependency tree, and exact
agreement with the published estimators on the methods that have a unique
solution are properties that justify that cost at those sizes.

## 1. Background and claim

A psychometric network analysis usually runs through `qgraph`, `IsingFit`,
`mgm`, or `bootnet`. These packages are correct and widely used. They are also
the visible end of a large dependency tree. Installing `bootnet` pulls more than
a hundred recursive hard dependencies, most of them compiled through `Rcpp` and
`RcppArmadillo`. The estimation itself runs in a
compiled kernel: `qgraph`'s EBIC graphical lasso calls the `glasso` Fortran
package, `IsingFit` and `mgm` call `glmnet`. The wrapper around that kernel, the
penalty path, the information criterion, and the symmetrization rule, is the part
each reference paper specifies, and it is short.

`psychnet` makes a single claim. It computes the same network models, in base R,
with no compiled dependency and no call to the reference packages, and it reports
for each regularized fit a self-contained measure of how far that fit is from the
optimum of its own objective. Where the model has a unique solution, `psychnet`
reaches it; where it does not, `psychnet` reaches the published one or states why
it cannot.

This goes further than the sibling package `Nestimate`, which
reconstructs the same wrappers but calls the same compiled kernels (`glasso`,
`glmnet`, `lme4`) and so agrees with the references to machine precision.
`psychnet` reimplements the kernels as well. That choice removes the last
compiled dependency, and it forces a different standard of proof: a from-scratch
solver cannot be trusted because it matches Fortran, since it might match a bug
in the Fortran; it has to be checked against the mathematics. Section 3 describes
how.

## 2. What clean-room means here

A grep over the shipped source settles the first part of the claim. The package
code calls no reference package and no compiled estimation kernel:

```
grep -rE "qgraph::|IsingFit::|mgm::|huge::|bootnet::|glasso::|glmnet::|Matrix::" R/
# (matches occur only inside roxygen comments, never in executable code)
```

The `DESCRIPTION` lists one import, `stats`, and one suggested package,
`testthat`. The recursive hard-dependency count is therefore the count for base
R, against seventy-one for `qgraph`, seventy-seven for `IsingFit` and `mgm`, and
more than a hundred for `bootnet` (counts as of 2026-06; the reference figures
are reproduced from the `Nestimate` technical report, which measured them on the
same CRAN snapshot, and may move by a few packages with CRAN updates). The core
estimators install and run with no C or Fortran toolchain.

The numerical work is done by two kernels written in the package. The first is a
covariance block-coordinate-descent graphical lasso (Friedman, Hastie and
Tibshirani, 2008), used by the graphical-lasso family. The second is a penalized
iteratively reweighted least squares nodewise generalized linear model (Friedman,
Hastie and Tibshirani, 2010), used by the Ising and mixed models with a logistic
or Gaussian link. Two estimators use neither kernel: the Triangulated Maximally
Filtered Graph is a greedy planar construction, and the relative-importance
network is a Shapley decomposition of nodewise R-squared.

The eleven estimator verbs, their objective, their selection rule, and their
certificate are listed below. They are distinct from the framework verbs
(`centrality`, `predictability`, `bootstrap_network`, `centrality_stability`,
`nct`), which take a fitted network and are not themselves estimators.

| verb | model | selection | certificate |
|---|---|---|---|
| `cor_network` | marginal correlation | significance threshold (optional) | none (closed form) |
| `pcor_network` | partial correlation | significance threshold (optional) | none (closed form) |
| `ebic_glasso` | graphical lasso | EBIC over a penalty path | `glasso_kkt` |
| `huge_network` | nonparanormal graphical lasso | EBIC over a penalty path | `glasso_kkt` |
| `ggm_modselect` | unregularized graph-restricted MLE | EBIC over candidate graphs | `ggm_support_kkt` |
| `tmfg_network` | maximally filtered graph | greedy planar construction | `tmfg_certificate` |
| `logo_network` | chordal Markov random field | filtered-graph support | `ggm_support_kkt` |
| `relimp_network` | relative importance (directed) | full subset enumeration | `lmg_certificate` |
| `ising_fit` | Ising, L1-penalized | per-node EBIC | `glm_lasso_kkt` |
| `ising_sampler` | Ising, unregularized | Wald pruning (optional) | `glm_lasso_kkt` |
| `mgm_fit` | mixed Gaussian and binary | per-node EBIC | `glm_lasso_kkt` |

All eleven accept the `na_method` argument; the seven Gaussian-graphical-model
verbs also accept a precomputed correlation matrix in place of raw data.

## 3. Self-certification

The graphical lasso minimizes a strictly convex function of the precision matrix,

    min over Theta positive definite of  -log det Theta + tr(S Theta)
                                          + rho * sum over i != j of |Theta_ij|.

Strict convexity means the minimizer is unique, and it means the subgradient
optimality conditions characterize it exactly. Writing W for the inverse of
Theta, those conditions are: the diagonal of W equals the diagonal of S; on an
edge that is present, W_ij minus S_ij equals rho times the sign of Theta_ij; on
an edge that is absent, the absolute value of W_ij minus S_ij is at most rho. The
exported function `glasso_kkt(theta, S, rho)` evaluates the largest violation of
these conditions and returns it. A return near zero certifies that the supplied
precision matrix is the global optimum, and it certifies this against the
objective itself, with no reference solver in the computation.

Every `ebic_glasso` result carries this number in its `$kkt` field. On the fits
in this paper it is near zero, below 1e-9, and reaches exact zero when the
selected graph is empty. The nodewise generalized
linear model has an analogous certificate, `glm_lasso_kkt`, built from the
penalized-likelihood score conditions; the Ising and mixed models store the worst
nodewise value. The unregularized estimators that have a closed form carry a
structural certificate instead: the constrained Gaussian Markov random field
behind `ggm_modselect` and `logo_network` is checked by `ggm_support_kkt`, which
verifies that W matches S on the retained edges and that the precision is exactly
zero off them; the filtered graph is checked for the planar edge count and
chordality; the relative-importance shares are checked against the efficiency
identity that they sum to each node's full-model R-squared.

The certificate is the device that makes a base-R reimplementation auditable. It
replaces the argument "trust this because it matches the reference" with the
argument "trust this because its distance from the unique optimum is printed and
is near zero."

## 4. Where psychnet may return the better answer

The reference packages are well established, and on most data `psychnet` agrees
with them closely. There are specific, reproducible situations where `psychnet` is the
more accurate of the two. They follow from the same source: `psychnet` solves to
the optimum and uses the textbook information criterion, while a reference solver
stops at a loose tolerance.

### 4.1 Optimality

`qgraph` runs `glasso` at its default convergence threshold of 1e-4, so its
returned precision leaves a residual consistent with that threshold. Applying
`glasso_kkt` to `qgraph`'s own returned precision confirms this. `psychnet`
refits the selected penalty to a stationarity residual near 1e-11. On a fixed
penalty the two agree to roughly 1e-6, and that residual is `qgraph`'s, not
`psychnet`'s, since the `psychnet` fit has the smaller certificate and is by that
measure closer to the optimum.

### 4.2 The textbook EBIC and a false positive

The extended Bayesian information criterion of the empty Gaussian graphical model
on a correlation matrix has a closed form. The empty model has precision equal to
the identity, so its log-likelihood is minus n times p over two, no parameters
are penalized, and the criterion equals n times p exactly. For a six-node problem
at n = 120 that is 720. `psychnet` returns 720.00. In our benchmark `qgraph`'s
default returns 728.37, because its default criterion is computed from a
model-fit object rather than from the closed-form Gaussian likelihood, and
because its glasso leaves a small residual at the boundary where the path
empties. The constant offset between the two criteria does not change which
penalty is selected, but the solver residual does.

We examined the one dataset in a hundred-dataset comparison where this changed
the result by a visible amount. It is a sparse six-node graph at n = 120. The
true conditional graph has six edges, all with partial correlations below 0.18,
which is below the detection threshold at that sample size; neither estimator
recovers them. The two disagree on a single edge, V2 to V5. The true partial
correlation of that edge is zero. The sample correlation is minus 0.24, which is
noise. `psychnet` returns the empty graph. `qgraph` returns the edge at minus
0.041. On this dataset `psychnet` is correct and `qgraph` reports a false
positive. The gap that keeps EBICglasso from being labeled "exact" against
`qgraph` in our comparison is, on inspection, `qgraph` keeping a spurious edge
that `psychnet` declined to keep.

### 4.3 The canonical filtered graph

The Triangulated Maximally Filtered Graph is a deterministic greedy
construction. `psychnet` reproduces the published algorithm (Massara, Di Matteo
and Aste, 2016) exactly: a faithful re-derivation of the textbook procedure
agrees with `psychnet` on every edge across the test datasets. The same datasets
show `NetworkToolbox`'s implementation differing from the textbook construction by
about a fifth of the edges, which traces to a different tie-handling rule in its
incremental gain table. Here `psychnet` matches the specification and the
reference package is the one that departs from it.

### 4.4 Missing data

The reference correlation routines retain incomplete cases through pairwise or
full-information correlations. `psychnet` does the same and makes pairwise the
default. The estimators that drop incomplete rows fail abruptly when missingness
is spread across many columns, because listwise deletion removes any row with one
missing cell, and the surviving sample falls below the number of nodes. In a
simulation at n = 150, p = 10, and fifteen percent missing completely at random,
listwise deletion left fewer than thirty complete rows and the graphical lasso
returned an empty graph, an F-measure of zero against the true chain. The pairwise
correlation retained the information in all rows and recovered the structure at an
F-measure near 0.86. With complete data the two settings return identical results,
so the default carries no cost on data that has no missing values.

### 4.5 Dependency footprint, reproducibility, portability

The absence of a compiled dependency tree is itself a property of the answer.
`psychnet` installs without a toolchain, runs where base R runs, and produces a
result that does not depend on the version of a Fortran library linked at build
time. The estimators that have a unique solution are deterministic and carry a
printed certificate, so a result can be rechecked from the object alone.

## 5. The cost: run time

Reimplementing the kernel in interpreted R is slower than calling Fortran. The
table reports median elapsed seconds over repeated fits on an Apple-silicon
laptop, from the benchmark script in the package's local equivalence directory.
The reference is the package each estimator reproduces. The ratios depend on the
machine and the reference solver's own settings; the absolute figures and the
orders of magnitude are the stable part.

| estimator | size (p, n) | psychnet (s) | reference (s) | ratio |
|---|---|---:|---:|---:|
| EBICglasso | 10, 250 | 0.29 | 0.002 | 146 |
| EBICglasso | 20, 500 | 1.33 | 0.007 | 190 |
| EBICglasso | 30, 500 | 3.22 | 0.017 | 190 |
| Ising | 6, 500 | 0.26 | 0.013 | 20 |
| mgm | 4, 500 | 0.14 | 0.23 | 0.6 |

The graphical lasso carries the largest penalty, about two orders of magnitude,
because the pure-R block-coordinate descent competes against a Fortran kernel.
The Ising model is about twenty times slower than `IsingFit`, which calls
`glmnet`. The mixed model is faster than `mgm`, because `mgm`'s own wrapper
overhead exceeds the cost of the base-R nodewise loop at this size.

The absolute figures matter more than the ratios. A thirty-node graphical lasso
fits in about three seconds. Psychometric networks are usually smaller than
fifty nodes, where the fit is a few seconds at most, and the bottleneck in
applied work is the bootstrap, which `psychnet` runs by resampling and refitting
in the same base-R path. Where the problem is large, hundreds of nodes, the
Fortran kernel is the appropriate tool and `psychnet` is not. The package trades
run time for the absence of a dependency tree, a printed proof of optimality, and
the correctness properties of Section 4. For typical psychometric problem sizes
that trade is reasonable; for large problems it is not, as stated here.

## 6. Validation

Agreement with the reference packages was measured on one hundred datasets
generated with the `Saqrlab` simulation tools, of which ninety-nine retained
enough complete cases to estimate and were compared. The datasets vary the
covariance structure (autoregressive, compound, two-cluster, factor, sparse),
the sample size, the node type (continuous, binary, mixed), and the missingness
mechanism (none, MCAR, MAR, MNAR). Each was estimated with `psychnet` and with
the matching reference package on the same complete-case input, and we recorded
the largest absolute edge difference and the agreement of the zero/nonzero
pattern. The per-estimator row counts in the table below are the datasets to
which each estimator applies.

| estimator | reference | mean delta | max delta | structure agreement |
|---|---|---:|---:|---:|
| pcor | base solve | 0 | 0 | 1.000 |
| relimp | relaimpo (LMG) | 0 | 0 | 1.000 |
| TMFG | textbook construction | 0 | 0 | 1.000 |
| LoGo | NetworkToolbox | 0 | 0 | 1.000 |
| EBICglasso | qgraph | 0.0008 | 0.041 | 0.997 |
| Ising | IsingFit | 0.029 | 0.226 | 0.995 |
| huge | huge | 0.034 | 0.107 | 0.964 |
| mgm | mgm | 0.075 | 0.204 | 0.958 |

The four estimators that have a unique closed-form or deterministic solution
agree exactly. The regularized estimators agree on the edge set in 96 to 99.7
percent of comparisons; the residual differences are penalty-selection flips on
small or near-empty datasets, the worst of which was the false positive examined
in Section 4.2. For the graphical lasso, sixty-one of the sixty-three datasets
agree to solver precision or exactly, and the two visible differences are both
sparse six-node graphs at n = 120.

Two convention differences were found and resolved during this work. The Gaussian
nodewise EBIC in the mixed model used the profiled-variance deviance n times log
of the residual sum of squares over n, while `glmnet` and `mgm` use the residual
sum of squares itself; the two differ by a logarithm that compresses the penalty
and selects denser graphs. Aligning the deviance to the residual sum of squares
brought the mixed model into edge-for-edge agreement with `mgm` on continuous
data and preserved exact ground-truth recovery on a known chain. The mixed model
also gained the Loh-Wainwright post-selection threshold that `mgm` applies by
default. Both changes are documented in the package change log.

One estimator was removed rather than reconciled. A non-convex graphical lasso
(SCAD, MCP, atan penalty) has no unique solution path, so its one-step local
linear approximation depends on the warm start, the penalty grid, and the
derivative parameterization; it differed from the `GGMncv` package by about 0.2
even on identical input, while recovering structure at least as well. Because it
could not be made reproducible across implementations by construction, and
because the reference package is little used, it was dropped. Convexity is what
makes cross-package agreement possible; the one estimator that lacked it was the
one estimator that could not agree.

## 7. What is not claimed

`psychnet` is not a new statistical method and claims no improved estimator
behavior beyond the specific reproducible cases in Section 4, which follow from
solving to the optimum and using the textbook criterion rather than from any new
model. It computes the established quantities, in base R, with a printed
certificate and without a compiled dependency, and it does so more slowly than
the compiled references. The graphical lasso cannot be made byte-identical to
`qgraph` in dependency-free code, because the residual difference is `qgraph`'s
Fortran solver tolerance and not a formula; reproducing it would mean depending
on `glasso` or reproducing its imprecision, and the package does neither. The
temporal models (`graphicalVAR`, `mlVAR`) are out of scope. The mixed model
supports Gaussian and binary nodes only, not categorical variables with more than
two levels or count variables. The network comparison test is defined for
Gaussian graphical models only. The package provides no plotting methods; a fitted
network is returned as a tidy edge list for plotting elsewhere.

## 8. Reproducibility

The package ships 180 test expectations across fifty-four test cases and passes
`R CMD check` with no errors, warnings, or notes. The equivalence comparison and
the timing benchmark are scripts in a
build-ignored, local-only directory, run against the installed reference
packages (`qgraph` 1.9.8, `huge` 1.3.5, `NetworkToolbox` 1.4.4, `relaimpo`
2.2.7, `IsingFit` 0.4, `mgm` 1.2.15). The certificate functions are exported, so
any fit can be rechecked: `glasso_kkt`, `glm_lasso_kkt`, `ggm_support_kkt`,
`tmfg_certificate`, and `lmg_certificate`.

## References

Chen, J., and Chen, Z. (2008). Extended Bayesian information criteria for model
selection. Biometrika, 95(3), 759-771.

Foygel, R., and Drton, M. (2010). Extended Bayesian information criteria for
Gaussian graphical models. Advances in Neural Information Processing Systems, 23.

Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse covariance
estimation with the graphical lasso. Biostatistics, 9(3), 432-441.

Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths for
generalized linear models via coordinate descent. Journal of Statistical
Software, 33(1), 1-22.

Haslbeck, J. M. B., and Waldorp, L. J. (2020). mgm: Estimating time-varying
mixed graphical models in high-dimensional data. Journal of Statistical
Software, 93(8), 1-46.

Massara, G. P., Di Matteo, T., and Aste, T. (2016). Network filtering for big
data: Triangulated Maximally Filtered Graph. Journal of Complex Networks, 5(2),
161-178.

Barfuss, W., Massara, G. P., Di Matteo, T., and Aste, T. (2016). Parsimonious
modeling with information filtering networks. Physical Review E, 94(6), 062306.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L.,
Schoevers, R. A., and Waldorp, L. J. (2014). A new method for constructing
networks from binary data. Scientific Reports, 4, 5918.

Grömping, U. (2006). Relative importance for linear regression in R: the package
relaimpo. Journal of Statistical Software, 17(1), 1-27.

Liu, H., Lafferty, J., and Wasserman, L. (2009). The nonparanormal: semiparametric
estimation of high-dimensional undirected graphs. Journal of Machine Learning
Research, 10, 2295-2328.

Epskamp, S., Borsboom, D., and Fried, E. I. (2018). Estimating psychological
networks and their accuracy: a tutorial paper. Behavior Research Methods, 50(1),
195-212.
