# psychnet validation

This directory holds a reproducible comparison of `psychnet`'s network
estimators against the established reference packages, on real questionnaire
data and on synthetic data with a known generating graph. It is not part of the
package build (the directory is build-ignored); it is the evidence behind the
equivalence claims.

## What it does

`run_validation.R` has three parts.

**Part A. EBIC graphical lasso, psychnet vs qgraph, on real questionnaire data.**
Each dataset is a set of item responses from a published instrument (Big Five,
depression, PTSD, state and trait anxiety, NEO openness, intelligence batteries,
and others). Both estimators receive the same pairwise correlation matrix and
effective sample size, so the comparison isolates the penalty path and the EBIC
selection from the upstream choice of correlation. We record the largest absolute
edge-weight difference and the agreement of the zero/nonzero pattern.

**Part B. Ising model, psychnet vs IsingFit, on real binary data.** For the
binary ability and intelligence items, both estimators fit nodewise from the same
binarized responses.

**Part C. Synthetic ground truth.** Real data cannot say which estimator is
closer to the truth, because the true network is unknown. Part C generates data
from a known sparse precision matrix and reports both the agreement with qgraph
and the F-measure of each estimator against the true edge set.

## How to run

From the package root:

```
Rscript validation/run_validation.R
```

It loads `psychnet` with `devtools::load_all()` if available, otherwise the
installed package. The comparison needs these packages installed (they are not
dependencies of `psychnet` itself): `qgraph`, `psych`, `psychTools`, `mgm`,
`NetworkToolbox`, `networktools`, `EGAnet`, `IsingFit`.

## Output

- `results_realdata.csv` — Part A, one row per dataset, with the source citation.
- `results_ising.csv` — Part B.
- `results_synthetic.csv` — Part C, with the recovery F-measures.
- `RESULTS.md` — a short text summary.

## Datasets and sources

The instruments and their sources are listed in the `citation` column of
`results_realdata.csv` and `results_ising.csv`. Each dataset ships with the R
package named in its label; see that package's documentation for the full
reference and licence. No data is downloaded; everything is loaded from installed
CRAN packages, which is what makes the comparison reproducible from a fixed set
of package versions. The exact versions used are printed at the end of the run.

## Interpreting the result

On real questionnaire data the two estimators select the same edge set and agree
on the weights to the reference solver's own convergence tolerance. The synthetic
part shows that where the truth is known, they recover it equally well. The
documented exceptions, both on synthetic near-empty graphs at small sample size,
are described in the package paper (`paper/psychnet-clean-room.md`).
