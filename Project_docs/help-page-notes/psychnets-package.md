# Implementation notes removed from `?psychnets-package` (2026-11-17)

Still true of the package; kept for maintainers.

- Old title: "psychnets: Clean-Room Base-R Psychometric Network Estimation".
  Old text: "Estimates psychometric network models in pure base R, with no
  compiled dependencies and dependency-free correctness certificates."
  (Imports are base packages only; `glasso` and `glmnet` are Suggests used by
  the `native = FALSE` engines.)
- Old "Certification" section: "Regularized estimators are graded against
  their own objective rather than an external solver. For the Gaussian
  graphical model, glasso_kkt() gives the stationarity (KKT) residual: zero
  certifies the unique global optimum of the convex objective." Now a
  "Certificates" section in plain words.
