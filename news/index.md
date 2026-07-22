# Changelog

## psychnets 0.4.3

- First CRAN release.

- `mgm_fit(..., moderators = )` (the moderated mixed graphical model,
  read with
  [`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md))
  now runs on psychnets’ own base-R penalized-IRLS kernel and is
  KKT-certified. It previously required the compiled `glmnet` package
  and silently ignored `native =`. Every estimator in the package is now
  pure base R by default; `glmnet` and `glasso` remain optional
  reference engines selected with `native = FALSE`.

- The base moderated kernel covers gaussian and binary nodes (the
  documented
  [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  scope) and errors explicitly on a categorical node with more than two
  levels, which needs `native = FALSE`.

- The two resampling-heavy long-form guides, “Network reliability” and
  “Visualizing networks”, are now website articles rather than installed
  vignettes, so they keep their full iteration counts while the package
  check stays inside CRAN’s time budget. Read them at
  <https://pak.dynasite.org/psychnets>. The eight remaining vignettes
  are unchanged and still ship with the package.
