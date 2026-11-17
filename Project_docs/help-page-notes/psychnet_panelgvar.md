# Implementation notes removed from `?psychnet_panelgvar` (2026-11-17)

Still true of R/panelgvar.R; kept for maintainers.

- **Storage**: the object is a list of three `psychnet` networks (method
  `"panel_gvar"`); everything else lives in attributes (`parameters`, `fit`,
  `kkt`, `n`, `starts`, `boundary`, `variance`, `waves`, `settings`,
  `group_col = "network"`), so `cograph::splot()` can grid the list. The
  temporal network stores `t(beta)` as weights (edge col -> row) and the
  full `beta` (with the autoregressive diagonal); the other two store their
  covariance (`sigma`) on the data scale.
- **Certificate access**: the individual networks carry no `kkt`, so
  `certificate()` on one of them returns `NA`; the certificate is the
  `score` column of `as.data.frame(x, what = "diagnostics")` (the `kkt`
  attribute). The old help described it as "the ML `score` certificate".
- **Diagnostics**: `starts_agree` counts starts within 1e-6 (relative) of
  the best standardized objective; `se_available` is `TRUE` when no free
  parameter has an `NA` standard error.
- **Variance split**: computed at fit time on the standardized scale
  (`Sigma_between` diagonal and stationary `Sigma0` diagonal) and rescaled
  by s^2.
- **Path `selected`**: the row whose lambda equals the stored lambda within
  `sqrt(.Machine$double.eps)` relative.
- **Plot**: circle layout with angle `pi/2 - 2*pi*(i-1)/p`; edge colours
  `#0072B2` (positive) and `#D55E00` (negative), node fill `#56B4E9`; edge
  style 1/2 for significant/not; plot-argument abbreviations are resolved
  against `cograph::splot()` formals before merging with the defaults.
