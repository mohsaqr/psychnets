# Implementation notes removed from `?logo_network` (2026-11-17)

Still true of R/logo.R and R/tmfg.R; kept for maintainers.

- **Old framing**: "Equivalent in purpose to `NetworkToolbox::LoGo()` /
  `bootnet`'s `"LoGo"` default, pure base R and self-certified via
  ggm_support_kkt() (the precision reproduces `S` exactly on the TMFG
  support)."
- **TMFG input**: the graph is built by `.tmfg_build(abs(S))`; the precision
  is `.logo_precision()`. A pseudo-inverse is deliberately not used for a
  singular block (it would not give the decomposable MLE); a block with
  `rcond < 1e-10` is refused because a near-duplicate column once gave a
  certificate of 0.14 without that check.
- **Stored fields** not named on the new page: `fit_kkt`,
  `certificate_target`.
