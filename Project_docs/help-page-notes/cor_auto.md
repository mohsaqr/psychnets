# Implementation notes removed from `?cor_auto` (2026-11-17)

Still true of R/cor_auto.R; kept for maintainers.

- **Old framing**: "The base-R counterpart of `qgraph::cor_auto()`; this is
  the correlation `bootnet`/`qgraph` use by default for Likert data."
- **Parity claim** (old help): the polychoric "matches `psych`/`qgraph` to
  about 1e-5". The polyserial difference (up to ~0.006) is kept on the page
  because it changes results.
- **Bivariate normal CDF**: Gauss-Legendre quadrature (60 nodes, Golub-Welsch)
  along the correlation path (Drezner & Wesolowsky 1990); `optimize()` over
  rho in (-0.999, 0.999), `tol = 1e-6`. Cell probabilities are floored at
  `1e-12` in the log-likelihood. Polyserial is clamped to [-0.999, 0.999];
  non-finite pair values become 0. The matrix is projected with
  `.nearest_pd_cor()`.
