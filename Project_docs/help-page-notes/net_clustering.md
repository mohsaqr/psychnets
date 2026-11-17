# Implementation notes removed from `?net_clustering` (2026-11-17)

Still true of R/graph_metrics.R; kept for maintainers.

- Ported from `qgraph::clustcoef_auto()`; matches it to about 1e-8.
- Rescaling: weights are divided by `max(abs(W))` only when some weight lies
  outside \[-1, 1\].
- Formulas use `.psn_diag_cube(A)` = `diag(A %*% A %*% A)`: Watts-Strogatz on
  `sign(W)` with degree `k(k-1)`; Zhang-Horvath on `W` with
  `colSums(|W|)^2 - colSums(W^2)`; Onnela on sign-preserving cube roots
  `|w|^(1/3)` with `k(k-1)`; Barrat on `|W|` with `s (k - 1)`.
- 0/0 (fewer than two neighbours) is reported as 0, as qgraph does.
- Symmetry tolerance `1e-8`.
