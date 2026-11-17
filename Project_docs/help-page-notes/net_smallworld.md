# Implementation notes removed from `?net_smallworld` (2026-11-17)

Still true of R/graph_metrics.R; kept for maintainers.

- Ported from `qgraph::smallworldness()`: the observed transitivity and path
  length match exactly; the random reference is stochastic.
- Random graphs (`.psn_rewire()`): double-edge swaps with a fixed number of
  proposals, `max(10, 10 * m)` for `m` edges; rejected moves are retained
  (stopping after a fixed number of accepted swaps would sample the jump chain
  and overweight graphs with more available swaps). One endpoint orientation
  is randomized per proposal so both reconnections can occur. Swaps that
  share an endpoint or would duplicate an edge are rejected.
- Path length (`.psn_aspl()`): Floyd-Warshall on unit weights, mean over
  finite upper-triangle distances; `NA` when there is none.
- Reference means use `na.rm = TRUE`; sigma is `NA` (with
  `psychnet_undefined_statistic`) when the reference transitivity or path
  length is non-finite or below `sqrt(.Machine$double.eps)`.
- `seed`: the caller's `.Random.seed` is saved and restored on exit (removed
  again if it did not exist). `seed` itself is not validated.
