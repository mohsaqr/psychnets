# Implementation notes removed from `?net_edge_betweenness` (2026-11-17)

Still true of R/centrality.R; kept for maintainers.

- Matches `igraph::edge_betweenness()` on unweighted graphs to about 1e-16.
- An edge a->b lies on a shortest s->t path when
  `|d(s,a) + len(a,b) + d(b,t) - d(s,t)| < 1e-9`; it carries
  `sigma(s,a) * sigma(b,t) / sigma(s,t)` of the geodesics.
- Name shared with Nestimate (same contract, `x, invert, labels`); Nestimate
  delegates to this function, so it was not renamed in the 2026-08-07 name
  reconciliation.
