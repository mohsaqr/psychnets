# Implementation notes removed from `?net_bridge` (2026-11-17)

Still true of R/bridge.R; kept for maintainers.

- Ported from `networktools::bridge()` for the undirected weighted case;
  matches it to about 1e-8.
- Symmetry is checked with tolerance `1e-8`; the diagonal is zeroed.
- Shortest paths run on the signed matrix through `.psn_floyd_warshall()`,
  which keeps only positive weights (networktools' `g2`), distances `1 / w`.
- Bridge betweenness counts geodesics (`sigma(s,i) * sigma(i,t)`), not
  fractions, over unordered pairs of nodes in different communities (halved
  from ordered pairs); on-path tolerance `1e-10`.
- Bridge closeness is 0 for a node with no reachable other-community node.
- Communities are resolved by `.psn_resolve_communities()` (shared with
  `net_aggregate()`); unnamed list elements are named by position.
