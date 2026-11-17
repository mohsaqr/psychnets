# Implementation notes removed from `?tmfg_certificate` (2026-11-17)

Still true of R/tmfg.R and R/certificate.R; kept for maintainers.

- **Old framing**: "A TMFG has no convex objective, so its correctness is
  certified structurally".
- The edge-count test applies only for `p >= 3` (always true for a network
  from `tmfg_network()`, which refuses `p < 4`). Chordality uses
  maximum-cardinality search (Tarjan & Yannakakis 1984, `.is_chordal()`);
  connectivity is a breadth-first reach from node 1; planarity is
  `.tmfg_planar_structure()`.
- `certificate()` and `print()` recompute this value from the stored
  `adjacency` (TMFG objects have no stored `kkt`).
