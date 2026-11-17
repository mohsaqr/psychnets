# Implementation notes removed from `?net_centralities` (2026-11-17)

Still true of R/centrality.R and R/plot.R; kept for maintainers.

- **Shortest paths**: `.psn_floyd_warshall()` (igraph-free, from Nestimate),
  distances `1 / |w|` over edges with `|w| > 0`, path counts accumulated for
  ties; a node lies on a geodesic when `|d(s,v) + d(v,t) - d(s,t)| < 1e-10`.
- **Betweenness** (`.psn_betweenness()`): fraction of geodesics through the
  node, summed over ordered pairs and halved, divided by `(p-1)(p-2)/2`; 0 for
  fewer than three nodes. The halving assumes an undirected graph; on a
  directed network (relimp) the shortest paths are directed but the halving is
  still applied (observation, not documented on the page).
- **Closeness** (`.psn_closeness()`): reachable count / sum of distances;
  isolated nodes score 0.
- **centrality_fn alignment** (`.psn_align_custom()`): numeric or logical
  values; a named vector must carry exactly the node labels (any order,
  no duplicates).
- **Plot**: `type = "line"` orders nodes by mean z-score over the measures (a
  missing value is a gap), ties by node name; panels keep their own axes.
