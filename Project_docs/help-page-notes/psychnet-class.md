# Implementation notes removed from `?psychnet-class` (2026-11-17)

Still true of R/result.R; kept for maintainers.

- Class `c("psychnet", "cograph_network")`, built by `.new_psychnet()`. The
  `cograph_network` class is what lets `cograph::splot()` dispatch on it; the
  node table also has `x`/`y` layout slots (NA), and the list carries `data =
  NULL`, `meta = list(source = "psychnet")` and `node_groups = NULL`.
- The stored `edges` use integer node indices (the cograph format) and drop
  weights with `abs(w) <= 1e-12`; `as.data.frame()` rebuilds a label-based
  edge list from `weights` with the same tolerance (upper triangle for
  undirected networks, all off-diagonal cells for directed ones).
- `$.psychnet` aliases: `graph` -> `weights`, `labels` -> `nodes$label`,
  `n_nodes` -> `nrow(nodes)`, `n_edges` -> `nrow(edges)`, `n_obs` -> `n`.
  Stored fields take precedence.
- `print()` labels: the structural certificate is printed as "violated TMFG
  invariants" (an integer count) for TMFG and "decomposition residual" for
  relimp; the others print "optimality (KKT residual)", or "unavailable for
  post-threshold graph (pre-threshold residual ...)" when the
  `certificate_target` is `"pre_threshold_fit"`.
- Other method-specific fields: `support`, `nodewise`, `fit_kkt`,
  `certificate_target`, `ebic`, `native`, `p_values`, `alpha`, `rule`.
- Node predictability is attached at fit time (`.attach_predictability()`) so
  `cograph::splot()` can draw the predictability ring.
