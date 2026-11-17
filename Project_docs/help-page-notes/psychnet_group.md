# Implementation notes removed from `?psychnet_group` (2026-11-17)

Still true of R/groups.R; kept for maintainers.

- Class `c("psychnet_group", "netobject_group")`, mirroring Nestimate's
  `netobject_group`; cograph::splot() grids it because it dispatches on
  `netobject_group`. The container keeps the per-level fit data
  (`attr(, "subsets")`), the settings (`attr(, "call")`), `group_col`,
  `source`, and with `within` the person ids (`persons`, `within`).
- Several grouping columns are combined into tuple levels joined by `"-"`
  (`.tuple_factor()`); with all-factor columns the levels follow lexical tuple
  order.
- Event-data groups: each level's log is reduced to its design matrix; the
  container's `source` is `"data"` only when every level reduced to independent
  rows (n equal to rows), otherwise `"eventdata"`, which the resampling verbs
  refuse (`.group_data_apply()`).
- The resampling verbs re-run on each stored subset with the original
  estimator settings (gamma, threshold for threshold methods, dots), merged
  with any `estimator_args` the caller supplies.
