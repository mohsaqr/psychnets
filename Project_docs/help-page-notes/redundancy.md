# Implementation notes removed from `?redundancy` (2026-11-17)

Still true of R/redundancy.R and R/aggregate.R; kept for maintainers.

- Ported from `networktools::goldbricker()`; the proportion matrix matches it
  exactly.
- The test of two dependent, overlapping correlations is that of Hittner, May
  and Silver (2003), reimplemented in base R (`.psn_hittner2003()`) and
  validated against `cocor`. Comparisons with |r| = 1 (an item with itself)
  give `NA` and are excluded from the proportion; their NaN warnings are
  suppressed inside the helper, as cocor/goldbricker do. (The page names the
  test only in words; Hittner et al. is not in the reference list.)
- Flagging uses strict inequalities: `proportion < threshold` and
  `correlation > cor_min`. Rows are ordered by increasing proportion, ties in
  upper-triangle order.
- Data preparation is `.psn_item_data()` (shared with `net_aggregate()`):
  labels, `.as_numeric_matrix()` (ordered factors become integer codes,
  other non-numeric columns are dropped with a `psychnet_dropped_columns`
  message), listwise deletion with `psychnet_rows_dropped`.
- Attributes stored: `proportion_matrix`, `p`, `threshold`, `cor_min`.
