# Implementation notes removed from `?event_frequencies` (2026-11-17)

Still true of R/multilevel.R; kept for maintainers.

- The conversion mirrors the "frequency" format of the Nestimate event
  pipeline; the session-from-gap rule is copied from `Nestimate::prepare()`
  (a new session starts when the gap to the previous event of the same
  actor/session exceeds `time_threshold`; the first event starts a session).
- `.as_time_numeric()`: numeric as is; POSIXct/Date via `as.numeric()` (so a
  Date column is in days, not seconds); character parsed with
  `as.POSIXct(tz = "UTC")`, any unparseable value is an error.
- Several `actor` / `session` columns are combined with `.tuple_factor()`.
- Rows with an NA or empty-string action are dropped silently before counting
  (not mentioned in the old help; now stated on the page).
- Counts are doubles (the old help said "numeric (double)").
