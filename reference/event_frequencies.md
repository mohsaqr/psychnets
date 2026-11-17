# Count the actions of each actor and session in an event log

`event_frequencies()` counts how often each action occurred in each
occasion of an event log, where an occasion is an actor, or an actor's
session. This table of action frequencies is what a network of actions,
which shows which actions tend to occur together in the same occasion,
is estimated from.

In an event log each row records one action by one actor, such as a
student reading, taking a quiz or writing a note; behavioural and
learning data often arrive in this form. Occasions are formed by
`actor`, by `actor` and `session` when sessions are given, and, with
`time` and `compute_sessions = TRUE`, by splitting an actor's events
into sessions wherever the gap between consecutive events exceeds
`time_threshold`.

The result is a data frame with one row per occasion and one column per
distinct action, holding the count.

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
performs this conversion itself when given an event log, so the function
is needed only to inspect or reuse the table. Events with a missing or
empty action are not counted.

## Usage

``` r
event_frequencies(
  data,
  actor = "Actor",
  action = "Action",
  session = NULL,
  time = NULL,
  compute_sessions = TRUE,
  time_threshold = 900
)
```

## Arguments

- data:

  A data frame with one row per event.

- actor:

  Column(s) identifying the actor. Default `"Actor"`.

- action:

  Column holding the action. Default `"Action"`.

- session:

  Optional column(s) identifying sessions within an actor.

- time:

  Optional column of event times (numbers, dates, date-times or
  timestamps as text), used only to split sessions.

- compute_sessions:

  Whether to split each actor's events into sessions at gaps in `time`
  longer than `time_threshold`. Default `TRUE`; has no effect without
  `time`.

- time_threshold:

  Longest gap between consecutive events within one session, in the
  units of `time` (seconds for date-times). Default `900` (15 minutes).

## Value

A data frame with one row per occasion and columns `actor`, `session`
(an occasion label joining the actor, the session and, when sessions are
computed from `time`, the session number, separated by `"."`) and one
numeric column of counts per action.

## Conditions

`psychnet_bad_input` when `data` is not a data frame, a named column is
missing, `action` has fewer than two distinct values, `actor` or
`session` has missing values, or `time` cannot be read as numbers or
timestamps.

## See also

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `actor` to estimate a network from an event log;
[psychnet_multilevel](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
for within-actor and between-actor networks.

## Examples

``` r
ev <- data.frame(
  Actor   = rep(c("a", "b"), each = 6),
  Session = rep(rep(1:2, each = 3), 2),
  Action  = c("read","quiz","read", "quiz","read","note",
              "note","note","read", "read","quiz","quiz"))
event_frequencies(ev, session = "Session")
#>   actor session note quiz read
#> 1     a     a.1    0    1    2
#> 2     a     a.2    1    1    1
#> 3     b     b.1    2    0    1
#> 4     b     b.2    0    2    1
```
