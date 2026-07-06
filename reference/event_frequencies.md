# Action frequencies from an event log

Converts a long event log – one row per event, with an actor, an action,
and optionally a session and time – into a frequency table: one row per
occasion, one column per distinct action, holding that action's count.
The conversion mirrors the "frequency" format of the Nestimate event
pipeline. An occasion is one (actor, session) cell; with
`compute_sessions = TRUE` the `time` column is used to split each actor
into sessions wherever the gap between consecutive events exceeds
`time_threshold`. This is the frequency input that
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
builds internally for event data (`source = "eventdata"`).

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

  A long event log (data frame), one row per event.

- actor:

  Column(s) naming the actor / subject. Default `"Actor"`.

- action:

  Column naming the action / state. Default `"Action"`.

- session:

  Optional column(s) naming an explicit session within an actor.

- time:

  Optional column used to compute sessions from gaps (see
  `compute_sessions`); it is not used for any temporal model.

- compute_sessions:

  If `TRUE`, split each actor into sessions from the `time` gaps (a new
  session starts when the gap exceeds `time_threshold`). Default `TRUE`.

- time_threshold:

  Maximum gap (in the units of `time`, seconds for a timestamp) between
  consecutive events before a new session begins. Default `900` (15
  minutes), as in Nestimate.

## Value

A `data.frame` with an `actor` column, a `session` index, and one
integer count column per action (one row per occasion).

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
