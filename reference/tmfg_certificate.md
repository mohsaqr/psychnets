# Structural certificate for a TMFG network

A TMFG has no convex objective, so its correctness is certified
structurally: a valid TMFG on `p >= 3` nodes has exactly `3(p - 2)`
edges, is connected, and is chordal (every cycle of length \>= 4 has a
chord). Returns a non-negative score that is `0` for a valid TMFG.

## Usage

``` r
tmfg_certificate(x)
```

## Arguments

- x:

  A [psychnet](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  object produced by
  [`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md).

## Value

Scalar; `0` certifies a valid TMFG (correct edge count, connected,
chordal), otherwise a positive integer counting the violated invariants.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 6), 200, 6)
tmfg_certificate(tmfg_network(x))
#> [1] 0
```
