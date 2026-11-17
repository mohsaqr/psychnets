# Check that a network has the structure of a TMFG

`tmfg_certificate()` checks that a network from
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md)
has the structure of a triangulated maximally filtered graph (TMFG). A
TMFG is defined by its structure rather than as the optimum of an
objective: for \\p\\ variables it is a connected planar graph with
exactly \\3(p - 2)\\ edges that is chordal (every cycle of four or more
nodes has a chord) and is built from four-variable cliques (Massara et
al., 2016).

These properties are what make the graph a TMFG and what the LoGo model
([`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md))
relies on, so checking them shows that a returned network is a valid
TMFG, for instance that no edge was lost or added. `tmfg_certificate()`
counts how many of the defining properties the network violates.

It returns that count: 0 confirms a valid TMFG.

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the same value for a fitted TMFG network. The check does not
assess which correlations the construction chose to retain.

## Usage

``` r
tmfg_certificate(x)
```

## Arguments

- x:

  A TMFG network from
  [`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md)
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

## Value

A single integer from 0 to 4, the number of violated properties; 0
confirms a valid TMFG.

## Details

Four properties are checked on the stored graph: the number of edges is
\\3(p - 2)\\; the graph is chordal, tested by maximum-cardinality
search; it is connected; and it is planar. For a connected chordal graph
with \\3(p - 2)\\ edges, planarity is checked as the structure of a
stacked triangulation: no triangle has more than two common neighbours,
and no triangle has two common neighbours that are adjacent to each
other.

## Conditions

`psychnet_bad_input` when `x` is not a network from
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md).

## References

Massara, G. P., Di Matteo, T., & Aste, T. (2016). Network filtering for
big data: Triangulated maximally filtered graph. *Journal of Complex
Networks*, 5(2), 161-178.
[doi:10.1093/comnet/cnw015](https://doi.org/10.1093/comnet/cnw015)

## See also

[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md);
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md);
[`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md),
the other structural check.

## Examples

``` r
set.seed(1)
x <- matrix(stats::rnorm(200 * 6), 200, 6)
tmfg_certificate(tmfg_network(x))
#> [1] 0
```
