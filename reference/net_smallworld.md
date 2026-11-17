# Compute the small-world index of a network

`net_smallworld()` computes the small-world index sigma of a network. A
small-world network combines the high clustering of a regular lattice
with the short paths between nodes of a random graph, and Humphries and
Gurney (2008) quantify this with sigma: the transitivity of the network
relative to that of comparable random graphs, divided by its average
shortest-path length relative to theirs,
`sigma = (C / C_rand) / (L / L_rand)`. A value of sigma above 1
indicates small-world structure.

The function computes the index on the unweighted network, in which an
edge is present when its weight is nonzero, and compares it with
`n_rand` random graphs that keep the degree of every node, obtained by
rewiring pairs of edges. Because the reference graphs are random, the
result varies between calls unless `seed` is given.

It returns one row with sigma, the transitivity and average
shortest-path length of the network, and their means over the random
graphs.

For the clustering of each node, use
[`net_clustering()`](https://pak.dynasite.org/psychnets/reference/net_clustering.md).

## Usage

``` r
net_smallworld(x, n_rand = 100L, seed = NULL)
```

## Arguments

- x:

  An undirected fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
  or a symmetric weighted adjacency matrix.

- n_rand:

  Number of random reference graphs. Default `100`.

- seed:

  Optional seed for the random graphs. The caller's random number stream
  is restored afterwards. Default `NULL` (the current stream).

## Value

A one-row data frame of class `psychnet_smallworld`, with columns
`smallworldness` (sigma; `NA` when undefined), `transitivity`, `aspl`
(average shortest-path length), and `transitivity_rand` and `aspl_rand`
(their means over the random graphs).

## Details

Transitivity is the global clustering coefficient: three times the
number of triangles divided by the number of connected triples. The
average shortest-path length is taken over the pairs of nodes joined by
a path. The index is undefined when the random graphs have no triangles
or no paths.

## Conditions

`psychnet_bad_input` when `x` is directed, is not a finite symmetric
square matrix of at least two nodes, or `n_rand` is not a positive whole
number; `psychnet_undefined_statistic` (warning) when the random graphs
have no triangles or no paths, so that sigma is `NA`.

## References

Humphries, M. D., & Gurney, K. (2008). Network 'small-world-ness': A
quantitative method for determining canonical network equivalence. *PLoS
ONE*, 3(4), e0002051.
[doi:10.1371/journal.pone.0002051](https://doi.org/10.1371/journal.pone.0002051)

## See also

[`net_clustering()`](https://pak.dynasite.org/psychnets/reference/net_clustering.md)
for clustering per node.

## Examples

``` r
S <- 0.4^abs(outer(1:8, 1:8, "-"))
net_smallworld(ebic_glasso(cor_matrix = S, n = 400), n_rand = 50, seed = 1)
#>   smallworldness transitivity aspl transitivity_rand aspl_rand
#> 1              0            0    3              0.04  2.735769
```
