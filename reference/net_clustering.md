# Compute weighted clustering coefficients of the nodes

`net_clustering()` computes weighted clustering coefficients for each
node of a network. The clustering coefficient of a node measures how
strongly its neighbours are connected to each other: in a psychometric
network, a node with high clustering lies within a closely related group
of variables, whereas a node with low clustering links variables that
are otherwise unrelated.

Several generalizations of the clustering coefficient to weighted
networks exist, and the function returns four of them: that of Watts and
Strogatz (1998), which uses only the presence of edges; of Zhang and
Horvath (2005), which uses the edge weights; of Onnela et al. (2005),
which uses the geometric mean of the weights in each triangle; and of
Barrat et al. (2004), which weights triangles by node strength. For the
first three it also returns the signed versions of Costantini and
Perugini (2014), in which a triangle with one or three negative edges
lowers the coefficient.

It returns one row per node, with one column per coefficient.

For a summary of clustering at the level of the whole network, compared
with random graphs, use
[`net_smallworld()`](https://pak.dynasite.org/psychnets/reference/net_smallworld.md);
for other node measures, use
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md).

## Usage

``` r
net_clustering(x, labels = NULL)
```

## Arguments

- x:

  An undirected fitted network
  ([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md))
  or a symmetric weighted adjacency matrix.

- labels:

  Optional node labels for a bare matrix, one unique label per node.

## Value

A data frame of class `psychnet_clustering`, one row per node, with
columns `node`, `clustWS`, `signed_clustWS`, `clustZhang`,
`signed_clustZhang`, `clustOnnela`, `signed_clustOnnela` and
`clustBarrat`.

## Details

Edge weights outside \[-1, 1\] are first divided by the largest absolute
weight. A node with fewer than two neighbours closes no triangle and has
coefficient 0. The Barrat coefficient uses absolute weights.

## Conditions

`psychnet_bad_input` when `x` is directed, is not a finite symmetric
square matrix, or `labels` is invalid.

## References

Costantini, G., & Perugini, M. (2014). Generalization of clustering
coefficients to signed correlation networks. *PLoS ONE*, 9(2), e88669.
[doi:10.1371/journal.pone.0088669](https://doi.org/10.1371/journal.pone.0088669)

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. *Nature*, 393, 440-442.
[doi:10.1038/30918](https://doi.org/10.1038/30918)

Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
co-expression network analysis. *Statistical Applications in Genetics
and Molecular Biology*, 4(1), Article 17.
[doi:10.2202/1544-6115.1128](https://doi.org/10.2202/1544-6115.1128)

Onnela, J.-P., Saramaki, J., Kertesz, J., & Kaski, K. (2005). Intensity
and coherence of motifs in weighted complex networks. *Physical Review
E*, 71(6), 065103.
[doi:10.1103/PhysRevE.71.065103](https://doi.org/10.1103/PhysRevE.71.065103)

Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A.
(2004). The architecture of complex weighted networks. *Proceedings of
the National Academy of Sciences*, 101(11), 3747-3752.
[doi:10.1073/pnas.0400087101](https://doi.org/10.1073/pnas.0400087101)

## See also

[`net_smallworld()`](https://pak.dynasite.org/psychnets/reference/net_smallworld.md)
for the small-world index;
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
for other node measures.

## Examples

``` r
S <- 0.4^abs(outer(1:6, 1:6, "-"))
net_clustering(ebic_glasso(cor_matrix = S, n = 400))
#>   node clustWS signed_clustWS clustZhang signed_clustZhang clustOnnela
#> 1   V1       0              0          0                 0           0
#> 2   V2       0              0          0                 0           0
#> 3   V3       0              0          0                 0           0
#> 4   V4       0              0          0                 0           0
#> 5   V5       0              0          0                 0           0
#> 6   V6       0              0          0                 0           0
#>   signed_clustOnnela clustBarrat
#> 1                  0           0
#> 2                  0           0
#> 3                  0           0
#> 4                  0           0
#> 5                  0           0
#> 6                  0           0
```
