# Weighted clustering coefficients

Per-node clustering coefficients for a weighted (optionally signed)
network: Watts-Strogatz, Zhang-Horvath, Onnela, and Barrat. For signed
networks the signed variants are also returned (a closed triangle of
three negative or one negative edge lowers the signed coefficient).
Ported from
[`qgraph::clustcoef_auto`](https://rdrr.io/pkg/qgraph/man/clustcoef_auto.html).

## Usage

``` r
net_clustering(x, labels = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix.

- labels:

  Optional node labels (used when `x` is a bare matrix).

## Value

A tidy `data.frame` (class `psychnet_clustering`), one row per node:
`node`, `clustWS`, `signed_clustWS`, `clustZhang`, `signed_clustZhang`,
`clustOnnela`, `signed_clustOnnela`, `clustBarrat`.

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
