# Small-world index

The Humphries & Gurney (2008) small-world index `sigma`: observed
transitivity and average shortest-path length compared with
degree-preserving random graphs (`sigma = (C/C_rand) / (L/L_rand)`;
`sigma > 1` indicates small-worldness). Computed on the binarised
network. Ported from
[`qgraph::smallworldness`](https://rdrr.io/pkg/qgraph/man/smallworldness.html).

## Usage

``` r
net_smallworld(x, n_rand = 100L, seed = NULL)
```

## Arguments

- x:

  A `psychnet` object or a square weighted adjacency matrix.

- n_rand:

  Number of degree-preserving random graphs. Default 100.

- seed:

  Optional integer for reproducibility of the random graphs.

## Value

A tidy one-row `data.frame` (class `psychnet_smallworld`):
`smallworldness`, `transitivity`, `aspl`, `transitivity_rand`,
`aspl_rand`.

## References

Humphries, M. D., & Gurney, K. (2008). *PLoS ONE*, 3(4), e0002051.

## Examples

``` r
S <- 0.4^abs(outer(1:8, 1:8, "-"))
net_smallworld(ebic_glasso(cor_matrix = S, n = 400), n_rand = 50, seed = 1)
#>   smallworldness transitivity aspl transitivity_rand aspl_rand
#> 1              0            0    3              0.13  2.375385
```
