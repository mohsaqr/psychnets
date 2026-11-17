# Estimate an unregularized Ising network for binary items

`ising_sampler()` estimates an unregularized Ising network from binary
items, in which edges can be kept or removed by their statistical
significance rather than by a penalty. The Ising model describes binary
items through a threshold for each item and an interaction for each
pair, the association between two items after conditioning on all others
(van Borkulo et al., 2014); estimated without a penalty, the
interactions are maximum-likelihood estimates, not shrunk toward zero.

The function regresses each item on all other items by ordinary
(unpenalized) logistic regression. With `alpha`, every coefficient whose
Wald p-value, optionally adjusted for multiplicity, exceeds `alpha` is
set to zero. Every pair is estimated twice, once in each item's
regression; the AND rule (default) keeps an edge when both estimates
remain nonzero, the OR rule when either does, and the edge weight is the
mean of the two. The certificate is the largest violation of the
maximum-likelihood optimality conditions (the score equations) across
the regressions, so a value near zero shows that every regression is at
its maximum. Pruning sets coefficients to zero after fitting, without
refitting; the pruned network is no longer the solution that was
checked, and its certificate is reported as missing.

The result is an undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).

Unpenalized estimates need many more observations than items. When an
item is perfectly predicted by the others (separation), its
maximum-likelihood estimate does not exist; the penalized
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
is the alternative then, and in small samples generally.

## Usage

``` r
ising_sampler(
  data,
  rule = c("AND", "OR"),
  alpha = NULL,
  adjust = "none",
  min_sum = NULL,
  weights = NULL,
  na_method = c("pairwise", "listwise"),
  labels = NULL
)
```

## Arguments

- data:

  Binary data: a data frame or matrix with one row per observation, with
  numeric 0/1, logical, or two-level factor columns, as in
  [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md).
  A constant column is dropped.

- rule:

  `"AND"` (default) keeps an edge when both nodewise estimates are
  nonzero; `"OR"` when either is.

- alpha:

  Significance level for pruning, a number in (0, 1). Default `NULL`
  keeps every edge.

- adjust:

  Multiplicity correction for the p-values used in pruning, one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html).
  Default `"none"`. Used only with `alpha`.

- min_sum:

  Minimum number of endorsed items per row; rows below it are dropped
  before fitting. Default `NULL` keeps every row.

- weights:

  Optional non-negative observation weights, one per row after
  missing-data handling (and filtered with `min_sum`). Only their
  relative sizes matter. Default `NULL` (unweighted).

- na_method:

  `"pairwise"` (default) fills each missing value with the item's more
  frequent value; `"listwise"` drops rows with a missing value.

- labels:

  Optional node labels: one per node, or one per input column.

## Value

An undirected network of class `psychnet`
([psychnet](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per edge (`from`, `to`, `weight`),
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
the maximum-likelihood certificate, and
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
the predictability of each node. The certificate is missing when pruning
removed a coefficient or when a regression did not converge. The object
also stores the node thresholds (`thresholds`), the unadjusted Wald
p-values of the nodewise coefficients (`p_values`, row = the regressed
item), the certificate of the unpruned fit (`fit_kkt`) and what the
certificate refers to (`certificate_target`: `"returned_network"`,
`"pre_pruning_fit"` or `"nonconverged_fit"`).

## Details

Edge weights and thresholds are on the log-odds scale of the 0/1 coding,
as in
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md).
The multiplicity correction (`adjust`) is applied over all nodewise
coefficients, two per pair, before the AND or OR rule. After pruning,
the thresholds are recomputed from the retained coefficients; the
remaining coefficients keep their unpruned estimates.

## Conditions

`psychnet_bad_input` when a column is not binary, fewer than two usable
columns remain, `rule`, `alpha`, `adjust`, `min_sum`, `weights` or
`labels` is invalid, or `min_sum` removes every row (or every row with
positive weight); `psychnet_constant_variable` (warning) when a constant
column is dropped; `psychnet_separation` (warning) when items are
perfectly predicted by the others, so their estimates do not exist;
`psychnet_no_converge` (warning) when a regression fails to converge for
another reason; `psychnet_singular_information` (warning) when an item's
Fisher information is singular, so its Wald p-values are unreliable.

## References

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
constructing networks from binary data. *Scientific Reports*, 4, 5918.
[doi:10.1038/srep05918](https://doi.org/10.1038/srep05918)

## See also

[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
for the penalized model;
[`dichotomize()`](https://pak.dynasite.org/psychnets/reference/dichotomize.md)
to prepare items;
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
and
[`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
for the optimality check;
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md);
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `method = "ising_sampler"`.

## Examples

``` r
b <- dichotomize(SRL_GPT)
ising_sampler(b)                 # every edge; certificate near zero
#> <psychnet> ising_sampler network
#>   nodes: 5   edges: 10   (undirected)
#>   optimality (KKT residual): 1.07e-08
pruned <- ising_sampler(b, alpha = 0.05)
as.data.frame(pruned)
#>   from to    weight
#> 1  CSU IV 2.4336970
#> 2  CSU SE 1.9776719
#> 3   IV SE 1.5285878
#> 4  CSU SR 1.7964616
#> 5   IV SR 2.0903313
#> 6   SE SR 0.9333181
certificate(pruned)              # missing: pruning changed the fit
#>          method certificate kind certified
#> 1 ising_sampler          NA  kkt        NA
```
