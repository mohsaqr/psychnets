# Gaussian graphical models: an extended tutorial

## What a psychological network represents

A psychological network is a statistical representation of relations
among measured variables. In a Gaussian graphical model the nodes are
variables: symptoms, questionnaire items, scale scores, behaviours, or
other measured constructs. The edges are partial correlations. An edge
states that two variables remain associated after conditioning on all
other variables in the network, a conditional dependence. The absence of
an edge states that two variables are conditionally independent given
the rest, so any observed bivariate relation between them is accounted
for through other variables in the model.

This extended tutorial assumes the reader has met the model in the
companion vignette. It concentrates on interpretation and on the
questions a careful analysis must answer past the point estimate: how
accurate the edge weights are, how stable the centrality ordering is,
and how the choice of estimation engine affects the result. The
`SRL_GPT` data set runs through the whole tutorial. It contains 300
observations of five Motivated Strategies for Learning Questionnaire
construct scores: cognitive strategy use (`CSU`), intrinsic value
(`IV`), self-efficacy (`SE`), self-regulation (`SR`), and test anxiety
(`TA`), each a mean score on a 1 to 7 scale.

``` r

head(SRL_GPT)
#>        CSU       IV       SE       SR   TA
#> 1 5.307692 5.666667 5.777778 5.333333 4.00
#> 2 5.846154 6.444444 6.000000 5.777778 4.00
#> 3 6.615385 6.666667 6.222222 6.333333 3.25
#> 4 5.692308 6.555556 6.333333 5.555556 4.50
#> 5 4.384615 5.555556 4.888889 4.777778 4.00
#> 6 4.846154 5.444444 5.666667 5.111111 3.50
```

## Fitting and reading the network

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
takes a numeric table and a `method` name and returns a fitted
`psychnet` object. Here `method = "glasso"` requests the
EBIC-regularized graphical lasso, which selects an L1 penalty by the
extended Bayesian information criterion and refits it to the certified
optimum.

``` r

net <- psychnet(SRL_GPT, method = "glasso")
net
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 2.21e-10
```

The network has five nodes and ten edges at a selected penalty of
`lambda` = 0.0086 and the default EBIC hyperparameter `gamma` = 0.5. The
edge table is read with
[`summary()`](https://rdrr.io/r/base/summary.html), which also reports
the range and mean of the weights.

``` r

summary(net)
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 2.21e-10
#>   edge weight: range [-0.350, 0.412], mean 0.174
```

Each row is a conditional association between two constructs. A positive
weight means that respondents higher on one construct tend to be higher
on the other after conditioning on the remaining constructs; a negative
weight means the opposite conditional relation. The learning and
motivation constructs form a positive cluster: `CSU` connects to `IV`
(0.41), `SE` (0.38), and `SR` (0.35), and `IV`, `SE`, and `SR` are
positively linked among themselves. The one substantial negative
association is between `SR` and `TA` (-0.35), which reads as lower test
anxiety among respondents with higher self-regulation, holding the other
constructs constant. The three positive edges involving `TA` are small
(0.05 to 0.11), so their substantive interpretation is cautious. These
weights are conditional relations among construct scores and carry no
causal claim.

[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
summarizes how each node is positioned in the network, returning a tidy
data frame with one row per node and a column per measure. Node strength
is the sum of the absolute edge weights at a node. Expected influence is
the signed sum, so positive and negative edges offset one another.

``` r

net_centralities(net)
#>   node  strength expected_influence
#> 1  CSU 1.1984512         1.19845117
#> 2   IV 1.0012765         1.00127649
#> 3   SE 0.8492185         0.84921854
#> 4   SR 1.2314317         0.53172852
#> 5   TA 0.6082238        -0.09147935
```

`SR` has high strength (1.23) because it connects to several learning
constructs and carries the negative edge with `TA`. Its expected
influence (0.53) is lower than its strength because the negative edge
offsets its positive connections. Strength is usually the most
interpretable and comparatively stable centrality index in psychological
networks. Betweenness and closeness are often unstable in these networks
and are available on request but are read with more caution.

[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
reports how well each node is predicted by its neighbours, returning a
tidy data frame with columns `node`, `type`, `metric`, `predictability`,
and `accuracy`. For a Gaussian graphical model the value is the
closed-form proportion of variance of each node explained by the rest of
the network, computed from the precision matrix.

``` r

net_predict(net)
#>   node     type metric predictability accuracy
#> 1  CSU gaussian     R2      0.8217840       NA
#> 2   IV gaussian     R2      0.7690346       NA
#> 3   SE gaussian     R2      0.7142811       NA
#> 4   SR gaussian     R2      0.7786477       NA
#> 5   TA gaussian     R2      0.1421163       NA
```

Predictability is an absolute complement to centrality: two nodes can
have similar centrality yet differ in how well the network explains
them. The learning constructs are highly predictable (`CSU` 0.82, `SR`
0.78, `IV` 0.77, `SE` 0.71), whereas `TA` is much less predictable
(0.14). That pattern fits the edge table: test anxiety has one
substantial edge and three weak ones.

``` r

cograph::splot(net, psych_styling = TRUE)
```

![](gaussian-graphical-models-v2_files/figure-html/ggm-plot-1.png)

The graph is a visual summary of the edge table and is useful for
orientation. The numerical edge table remains the primary object for
interpretation.

## Assessing accuracy of edge weights

A reported edge weight is one estimate from one sample. Before an edge
is read as substantive, its sampling variability is assessed.
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
resamples the observations with replacement, re-estimates the network on
each resample, and summarizes the sampling distribution of every edge
weight and centrality. It returns an object of class
`psychnet_bootstrap` whose edge table carries the observed weight, the
bootstrap mean, the lower and upper percentile confidence limits, the
inclusion proportion, and a `significant` flag set when the percentile
interval excludes zero. The following run uses 200 resamples for a
compact vignette; a full analysis uses more.

``` r

set.seed(1)
bs <- net_boot(SRL_GPT, method = "glasso", n_boot = 200, cores = 1)
bs
#> <psychnet_bootstrap> glasso, 200 resamples, 95% CI
#>   10 edges (7 significant), 5 nodes, measures: strength, expected_influence
```

The printed summary reports that seven of the ten edges are significant
at the 95% level. The full edge table gives the intervals.

``` r

summary(bs)
#>                 Length Class      Mode     
#> observed          19   psychnet   list     
#> edges              8   data.frame list     
#> centrality         7   data.frame list     
#> edge_boot       2000   -none-     numeric  
#> str_boot        1000   -none-     numeric  
#> ei_boot         1000   -none-     numeric  
#> centrality_boot    2   -none-     list     
#> edge_labels       10   -none-     character
#> node_labels        5   -none-     character
#> measures           2   -none-     character
#> n_boot             1   -none-     numeric  
#> ci                 1   -none-     numeric  
#> method             1   -none-     character
#> lambda_path      100   -none-     numeric  
#> lambda_selected    1   -none-     numeric
```

A percentile confidence interval is the range of edge values consistent
with the data at the stated level. An interval that excludes zero
indicates a non-zero edge. The strong edges of the positive cluster have
narrow intervals well away from zero: `CSU`-`IV` runs from 0.32 to 0.49,
and the negative `SR`-`TA` edge runs from -0.43 to -0.24. The three weak
positive `TA` edges have intervals whose lower limit is zero, so they
are not significant, and their inclusion proportions are below one
(`CSU`-`TA` at 0.72, `SE`-`TA` at 0.96, `IV`-`TA` at 0.97), meaning the
edge was set to zero in a share of the resamples. The bootstrap
therefore gives the same reading as the point estimate but with the
added evidence that the weak `TA` edges are uncertain.

## Assessing stability of centrality

Centrality indices are functions of many estimated edges, so their
ordering also needs a stability check.
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
drops random subsets of cases, re-estimates the network on each subset,
and correlates the subset centralities with the full-sample
centralities. It returns an object of class `psychnet_stability`
carrying the correlation-stability (CS) coefficient per measure and a
tidy table of the mean correlation at each drop proportion. The CS
coefficient is the largest proportion of cases that can be dropped while
the rank correlation with the full sample stays at or above a threshold
(default 0.7) with a stated certainty (default 0.95).

``` r

set.seed(1)
cs <- net_stability(SRL_GPT, method = "glasso",
                    drop_prop = c(0.3, 0.5, 0.7), iter = 50)
cs
#> <psychnet_stability> glasso, 50 subsets/proportion
#>   CS-coefficient (cor >= 0.70 with 95% certainty):
#>     strength             0.50
#>     expected_influence   0.70
```

The CS coefficient is 0.50 for strength and 0.70 for expected influence.
The strength ordering stays correlated at 0.7 or above with the
full-sample ordering, with 95% certainty, when up to 50% of the cases
are dropped; the expected-influence ordering holds up to 70%. A common
guideline treats a CS coefficient of at least 0.5 as adequate and 0.25
as a minimum, so both measures are stable enough to interpret here. The
supporting table shows how the mean correlation declines as more cases
are dropped.

``` r

cs$table
#>              measure drop_prop mean_cor     sd_cor prop_above
#> 1           strength       0.3    0.970 0.05050763       1.00
#> 2           strength       0.5    0.940 0.05714286       1.00
#> 3           strength       0.7    0.888 0.07182746       0.94
#> 4 expected_influence       0.3    0.998 0.01414214       1.00
#> 5 expected_influence       0.5    0.992 0.02740475       1.00
#> 6 expected_influence       0.7    0.964 0.06928203       0.96
```

For strength the mean correlation falls from 0.97 at a 30% drop to 0.89
at a 70% drop; for expected influence it falls from 0.998 to 0.96 over
the same range. The higher correlations for expected influence match its
larger CS coefficient.

## Native and external estimation engines

By default (`native = TRUE`) the graphical lasso is estimated with the
package’s own solver, written entirely in base R. This is a clean-room
implementation: the algorithm is reconstructed from first principles and
reproduces the estimate without wrapping an external library. It depends
only on base R, so it installs and runs anywhere with no compilation; it
is readable R code that can be inspected end to end; and each fit is
graded against the convex objective it optimizes, so its correctness is
read from the certificate.

Setting `native = FALSE` delegates each penalized solve to the compiled
`glasso` Fortran routine, the same code underlying
[`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html),
which is faster on large problems. The optional `glasso` package
supplies it.

``` r

net_fortran <- psychnet(SRL_GPT, method = "glasso", native = FALSE)
summary(net_fortran)
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 8.19e-05
#>   edge weight: range [-0.350, 0.412], mean 0.174
```

The two solvers estimate the same model and return the same network: the
edge weights agree with the native fit to about five decimal places
(`CSU`-`IV` is 0.4117 under both). The certificate reads the difference
between them.

``` r

certificate(net_fortran)
#>   method  certificate kind certified
#> 1 glasso 8.193488e-05  kkt     FALSE
```

The residual here is 8.2e-05, above the default tolerance of 1e-6, so
`certified` is `FALSE`. The Fortran routine converges to its own looser
tolerance (near 1e-4), and the certificate reports that looser
tolerance; the base-R refit reaches a tighter one. The native fit
certifies at 2.2e-10, so both solvers recover the same network while the
base-R solver returns it to a tighter numerical optimum. A study can use
the native solver and read the external one as an independent check.

## Non-normal and ordinal data

Gaussian graphical models are Gaussian, yet psychological data are often
skewed, bounded, and measured on ordinal response scales. Scale means
such as the five `SRL_GPT` construct scores are commonly treated as
approximately continuous, which is an approximation of the questionnaire
response format.

For ordinal item-level data, `cor_method = "auto"` bases the correlation
input on polychoric or polyserial associations when appropriate. This
choice is most relevant when the variables are individual Likert items
or mixed ordinal-continuous measures.

For continuous variables with pronounced monotone departures from
normality, `method = "huge"` estimates a nonparanormal graphical model.
The nonparanormal assumes the variables become jointly Gaussian after
unknown monotone transformations of their marginal distributions, then
applies the identical EBIC graphical lasso to the transformed
correlation matrix. It retains the conditional-dependence interpretation
of the Gaussian graphical model while relaxing strict multivariate
normality.

``` r

huge_net <- psychnet(SRL_GPT, method = "huge")
certificate(huge_net)
#>   method  certificate kind certified
#> 1   huge 2.471392e-10  kkt      TRUE
```

The nonparanormal fit certifies at 2.5e-10, at the order of machine
precision, so it reached the optimum of its convex objective on the
transformed correlation.

## Mathematical foundations

This section gives the statistical definitions behind the tutorial. It
can be skipped without losing the worked example.

### Precision matrix and the pairwise Markov property

Let $`\mathbf{X} = (X_1, \ldots, X_p)`$ follow a multivariate normal
distribution with covariance matrix $`\boldsymbol{\Sigma}`$ and
precision matrix $`\mathbf{K} = \boldsymbol{\Sigma}^{-1}`$. The model
represents conditional independences by zeros in the precision matrix.
For $`i \neq j`$,

``` math
K_{ij} = 0
\quad \Longleftrightarrow \quad
X_i \perp\!\!\!\perp X_j \mid \mathbf{X}_{-(i,j)} ,
```

so an absent edge is a conditional independence given all remaining
variables.

### Partial-correlation identity

The edge weights are reported as partial correlations. The partial
correlation between variables $`i`$ and $`j`$, conditional on all other
variables, is

``` math
\rho_{ij \cdot \mathrm{rest}}
= -\frac{K_{ij}}{\sqrt{K_{ii}K_{jj}}} ,
```

which puts the conditional association on the correlation scale: a
positive value is a positive conditional association, a negative value a
negative one, and zero is no edge.

### The L1-penalized Gaussian log-likelihood

Let $`\mathbf{S}`$ be the sample correlation matrix. The graphical lasso
estimates the precision matrix by maximizing the penalized Gaussian
log-likelihood,

``` math
\hat{\mathbf{K}}
= \arg\max_{\mathbf{K} \succ 0}
\left[
\log\det(\mathbf{K})
- \operatorname{tr}(\mathbf{S}\mathbf{K})
- \lambda \sum_{i \neq j} |K_{ij}|
\right] .
```

The penalty $`\lambda \geq 0`$ controls sparsity: larger values shrink
more off-diagonal entries to zero and yield fewer edges. The objective
is strictly concave in $`\mathbf{K}`$, so its maximizer is unique.

### EBIC with gamma

The extended Bayesian information criterion selects among candidate
sparse graphs by balancing fit against complexity. For a graph with
$`E`$ non-zero edges from $`n`$ observations and $`p`$ variables,

``` math
\mathrm{EBIC}_{\gamma}
= -2\ell(\hat{\mathbf{K}})
+ E \log(n)
+ 4 E \gamma \log(p) ,
```

where $`\ell(\hat{\mathbf{K}})`$ is the maximized Gaussian
log-likelihood and $`\gamma`$ controls the additional penalty for the
size of the graph space. When $`\gamma = 0`$ the criterion reduces to
ordinary BIC; larger $`\gamma`$ favours sparser networks.

### Bootstrap and case-dropping stability

The nonparametric bootstrap draws resamples of the observations with
replacement, re-estimates the network on each, and forms the percentile
confidence interval of each edge weight and centrality from the
resampling distribution. The case-dropping procedure removes random
fractions of the cases, re-estimates on each subset, and correlates the
subset centralities with the full-sample centralities. The
correlation-stability coefficient is the largest drop proportion at
which that correlation stays at or above a threshold with a stated
certainty. Both procedures re-estimate through the same
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
path, so they apply to any estimator with no method-specific code.
