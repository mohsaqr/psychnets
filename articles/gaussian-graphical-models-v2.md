# Gaussian graphical models: an extended tutorial

``` r

library(psychnets)
```

## What a psychological network represents

A psychological network is a statistical representation of relations
among measured variables. In a Gaussian graphical model (GGM), the nodes
are variables: symptoms, questionnaire items, scale scores, behaviours,
or other measured constructs. The edges are partial correlations. An
edge therefore does not mean that two variables are merely correlated in
the ordinary bivariate sense. It means that they remain associated after
conditioning on all other variables included in the network.

This distinction is central to the interpretation of psychological
networks (Epskamp & Fried, 2018; Borsboom et al., 2021; Saqr, Beck &
López-Pernas, 2024). Suppose that cognitive strategy use and
self-regulation are connected in a network that also contains intrinsic
value, self-efficacy, and test anxiety. The edge asks a ceteris paribus
question: are cognitive strategy use and self-regulation still
associated among respondents who are comparable on the other constructs
in the model? If so, the edge represents a conditional dependence. If
not, the pair is conditionally independent in the fitted graph.

The absence of an edge is consequently informative. It says that, given
the other variables in the network, the data do not support a direct
conditional association between the two nodes under the fitted model.
This does not prove that the constructs are unrelated in every
substantive sense. It means that any observed zero-order relation
between them can be accounted for through other variables in the
network, sampling variation, or both. Interpreting both present and
absent edges is what makes a GGM more than a visual display of ordinary
correlations.

## Why and how we regularize

An unregularized partial-correlation network is usually dense. When the
sample correlation matrix is inverted directly, almost every pair of
variables receives a non-zero partial correlation. Some of those partial
correlations may reflect stable conditional associations, but many small
estimates will arise because sample correlations fluctuate. The problem
becomes more serious as the number of variables grows relative to the
number of observations, because each edge is estimated while adjusting
for many other variables.

Regularization addresses this by shrinking weak edges. The graphical
lasso (Friedman, Hastie & Tibshirani, 2008) estimates a sparse precision
matrix by penalizing the absolute size of off-diagonal entries. In the
partial-correlation network, this has the practical consequence that
some weak conditional associations are set exactly to zero. The result
is a graph that is easier to read and less dominated by small
sample-specific estimates.

Regularization is not a guarantee of truth. It changes the error
tradeoff. A sparser graph typically reduces false-positive edges, but it
may also remove weak associations that are real in the population. The
substantive question is therefore not simply whether an edge appears,
but how strong and stable the edge is, how plausible the model
assumptions are, and how the result changes under reasonable analytic
choices.

The extended Bayesian information criterion (EBIC; Foygel & Drton, 2010)
is commonly used to choose the amount of regularization. Its
hyperparameter, usually called gamma, controls an additional penalty for
model complexity. A larger gamma favours sparser networks; a smaller
gamma is more permissive. Choosing gamma is therefore a modelling
decision about conservatism. In psychological-network applications,
gamma equal to 0.5 is often used as a default because it gives a
relatively cautious graph, but researchers should still treat weak and
borderline edges as uncertain until accuracy and stability have been
examined.

## A worked analysis with SRL_GPT

The `SRL_GPT` data set contains 300 observations of five Motivated
Strategies for Learning Questionnaire construct scores: cognitive
strategy use (`CSU`), intrinsic value (`IV`), self-efficacy (`SE`),
self-regulation (`SR`), and test anxiety (`TA`). Each variable is a mean
score on a 1-7 response scale.

``` r

head(x = SRL_GPT)
#>        CSU       IV       SE       SR   TA
#> 1 5.307692 5.666667 5.777778 5.333333 4.00
#> 2 5.846154 6.444444 6.000000 5.777778 4.00
#> 3 6.615385 6.666667 6.222222 6.333333 3.25
#> 4 5.692308 6.555556 6.333333 5.555556 4.50
#> 5 4.384615 5.555556 4.888889 4.777778 4.00
#> 6 4.846154 5.444444 5.666667 5.111111 3.50
```

The network is estimated with
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md).
Here `method = "glasso"` requests the EBIC-regularized graphical lasso.
The method name `"EBICglasso"` is accepted as an alias, but this
vignette uses the shorter spelling.

``` r

psychnet(data = SRL_GPT, method = "glasso")
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 2.21e-10
```

The edge list is the most useful first output because it states the
fitted partial correlations directly.

``` r

as.data.frame(x = net)
#>    from to      weight
#> 1   CSU IV  0.41166866
#> 2   CSU SE  0.38290223
#> 3    IV SE  0.15992605
#> 4   CSU SR  0.35481359
#> 5    IV SR  0.31909443
#> 6    SE SR  0.20767208
#> 7   CSU TA  0.04906668
#> 8    IV TA  0.11058735
#> 9    SE TA  0.09871818
#> 10   SR TA -0.34985157
```

Each row is a conditional association between two constructs. Positive
weights mean that respondents who are higher on one construct tend also
to be higher on the other, after conditioning on the remaining
constructs. Negative weights mean that respondents who are higher on one
construct tend to be lower on the other, again conditionally on the rest
of the network.

In this example, the learning-strategy and motivation constructs form a
positive cluster. `CSU` is positively connected to `IV`, `SE`, and `SR`;
`IV`, `SE`, and `SR` are also positively linked. These edges should be
read as conditional relations among construct scores, not as evidence
that one construct causes another. The strongest negative association is
between `SR` and `TA`, suggesting that test anxiety is lower among
respondents with higher self-regulation when cognitive strategy use,
intrinsic value, and self-efficacy are held constant. Several positive
edges involving `TA` are small; their substantive interpretation should
be correspondingly cautious.

Centrality summaries describe how a node is positioned in the fitted
network.

``` r

net_centralities(x = net)
#>   node  strength expected_influence
#> 1  CSU 1.1984512         1.19845117
#> 2   IV 1.0012765         1.00127649
#> 3   SE 0.8492185         0.84921854
#> 4   SR 1.2314317         0.53172852
#> 5   TA 0.6082238        -0.09147935
```

Node strength is the sum of the absolute edge weights incident on a
node. It is usually the most interpretable and comparatively stable
centrality index in psychological networks. Expected influence is the
signed analogue: positive and negative edges can offset one another. In
this data set, `SR` has high strength because it is connected to several
learning constructs and also has a sizeable negative edge with `TA`. Its
expected influence is lower than its strength because the negative edge
offsets its positive connections.

Other centrality indices, especially betweenness and closeness, are
often unstable in psychological networks and should not be
overinterpreted without dedicated stability evidence (Bringmann et al.,
2019). For this reason, this vignette emphasizes strength and expected
influence rather than treating every available graph-theoretic index as
equally meaningful.

Predictability provides an absolute complement to centrality. Instead of
asking how connected a node is, it asks how much of the variance in each
node can be accounted for by its neighbours in the network (Haslbeck &
Waldorp, 2018).

``` r

net_predict(x = net)
#>   node     type metric predictability accuracy
#> 1  CSU gaussian     R2      0.8217840       NA
#> 2   IV gaussian     R2      0.7690346       NA
#> 3   SE gaussian     R2      0.7142811       NA
#> 4   SR gaussian     R2      0.7786477       NA
#> 5   TA gaussian     R2      0.1421163       NA
```

Predictability is useful because two nodes can have similar centrality
but differ in how well the rest of the network explains them. Here the
learning-related constructs are highly predictable from their
neighbours, whereas `TA` is much less predictable. That pattern fits the
edge list: test anxiety has one substantial negative edge with
self-regulation and several smaller positive edges.

The same fitted object can be plotted with `cograph`. The chunk is
guarded so the vignette can still render on systems where the optional
plotting package is not available.

``` r

cograph::splot(x = net, psych_styling = TRUE)
```

![](gaussian-graphical-models-v2_files/figure-html/ggm-plot-1.png)

The graph is a visual summary of the edge list. It is useful for
orientation, but the numerical edge table remains the primary object for
interpretation. Layout algorithms place nodes to make the graph
readable; distances on the page should not be interpreted as direct
measurements unless the plotting method explicitly defines them that
way.

## Assessing accuracy and replicability

Point estimates are not enough. A reported edge weight is one estimate
from one sample, and centrality indices are functions of many estimated
edges. Researchers should therefore assess how much those estimates vary
under resampling.

For edge weights, nonparametric bootstrapping can be used to obtain
uncertainty intervals and to examine whether differences between edges
are large relative to sampling variability. For centrality,
case-dropping procedures ask whether the centrality ordering remains
similar when portions of the sample are removed. This is the logic
behind the bootstrap and stability procedures recommended by Epskamp,
Borsboom, and Fried (2018).

In `psychnets`,
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
and
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
provide these analyses. They are not run here because a thorough
bootstrap is computationally longer than a compact vignette should
require, but they should be part of a complete empirical analysis when
edge accuracy or centrality claims are substantively important.

## Non-normal and ordinal data

GGMs are Gaussian models, yet psychological data are often skewed,
bounded, and measured on ordinal response scales. Scale means, such as
the five `SRL_GPT` construct scores, are often treated as approximately
continuous, but this is a modelling approximation rather than a property
guaranteed by the questionnaire.

For ordinal item-level data, `cor_method = "auto"` can be used so that
the correlation input is based on polychoric or polyserial associations
when appropriate. This choice is most relevant when variables are
individual Likert items or mixed ordinal-continuous measures rather than
scale means.

For continuous variables with pronounced monotone departures from
normality, `method = "huge"` estimates a nonparanormal graphical model.
The nonparanormal approach assumes that variables become jointly
Gaussian after unknown monotone transformations of their marginal
distributions (Liu, Lafferty & Wasserman, 2009). It is a way to retain a
GGM-style conditional-dependence interpretation while relaxing strict
multivariate normality.

## Native estimation

By default (`native = TRUE`) `psychnets` estimates the graphical lasso
with its own solver, written entirely in base R. This is a *clean-room*
implementation: the algorithm is reconstructed from the original
description (Friedman, Hastie & Tibshirani, 2008) rather than wrapping
an existing library, and it reproduces the estimate from first
principles. Its advantages are:

- it depends only on base R, so it installs and runs anywhere with no
  compilation or external library;
- the estimator is plain, readable R code that can be inspected and
  audited end to end, not a compiled black box;
- each fit is graded against the convex objective it optimizes, so its
  correctness is established internally rather than by agreement with
  another program;
- the results are reproducible across platforms.

## The Fortran glasso estimator

The native solver is not the only option. Setting `native = FALSE`
delegates each penalized solve to the compiled `glasso` Fortran routine
— the same code underlying
[`qgraph::EBICglasso()`](https://rdrr.io/pkg/qgraph/man/EBICglasso.html)
— which is faster on large problems and returns numerically equivalent
estimates:

``` r

psychnet(data = SRL_GPT, method = "glasso", native = FALSE)
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 8.19e-05
```

The native (`native = TRUE`) and Fortran (`native = FALSE`) solvers
estimate the same model; they differ only in the implementation that
performs the computation, so a study can use the convenient native
solver and confirm it against the established external one.

## Mathematical foundations

This final section gives the statistical definitions behind the
tutorial. It is included for readers who want the mathematical form of
the model and estimator, and can be skipped without losing the worked
example above.

### Precision matrix and the pairwise Markov property

Let $`\mathbf{X} = (X_1, \ldots, X_p)`$ follow a multivariate normal
distribution with covariance matrix $`\boldsymbol{\Sigma}`$ and
precision matrix $`\mathbf{K} = \boldsymbol{\Sigma}^{-1}`$. A Gaussian
graphical model represents conditional independences among the variables
by zeros in the precision matrix. For $`i \neq j`$, the pairwise Markov
property states that

``` math
K_{ij} = 0
\quad \Longleftrightarrow \quad
X_i \perp\!\!\!\perp X_j \mid \mathbf{X}_{-(i,j)} .
```

Thus an absent edge corresponds to conditional independence given all
remaining variables in the model (Lauritzen, 1996).

### Partial-correlation identity

The edge weights in a GGM are usually reported as partial correlations
rather than raw precision-matrix entries. The partial correlation
between variables $`i`$ and $`j`$, conditional on all other variables,
is

``` math
\rho_{ij \cdot \mathrm{rest}}
= -\frac{K_{ij}}{\sqrt{K_{ii}K_{jj}}}.
```

This transformation puts the conditional association on the familiar
correlation scale. A positive value indicates a positive conditional
association, a negative value indicates a negative conditional
association, and zero corresponds to no edge in the graph.

### The L1-penalized Gaussian log-likelihood objective

Let $`\mathbf{S}`$ be the sample covariance or correlation matrix. The
graphical lasso estimates the precision matrix by maximizing the
penalized Gaussian log-likelihood

``` math
\hat{\mathbf{K}}
= \arg\max_{\mathbf{K} \succ 0}
\left[
\log\det(\mathbf{K})
- \operatorname{tr}(\mathbf{S}\mathbf{K})
- \lambda \sum_{i \neq j} |K_{ij}|
\right].
```

The penalty parameter $`\lambda \geq 0`$ controls sparsity. Larger
values shrink more off-diagonal precision entries to zero, producing
fewer edges; smaller values retain more edges. This is the statistical
source of the sparse network estimated by the graphical lasso (Friedman,
Hastie & Tibshirani, 2008).

### EBIC with gamma

The EBIC selects among candidate sparse graphs by balancing model fit
against model complexity. For a graph with $`E`$ non-zero edges
estimated from $`n`$ observations and $`p`$ variables, a common form is

``` math
\mathrm{EBIC}_{\gamma}
= -2\ell(\hat{\mathbf{K}})
+ E \log(n)
+ 4 E \gamma \log(p),
```

where $`\ell(\hat{\mathbf{K}})`$ is the maximized Gaussian
log-likelihood and $`\gamma`$ controls the additional penalty for the
size of the graph space (Foygel & Drton, 2010). When $`\gamma = 0`$,
EBIC reduces to the ordinary BIC penalty for the number of edges. Larger
values of $`\gamma`$ favour sparser networks.

## References

Borsboom, D., Deserno, M. K., Rhemtulla, M., Epskamp, S., Fried, E. I.,
McNally, R. J., Robinaugh, D. J., Perugini, M., Dalege, J., Costantini,
G., Isvoranu, A.-M., Wysocki, A. C., van Borkulo, C. D., van Bork, R., &
Waldorp, L. J. (2021). Network analysis of multivariate data in
psychological science. *Nature Reviews Methods Primers*, 1, Article 58.

Bringmann, L. F., Elmer, T., Epskamp, S., Krause, R. W., Schoch, D.,
Wichers, M., Wigman, J. T. W., & Snippe, E. (2019). What do centrality
measures measure in psychological networks? *Journal of Abnormal
Psychology*, 128(8), 892-903.

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods*, 50(1), 195-212.

Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
correlation networks. *Psychological Methods*, 23(4), 617-634.

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. In *Advances in Neural Information
Processing Systems* (Vol. 23, pp. 604-612).

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853-861.

Lauritzen, S. L. (1996). *Graphical Models*. Oxford University Press.

Liu, H., Lafferty, J., & Wasserman, L. (2009). The nonparanormal:
Semiparametric estimation of high dimensional undirected graphs.
*Journal of Machine Learning Research*, 10, 2295-2328.

Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological networks: A
modern approach to the analysis of learning and complex learning
processes. In M. Saqr & S. López-Pernas (Eds.), *Learning Analytics
Methods and Tutorials: A Practical Guide Using R* (Chapter 19).
Springer.
<https://lamethods.org/book1/chapters/ch19-psychological-networks/ch19-psych.html>
