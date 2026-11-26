# Estimating and interpreting Gaussian graphical models with the psychnets R package

Does the association between motivation and achievement persist once
effort, prior achievement and other measured constructs are taken into
account? This question arises when psychological and educational
researchers examine interrelated characteristics such as self-efficacy,
learning strategies, self-regulation and anxiety. Pairwise correlations
may reflect associations shared with other scores, making them difficult
to interpret. Network analysis examines such interdependence without
assigning every construct a fixed role as predictor or outcome (Borsboom
et al., 2021). A Gaussian graphical model (GGM) represents variables as
nodes and conditional associations as edges, estimating partial
correlations between pairs of variables while accounting for the
remaining variables (Epskamp & Fried, 2018). An edge does not establish
causality, and a missing edge does not prove that two constructs are
unrelated.

GGMs are useful when researchers want to distinguish associations that
remain after adjustment from those visible in an ordinary correlation
matrix. For example, two learning strategies may be strongly correlated
but have a much weaker partial correlation once motivation and
self-regulation are included. Regularisation can set small estimated
edges to zero and make a network easier to interpret, although the
selected model need not be sparse. The resulting pattern depends on the
variables, measurement scale, sample and estimation procedure.
Researchers therefore need to consider the accuracy of estimated edges
and the consequences of methodological choices, rather than interpreting
a network diagram on its own (Epskamp et al., 2018; Isvoranu & Epskamp,
2023; Williams & Rast, 2020).

Psychological networks have also been applied to educational questions
about how learning processes operate together. Studies have examined
associations among monitoring events during self-regulated learning,
including how those events relate to phases of regulation (Malmberg,
Saqr, Järvenoja, Haataja, et al., 2022; Malmberg, Saqr, Järvenoja, &
Järvelä, 2022). Network approaches have likewise been used to study the
relational and temporal aspects of self-regulation in academic writing
(Saqr et al., 2021). These applications illustrate why the relationships
among learning processes may be of interest alongside the levels of
individual constructs (Saqr et al., 2024). The choice of network model
must nevertheless match the research design: a cross-sectional GGM does
not describe within-person change over time.

This article provides a worked example of a regularised GGM using five
self-regulated learning construct scores from language-model-generated
survey responses. It explains how to estimate and interpret
partial-correlation edges, compare them with ordinary correlations,
assess computational optimality, and describe node strength, expected
influence and predictability. It then compares Pearson and Spearman
estimation, applies a post-hoc edge threshold and fits a nonparanormal
model. The aim is to support careful interpretation and reporting by
applied researchers, not to draw conclusions about students from these
generated responses.

``` r

library(psychnets)
```

## The model

For variables $`X_1, \ldots, X_p`$ with covariance matrix
$`\boldsymbol\Sigma`$, the inverse covariance matrix, or precision
matrix, $`\mathbf{K} = \boldsymbol\Sigma^{-1}`$, describes their
conditional association structure. The partial correlation between
$`X_i`$ and $`X_j`$ given the remaining variables is

``` math
\rho_{ij \cdot \mathrm{rest}} = -\frac{K_{ij}}{\sqrt{K_{ii}K_{jj}}},
```

Under multivariate normality, $`K_{ij} = 0`$ means that $`X_i`$ and
$`X_j`$ are conditionally independent given the remaining variables. A
non-zero partial correlation is represented by an edge. Its sign
indicates the direction of the conditional association and its magnitude
indicates its strength on the partial-correlation scale.

Sample partial correlations are subject to sampling variation, and small
non-zero estimates may not reflect a reproducible association. The
graphical lasso (Friedman et al., 2008) estimates the precision matrix
while penalising the absolute values of its off-diagonal entries,

``` math
\widehat{\mathbf{K}} = \underset{\mathbf{K} \succ 0}{\operatorname{argmin}}
\left\{ \operatorname{tr}(\mathbf{S}\mathbf{K}) - \log\det(\mathbf{K})
+ \lambda \sum_{i \ne j}|K_{ij}| \right\},
```

where $`\mathbf{S}`$ is the sample correlation matrix and $`\lambda`$
controls the amount of shrinkage. Larger penalties can set more
off-diagonal entries to zero, removing their corresponding edges. The
penalty $`\lambda`$ is selected from candidate values using the extended
Bayesian information criterion (EBIC; Foygel & Drton, 2010),

``` math
\operatorname{EBIC}_{\gamma} = -2\ell(\widehat{\mathbf{K}}) + E\log n
+ 4E\gamma\log p,
```

Here $`\ell`$ is the Gaussian log-likelihood, $`E`$ is the number of
edges, $`n`$ is the sample size and $`p`$ is the number of nodes. The
hyperparameter $`\gamma`$ (0.5 by default) adds a complexity penalty
beyond the usual BIC term. Larger values favour sparser networks, but
the selected number of edges still depends on the data.

## Data

`SRL_GPT` contains 300 sets of responses to the Motivated Strategies for
Learning Questionnaire. Each variable is a construct score calculated as
the mean of its items on a scale from 1 to 7. The responses were
generated by a large language model (GPT) for a study of how these
models complete psychological questionnaires (Vogelsmeier et al., 2025).
They provide an illustrative dataset for learning the method and must
not be interpreted as observations of students.

- **Nodes:** five construct scores, each treated as continuous:
  - `CSU`, cognitive strategy use;
  - `IV`, intrinsic value;
  - `SE`, self-efficacy;
  - `SR`, self-regulation;
  - `TA`, test anxiety.
- **Observations:** one row per respondent; the data are complete.

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

Because the construct scores are averages of several items and have many
distinct values, treating them as continuous is reasonable for this
example. For analyses of item-level ordinal responses, polychoric
correlations may be more appropriate (`cor_method = "auto"`).

## Which constructs are directly associated?

[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
with `method = "glasso"` estimates the EBIC graphical lasso. The
defaults are Pearson correlations, pairwise-complete handling of missing
values and $`\gamma = 0.5`$.

``` r

net <- psychnet(data = SRL_GPT, method = "glasso")
summary(net)
#> <psychnet> glasso network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality (KKT residual): 2.21e-10
#>   edge weight: range [-0.350, 0.412], mean 0.174
```

The selected penalty is $`\lambda = 0.00861`$. The fitted network
retains all 10 possible edges among the five nodes, so regularisation
has not removed any pairwise connection in this example. The partial
correlations range from -0.350 to 0.412. The edge table gives the
estimated weight of each pair:

``` r

as.data.frame(net)
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

The four learning constructs are joined by positive edges. The largest
positive partial correlation is between cognitive strategy use and
intrinsic value (0.412), followed by cognitive strategy use and
self-efficacy (0.383), and cognitive strategy use and self-regulation
(0.355). These estimates indicate positive associations after adjustment
for the other three constructs, rather than unadjusted pairwise
relationships. Self-regulation and test anxiety have the only negative
edge (-0.350), indicating that higher self-regulation is associated with
lower test anxiety when the other variables are held constant. The
remaining edges involving test anxiety are small and positive: `CSU`
(0.049), `SE` (0.099) and `IV` (0.111). Their substantive importance
should not be inferred from the presence of an edge alone.

## How do partial correlations differ from correlations?

A marginal correlation describes the association between two variables
without controlling for the others. The correlation network estimated
with `method = "cor"` provides a useful comparison with the adjusted
associations in the GGM.

``` r

as.data.frame(psychnet(data = SRL_GPT, method = "cor"))
#>    from to      weight
#> 1   CSU IV  0.86099753
#> 2   CSU SE  0.83458495
#> 3    IV SE  0.78460033
#> 4   CSU SR  0.84743518
#> 5    IV SR  0.82163456
#> 6    SE SR  0.77913881
#> 7   CSU TA -0.13851086
#> 8    IV TA -0.11372136
#> 9    SE TA -0.09672104
#> 10   SR TA -0.30196344
```

The six correlations among the four learning constructs range from 0.779
to 0.861, considerably larger than their partial correlations, which
range from 0.160 to 0.412. This reduction indicates that their pairwise
relationships overlap with associations involving the other constructs.
Test anxiety has a negative marginal correlation with each learning
construct (-0.097 to -0.302), whereas its partial correlations with
`CSU`, `IV` and `SE` are small and positive. Its partial correlation
with self-regulation remains negative (-0.350). Thus, the adjusted
results give a different account of how test anxiety relates to the
constructs. They do not establish that self-regulation causes the other
marginal associations.

## Is the fitted network the optimum of its objective?

The graphical lasso uses numerical optimisation to estimate the network.
Before interpreting the edges, it is useful to check whether the
algorithm has converged sufficiently closely to the solution of its
objective. The Karush–Kuhn–Tucker optimality conditions provide such a
check.
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the largest numerical violation of these conditions and returns
`certified = TRUE` when that value is below $`10^{-6}`$. This is a check
on computational accuracy, not on the reliability of the estimated
associations.

``` r

certificate(net)
#>   method  certificate kind certified
#> 1 glasso 2.213389e-10  kkt      TRUE
```

The reported violation is $`2.21 \times 10^{-10}`$, well below the
stated tolerance. The estimated network therefore satisfies the
numerical optimality check for the selected penalty. This result does
not assess sampling uncertainty, model fit or the plausibility of the
distributional assumptions.

## Which constructs are most connected?

[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
summarizes each node’s edges. Strength is the sum of absolute edge
weights; expected influence is the signed sum, in which negative edges
offset positive ones (Robinaugh et al., 2016).

``` r

net_centralities(net)
#>   node  strength expected_influence
#> 1  CSU 1.1984512         1.19845117
#> 2   IV 1.0012765         1.00127649
#> 3   SE 0.8492185         0.84921854
#> 4   SR 1.2314317         0.53172852
#> 5   TA 0.6082238        -0.09147935
```

Self-regulation has the highest strength (1.231), closely followed by
cognitive strategy use (1.198). Cognitive strategy use has only positive
edges, so its expected influence also equals 1.198. By contrast, the
negative edge between self-regulation and test anxiety reduces the
expected influence of self-regulation to 0.532. Test anxiety has the
lowest strength (0.608) and a slightly negative expected influence
(-0.091). Strength therefore identifies self-regulation as the most
connected node in terms of absolute edge weights, while expected
influence shows why the signs of those edges matter. These are
descriptions of this fitted network, not evidence that changing a highly
connected construct would change the others.

## How well is each construct predicted by the others?

Predictability is the proportion of a node’s variance explained by its
neighbours, $`R^2`$ (Haslbeck & Waldorp, 2018). For a GGM it follows
from the precision matrix, $`R_j^2 = 1 - 1/(K_{jj}S_{jj})`$, so
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
needs no further data.

``` r

net_predict(net)
#>   node     type metric predictability accuracy
#> 1  CSU gaussian     R2      0.8217840       NA
#> 2   IV gaussian     R2      0.7690346       NA
#> 3   SE gaussian     R2      0.7142811       NA
#> 4   SR gaussian     R2      0.7786477       NA
#> 5   TA gaussian     R2      0.1421163       NA
```

In this fitted model, the other nodes account for approximately 82% of
the variance in cognitive strategy use ($`R^2 = 0.822`$), 78% in
self-regulation, 77% in intrinsic value and 71% in self-efficacy. The
corresponding value for test anxiety is much lower, at 14%
($`R^2 = 0.142`$). This suggests that the other four scores contain
comparatively little information about variation in test anxiety in
these generated responses, despite the negative self-regulation edge.
Predictability is calculated in sample and should not be presented as
out-of-sample predictive performance without validation.

## What does the network look like?

[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws the fitted network. With `psych_styling = TRUE`, positive edges
are green, negative edges red, edge width follows the absolute weight,
and the ring around each node shows its predictability.

``` r

cograph::splot(net, psych_styling = TRUE)
```

![](gaussian-graphical-models_files/figure-html/plot-1.png)

The figure shows the dense positive cluster of the learning constructs
and the single negative edge to test anxiety. Node positions come from a
layout algorithm, and the distance between two nodes carries no
statistical meaning.

## Do the conclusions depend on the correlation estimator?

Spearman correlations use ranks and are less affected by skewed
distributions and extreme values. Refitting with
`cor_method = "spearman"` tests whether the structure depends on the
Pearson estimator.

``` r

spearman_net <- psychnet(data = SRL_GPT, method = "glasso",
                         cor_method = "spearman")
summary(spearman_net)
#> <psychnet> glasso network
#>   nodes: 5   edges: 9   (undirected)
#>   lambda: 0.01882   gamma: 0.5
#>   optimality (KKT residual): 2.67e-10
#>   edge weight: range [-0.237, 0.424], mean 0.192
as.data.frame(spearman_net)
#>   from to      weight
#> 1  CSU IV  0.41037571
#> 2  CSU SE  0.42434529
#> 3   IV SE  0.12821527
#> 4  CSU SR  0.32164519
#> 5   IV SR  0.34602908
#> 6   SE SR  0.18400655
#> 7   IV TA  0.07280225
#> 8   SE TA  0.07509761
#> 9   SR TA -0.23654039
```

The Spearman-based model selects a larger penalty
($`\lambda = 0.01882`$) and retains 9 edges. All six positive edges
among the learning constructs remain, as does the negative edge between
self-regulation and test anxiety, although its weight changes from
-0.350 to -0.237. The smallest Pearson-based edge, `CSU`–`TA` (0.049),
is absent. These results suggest that the main pattern is similar under
the two correlation estimators, while the weaker connections involving
test anxiety are more sensitive to the choice of estimator.

## What does a post-hoc threshold change?

A post-hoc threshold removes estimated edges whose absolute weights fall
below a specified value. This may simplify a display but changes the
reported graph without refitting the model. The threshold should be
justified and, preferably, specified before examining the results; it is
not part of the EBIC selection procedure.

``` r

threshold_net <- psychnet(data = SRL_GPT, method = "glasso", threshold = 0.1)
threshold_net
#> <psychnet> glasso network
#>   nodes: 5   edges: 8   (undirected)
#>   lambda: 0.00861   gamma: 0.5
#>   optimality: unavailable for post-threshold graph  (pre-threshold residual 2.21e-10)
as.data.frame(threshold_net)
#>   from to     weight
#> 1  CSU IV  0.4116687
#> 2  CSU SE  0.3829022
#> 3   IV SE  0.1599260
#> 4  CSU SR  0.3548136
#> 5   IV SR  0.3190944
#> 6   SE SR  0.2076721
#> 7   IV TA  0.1105874
#> 8   SR TA -0.3498516
```

At a threshold of 0.1, 8 edges remain. The `CSU`–`TA` (0.049) and
`SE`–`TA` (0.099) edges are omitted, while the selected penalty remains
unchanged. Because the thresholded graph is not itself the solution of
the graphical-lasso objective, its optimality is reported as
unavailable. The printout gives the pre-threshold residual
($`2.21 \times 10^{-10}`$), and `certificate(threshold_net)` returns
`NA`. Reports should distinguish these post-hoc exclusions from edges
set to zero during model estimation.

## Do the conclusions depend on normality?

The GGM’s conditional-independence interpretation rests on multivariate
normality. The nonparanormal model (Liu et al., 2009) relaxes it: each
variable is mapped to normal scores through its ranks, and the same EBIC
graphical lasso is fitted to the transformed data. `method = "huge"`
fits this model.

``` r

npn_net <- psychnet(data = SRL_GPT, method = "huge")
summary(npn_net)
#> <psychnet> huge network
#>   nodes: 5   edges: 10   (undirected)
#>   lambda: 0.008559   gamma: 0.5
#>   optimality (KKT residual): 2.47e-10
#>   edge weight: range [-0.343, 0.426], mean 0.172
as.data.frame(npn_net)
#>    from to      weight
#> 1   CSU IV  0.42556520
#> 2   CSU SE  0.40615961
#> 3    IV SE  0.12908167
#> 4   CSU SR  0.34672549
#> 5    IV SR  0.31792725
#> 6    SE SR  0.19589994
#> 7   CSU TA  0.02390372
#> 8    IV TA  0.12511890
#> 9    SE TA  0.09634264
#> 10   SR TA -0.34267471
certificate(npn_net)
#>   method  certificate kind certified
#> 1   huge 2.471392e-10  kkt      TRUE
```

The nonparanormal fit retains all 10 edges, with weights from -0.343 to
0.426. Its certificate ($`2.47 \times 10^{-10}`$) satisfies the
numerical optimality criterion. The six positive edges among the
learning constructs remain, and the self-regulation to test anxiety edge
is still negative (-0.343). The `CSU`–`TA` edge decreases from 0.049 in
the Pearson-based GGM to 0.024. The principal pattern is therefore
similar across these fits, but this comparison alone cannot establish
that the assumptions hold or quantify the sampling stability of
individual edges.

## Reporting the analysis

A report states the variables and their measurement scale, the sample
size, the correlation estimator, the treatment of missing data, the
estimator and the value of $`\gamma`$, and any threshold. The edges are
reported as partial correlations, with centrality and predictability
defined before they are interpreted.

For these 300 complete, language-model-generated responses, an EBIC
graphical lasso ($`\gamma = 0.5`$) based on Pearson correlations
selected $`\lambda = 0.00861`$ and retained all 10 possible edges. The
largest positive partial correlation was between cognitive strategy use
and intrinsic value (0.412), and the only negative partial correlation
was between self-regulation and test anxiety (-0.350). The
Spearman-based model retained 9 edges and the nonparanormal model
retained 10; both preserved the positive learning-construct cluster and
the negative self-regulation to test anxiety association. Changes in
individual edge weights, particularly those involving test anxiety,
should be reported rather than treated as evidence of exact replication.

## Assumptions and when to use the GGM

- The variables are continuous, or ordinal with enough categories to be
  treated as continuous; otherwise use polychoric correlations
  (`cor_method = "auto"`), the Ising model for binary data, or a mixed
  graphical model when continuous and binary variables are combined.
- The joint distribution is approximately multivariate normal; the
  nonparanormal model is the sensitivity check when it is not.
- Observations are independent. Repeated measures of the same persons
  call for a model that separates within- and between-person variation,
  such as the panel graphical VAR.
- The regularized GGM suits analyses where a sparse, readable network is
  the goal. When the sample is large relative to the number of variables
  and unshrunk partial correlations are wanted, the stepwise model
  selection of `method = "ggm"` is the alternative.
- An edge represents a conditional association, not a causal effect. A
  missing edge is an estimated zero in the selected model, not proof
  that the variables are unrelated in the population. Sampling
  uncertainty and the choice of estimation method remain important.

## Mathematical details

### Optimality conditions

Write $`\mathbf{W} = \widehat{\mathbf{K}}^{-1}`$. At the minimum of the
graphical-lasso objective the diagonal satisfies $`W_{ii} = S_{ii}`$,
and for the off-diagonal entries

``` math
W_{ij} - S_{ij} = \lambda\,\operatorname{sign}(\widehat{K}_{ij})
\quad \text{if } \widehat{K}_{ij} \ne 0, \qquad
|W_{ij} - S_{ij}| \leq \lambda \quad \text{if } \widehat{K}_{ij} = 0.
```

[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports the largest violation of these conditions. For the strictly
convex graphical-lasso objective, sufficiently small violations indicate
that the numerical solution is close to the unique optimum. This
assessment concerns optimisation rather than statistical uncertainty.

### Selection of the penalty

The path holds candidate values of $`\lambda`$ on a logarithmic grid.
Each is fitted, the EBIC is computed, and the value with the smallest
EBIC is refitted at a tighter numerical tolerance before it is returned.

### Predictability

The residual variance of node $`j`$ given all other nodes is
$`1/K_{jj}`$. With marginal variance $`S_{jj}`$, the explained
proportion is

``` math
R_j^2 = 1 - \frac{1}{K_{jj}S_{jj}}.
```

## How to cite this article

Saqr, M., & López-Pernas, S. (2026). *Estimating and interpreting
Gaussian graphical models with the psychnets R package*. psychnets.
<https://doi.org/10.5281/zenodo.24569494> ·
[HTML](https://pak.dynasite.org/psychnets/articles/gaussian-graphical-models.html)
·
[PDF](https://pak.dynasite.org/psychnets/pdf/gaussian-graphical-models.pdf)

## References

Borsboom, D., Deserno, M. K., Rhemtulla, M., Epskamp, S., Fried, E. I.,
McNally, R. J., Robinaugh, D. J., Perugini, M., Dalege, J., Costantini,
G., Isvoranu, A.-M., Wysocki, A. C., van Borkulo, C. D., van Bork, R., &
Waldorp, L. J. (2021). Network analysis of multivariate data in
psychological science. *Nature Reviews Methods Primers*, 1, 58.
[doi:10.1038/s43586-021-00055-w](https://doi.org/10.1038/s43586-021-00055-w)

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods*, 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

Epskamp, S., & Fried, E. I. (2018). A tutorial on regularized partial
correlation networks. *Psychological Methods*, 23(4), 617-634.
[doi:10.1037/met0000167](https://doi.org/10.1037/met0000167)

Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems*, 23, 604-612.

Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
covariance estimation with the graphical lasso. *Biostatistics*, 9(3),
432-441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853-861.
[doi:10.3758/s13428-017-0910-x](https://doi.org/10.3758/s13428-017-0910-x)

Isvoranu, A.-M., & Epskamp, S. (2023). Which estimation method to choose
in network psychometrics? Deriving guidelines for applied researchers.
*Psychological Methods*, 28(4), 925-946.
[doi:10.1037/met0000439](https://doi.org/10.1037/met0000439)

Liu, H., Lafferty, J., & Wasserman, L. (2009). The nonparanormal:
Semiparametric estimation of high dimensional undirected graphs.
*Journal of Machine Learning Research*, 10, 2295-2328.

Malmberg, J., Saqr, M., Järvenoja, H., Haataja, E., Pijeira-Díaz, H. J.,
& Järvelä, S. (2022). Modeling the complex interplay between monitoring
events for regulated learning with psychological networks. In M.
Giannakos, D. Spikol, D. Di Mitri, K. Sharma, X. Ochoa, & R. Hammad
(Eds.), *The multimodal learning analytics handbook* (pp. 79-104).
Springer.
[doi:10.1007/978-3-031-08076-0_4](https://doi.org/10.1007/978-3-031-08076-0_4)

Malmberg, J., Saqr, M., Järvenoja, H., & Järvelä, S. (2022). How the
monitoring events of individual students are associated with phases of
regulation: A network analysis approach. *Journal of Learning
Analytics*, 9(1), 77-92.
[doi:10.18608/jla.2022.7429](https://doi.org/10.18608/jla.2022.7429)

Robinaugh, D. J., Millner, A. J., & McNally, R. J. (2016). Identifying
highly influential nodes in the complicated grief network. *Journal of
Abnormal Psychology*, 125(6), 747-757.
[doi:10.1037/abn0000181](https://doi.org/10.1037/abn0000181)

Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological networks: A
modern approach to analysis of learning and complex learning processes.
In M. Saqr & S. López-Pernas (Eds.), *Learning analytics methods and
tutorials: A practical guide using R* (pp. 639-671). Springer.
[doi:10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)

Saqr, M., Peeters, W., & Viberg, O. (2021). The relational, co-temporal,
contemporaneous, and longitudinal dynamics of self-regulation for
academic writing. *Research and Practice in Technology Enhanced
Learning*, 16, 29.
[doi:10.1186/s41039-021-00175-7](https://doi.org/10.1186/s41039-021-00175-7)

Vogelsmeier, L. V. D. E., Oliveira, E., Misiejuk, K., López-Pernas, S.,
& Saqr, M. (2025). Delving into the psychology of machines: Exploring
the structure of self-regulated learning via LLM-generated survey
responses. *Computers in Human Behavior*, 173, 108769.
[doi:10.1016/j.chb.2025.108769](https://doi.org/10.1016/j.chb.2025.108769)

Williams, D. R., & Rast, P. (2020). Back to the basics: Rethinking
partial correlation network methodology. *British Journal of
Mathematical and Statistical Psychology*, 73(2), 187-212.
[doi:10.1111/bmsp.12173](https://doi.org/10.1111/bmsp.12173)
