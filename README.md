# psychnet

Clean-room, base-R psychometric network estimation. The R counterpart to the
[`psychaj`](https://github.com/mohsaqr/psychaj) TypeScript library.

`psychnet` estimates the cross-sectional network models used in psychometrics —
correlation, partial correlation, EBIC-regularized Gaussian graphical models
(graphical lasso), the Ising model for binary data, and mixed graphical
models — **reimplemented from first principles in base R, with no compiled
dependencies**.

## Why "clean room"

Each regularized estimator ships a **dependency-free correctness certificate**.
Rather than trusting a result because it matches an external solver, `psychnet`
grades it against its own convex objective:

- For the Gaussian graphical model, `glasso_kkt()` returns the **stationarity
  (KKT) residual** of the fitted precision matrix. Because the graphical-lasso
  objective is strictly convex, its minimiser is unique, so a near-zero residual
  certifies the global optimum — with no reference solver involved.
- For the nodewise lasso behind `ising_fit()` and `mgm_fit()`, `glm_lasso_kkt()`
  plays the same role for the penalized-likelihood optimum.

Every fitted network carries its certificate in `$kkt`. External packages
(`qgraph`, `IsingFit`, `mgm`, `bootnet`) are used only as *cross-checks at
independent-solver precision*, never as the definition of correct. In practice
`psychnet`'s graphical lasso is **provably no further from the optimum** than
`qgraph::EBICglasso()`, which stops at glasso's default `thr = 1e-4`.

## Estimators (v0.1)

Full cross-sectional `bootnet` / `psychaj` estimator parity (temporal models
excluded), each pure base R and self-certified.

| Function | Model | Data |
|---|---|---|
| `cor_network()` | marginal correlations (+ significance threshold) | continuous |
| `pcor_network()` | partial correlations (+ significance threshold) | continuous |
| `ebic_glasso()` | EBIC graphical lasso (GGM) | continuous |
| `huge_network()` | nonparanormal graphical model | continuous |
| `ggmncv_network()` | non-convex GGM (SCAD / MCP / atan) | continuous |
| `ggm_modselect()` | unregularized stepwise GGM | continuous |
| `tmfg_network()` | Triangulated Maximally Filtered Graph | continuous |
| `logo_network()` | Local-Global sparse inverse covariance | continuous |
| `relimp_network()` | relative importance (LMG / Shapley, directed) | continuous |
| `ising_fit()` | Ising model (L1-penalized) | binary |
| `ising_sampler()` | Ising model (unregularized + Wald pruning) | binary |
| `mgm_fit()` | mixed graphical model | gaussian + binary |
| `estimate_network()` | unified front door (à la `bootnet`) | — |
| `centrality()` | strength + expected influence | — |
| `predictability()` | per-node R² / classification accuracy | — |
| `bootstrap_network()` | edge / centrality accuracy CIs | — |
| `centrality_stability()` | CS-coefficient (case-dropping) | — |
| `nct()` | network comparison test | continuous |

## Example

```r
library(psychnet)

# Continuous data -> EBIC graphical lasso, self-certified
S   <- 0.4^abs(outer(1:8, 1:8, "-"))      # an AR(1) correlation matrix
fit <- ebic_glasso(cor_matrix = S, n = 250)
fit                                        # prints nodes, edges, lambda, KKT residual

as.data.frame(fit)                         # tidy edge list (from, to, weight)
centrality(fit)                            # tidy per-node strength / expected influence
glasso_kkt(fit$precision, S, fit$lambda)   # the certificate, directly
```

## Design

- **Base R only.** `Imports: stats`. No `glasso`, `glmnet`, or `qgraph`.
- **Tidy surface.** Verbs take simple named arguments and return a `psychnet`
  object with `print`/`summary`/`as.data.frame` methods; `centrality()` returns
  a one-row-per-node data frame.
- **Self-verifying.** Correctness is certified by the mathematics, not by an
  external dependency.

Roadmap: external cross-check suite at solver precision; multinormal / Poisson
`mgm`; the `psychaj` framework-tier extras (bridge / betweenness / closeness
centrality, community detection, graph metrics); plot methods.
