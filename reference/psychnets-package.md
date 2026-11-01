# psychnets: Clean-Room Base-R Psychometric Network Estimation

Estimates psychometric network models in pure base R, with no compiled
dependencies and dependency-free correctness certificates. See
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
for the unified entry point to the cross-sectional estimators
([`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md),
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md),
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md),
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md),
[`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md),
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md),
[`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md),
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md),
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md))
and
[`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md)
for the panel graphical VAR. Fitted networks are read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md),
and their accuracy and stability are assessed with
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
and
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md).

## Certification

Regularized estimators are graded against their own objective rather
than an external solver. For the Gaussian graphical model,
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
gives the stationarity (KKT) residual: zero certifies the unique global
optimum of the convex objective. Every fitted network stores its
residual, and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports it as a one-row table for any estimator.

## See also

Useful links:

- <https://pak.dynasite.org/psychnets>

- <https://github.com/mohsaqr/psychnets>

- Report bugs at <https://github.com/mohsaqr/psychnets/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me> \[copyright holder\]

Authors:

- Sonsoles López-Pernas
