# psychnets: Clean-Room Base-R Psychometric Network Estimation

Estimates cross-sectional psychometric network models in pure base R,
with no compiled dependencies and dependency-free correctness
certificates. See
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
for the unified entry point and
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md),
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md),
[`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
and
[`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
for the individual estimators.

## Certification

Regularized estimators are graded against their own convex objective
rather than an external solver. For the Gaussian graphical model,
[`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
returns the stationarity (KKT) residual: zero certifies the unique
global optimum. Every
[`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
result carries this value in `$kkt`.

## See also

Useful links:

- <https://github.com/mohsaqr/psychnet>

- Report bugs at <https://github.com/mohsaqr/psychnet/issues>

## Author

**Maintainer**: Mohammed Saqr <ueflaunit@gmail.com>
