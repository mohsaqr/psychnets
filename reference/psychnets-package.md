# psychnets: Psychometric Network Estimation

`psychnets` estimates psychometric networks and assesses their accuracy,
stability and differences. A psychometric network represents
psychological variables, such as items, symptoms or behaviours, as nodes
and their associations as edges; it shows which variables are directly
related, which are central, and how the structure differs between
groups.

The cross-sectional estimators are reached through
[`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
or called directly: Gaussian graphical models
([`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md),
[`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md),
[`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)),
information-filtering networks
([`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md),
[`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)),
relative-importance networks
([`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md)),
correlation and partial-correlation networks
([`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md),
[`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)),
Ising models for binary data
([`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md),
[`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md))
and mixed graphical models
([`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)).
Panel data from a few waves are modelled with the panel graphical VAR
([`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md)).
Fitted networks are read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md),
[`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
and
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md);
their accuracy, stability and group differences are assessed with
[`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md),
[`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
and
[`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md).

## Certificates

Each regularized or constrained estimator stores how far its network is
from satisfying the optimality conditions of its objective, so that the
fitted network can be checked without a second solver.
[`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
reports this value for any fitted network.

## See also

Useful links:

- <https://pak.dynasite.org/psychnets>

- <https://github.com/mohsaqr/psychnets>

- Report bugs at <https://github.com/mohsaqr/psychnets/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

Authors:

- Mohammed Saqr <saqr@saqr.me>
  ([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]

- Sonsoles López-Pernas ([ORCID](https://orcid.org/0000-0002-9621-1392))
