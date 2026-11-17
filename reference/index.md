# Package index

## Overview

- [`psychnets`](https://pak.dynasite.org/psychnets/reference/psychnets-package.md)
  [`psychnets-package`](https://pak.dynasite.org/psychnets/reference/psychnets-package.md)
  : psychnets: Psychometric Network Estimation
- [`psychnet()`](https://pak.dynasite.org/psychnets/reference/psychnet.md)
  : Estimate a psychometric network with any of the package's estimators

## Gaussian graphical models

Partial-correlation networks for continuous or ordinal indicators,
regularized, selected stepwise, or nonparanormal.

- [`ebic_glasso()`](https://pak.dynasite.org/psychnets/reference/ebic_glasso.md)
  : Estimate a Gaussian graphical model with the EBIC graphical lasso
- [`ggm_modselect()`](https://pak.dynasite.org/psychnets/reference/ggm_modselect.md)
  : Select an unregularized Gaussian graphical model by stepwise search
- [`pcor_network()`](https://pak.dynasite.org/psychnets/reference/pcor_network.md)
  : Estimate a partial correlation network
- [`cor_network()`](https://pak.dynasite.org/psychnets/reference/cor_network.md)
  : Estimate a correlation network
- [`huge_network()`](https://pak.dynasite.org/psychnets/reference/huge_network.md)
  : Estimate a nonparanormal Gaussian graphical model

## Information-filtering and relative-importance networks

Sparse networks built from the correlation structure (TMFG, LoGo) and
networks of shared explained variance.

- [`tmfg_network()`](https://pak.dynasite.org/psychnets/reference/tmfg_network.md)
  : Filter a correlation network with the triangulated maximally
  filtered graph
- [`logo_network()`](https://pak.dynasite.org/psychnets/reference/logo_network.md)
  : Estimate a sparse Gaussian graphical model on a filtered graph
  (LoGo)
- [`relimp_network()`](https://pak.dynasite.org/psychnets/reference/relimp_network.md)
  : Estimate a directed relative-importance network

## Binary and mixed data

Ising networks for binary indicators and mixed graphical models for
variables of different types, including moderation.

- [`ising_fit()`](https://pak.dynasite.org/psychnets/reference/ising_fit.md)
  : Estimate a regularized Ising network for binary items
- [`ising_sampler()`](https://pak.dynasite.org/psychnets/reference/ising_sampler.md)
  : Estimate an unregularized Ising network for binary items
- [`mgm_fit()`](https://pak.dynasite.org/psychnets/reference/mgm_fit.md)
  : Estimate a mixed graphical model for continuous and binary variables
- [`condition()`](https://pak.dynasite.org/psychnets/reference/condition.md)
  [`as.data.frame(`*`<psychnet_moderated>`*`)`](https://pak.dynasite.org/psychnets/reference/condition.md)
  [`print(`*`<psychnet_moderated>`*`)`](https://pak.dynasite.org/psychnets/reference/condition.md)
  : Read a moderated mixed graphical model at a value of the moderator
- [`dichotomize()`](https://pak.dynasite.org/psychnets/reference/dichotomize.md)
  : Split numeric items into binary 0/1 variables

## Panel and within-person networks

Temporal, contemporaneous and between-person networks from repeated
measurements.

- [`panel_gvar()`](https://pak.dynasite.org/psychnets/reference/panel_gvar.md)
  : Estimate temporal, contemporaneous and between-person networks from
  panel data
- [`as.data.frame(`*`<psychnet_panelgvar>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)
  [`print(`*`<psychnet_panelgvar>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)
  [`summary(`*`<psychnet_panelgvar>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)
  [`plot(`*`<psychnet_panelgvar>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_panelgvar.md)
  : Read the networks, parameters and fit of a panel graphical VAR
- [`print(`*`<psychnet_multilevel>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
  [`as.data.frame(`*`<psychnet_multilevel>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
  [`summary(`*`<psychnet_multilevel>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
  [`plot(`*`<psychnet_multilevel>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_multilevel.md)
  : Within-actor and between-actor networks from event data

## Groups and event data

One network per group, and networks from logs of coded events.

- [`print(`*`<psychnet_group>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  [`as.data.frame(`*`<psychnet_group>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  [`summary(`*`<psychnet_group>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_group.md)
  : Networks estimated separately for each group
- [`print(`*`<psychnet_result_group>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md)
  [`as.data.frame(`*`<psychnet_result_group>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet_result_group.md)
  : Results of an analysis applied to networks per group
- [`event_frequencies()`](https://pak.dynasite.org/psychnets/reference/event_frequencies.md)
  : Count the actions of each actor and session in an event log

## Preparing data

- [`cor_auto()`](https://pak.dynasite.org/psychnets/reference/cor_auto.md)
  : Compute correlations suited to ordinal and continuous variables
- [`redundancy()`](https://pak.dynasite.org/psychnets/reference/redundancy.md)
  [`print(`*`<psychnet_redundancy>`*`)`](https://pak.dynasite.org/psychnets/reference/redundancy.md)
  : Detect redundant pairs of items
- [`net_aggregate()`](https://pak.dynasite.org/psychnets/reference/net_aggregate.md)
  : Aggregate communities of items into a network between communities

## Node and network measures

- [`net_centralities()`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
  [`plot(`*`<psychnet_centrality>`*`)`](https://pak.dynasite.org/psychnets/reference/net_centralities.md)
  : Compute the centrality of each node in a network
- [`net_bridge()`](https://pak.dynasite.org/psychnets/reference/net_bridge.md)
  [`plot(`*`<psychnet_bridge>`*`)`](https://pak.dynasite.org/psychnets/reference/net_bridge.md)
  : Compute bridge centrality between communities of nodes
- [`net_edge_betweenness()`](https://pak.dynasite.org/psychnets/reference/net_edge_betweenness.md)
  : Compute the edge betweenness of a network
- [`net_clustering()`](https://pak.dynasite.org/psychnets/reference/net_clustering.md)
  : Compute weighted clustering coefficients of the nodes
- [`net_smallworld()`](https://pak.dynasite.org/psychnets/reference/net_smallworld.md)
  : Compute the small-world index of a network
- [`net_predict()`](https://pak.dynasite.org/psychnets/reference/net_predict.md)
  : Compute how well each node is predicted by the network
- [`node_predictability()`](https://pak.dynasite.org/psychnets/reference/node_predictability.md)
  : Return node predictability as a vector for plotting

## Accuracy and stability

Bootstrap intervals, difference tests, and the stability and reliability
of edges and centralities under case dropping and splitting.

- [`net_boot()`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  [`as.data.frame(`*`<psychnet_bootstrap>`*`)`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  [`print(`*`<psychnet_bootstrap>`*`)`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  [`plot(`*`<psychnet_bootstrap>`*`)`](https://pak.dynasite.org/psychnets/reference/net_boot.md)
  : Bootstrap the edges and centralities of a network
- [`difference_test()`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
  [`plot(`*`<psychnet_difference>`*`)`](https://pak.dynasite.org/psychnets/reference/difference_test.md)
  : Test differences between edges or between centralities within a
  network
- [`plot(`*`<psychnet_stability>`*`)`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  [`net_stability()`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  [`print(`*`<psychnet_stability>`*`)`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  [`as.data.frame(`*`<psychnet_stability>`*`)`](https://pak.dynasite.org/psychnets/reference/net_stability.md)
  : Assess the stability of centralities when cases are dropped
- [`net_casedrop_reliability()`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  [`print(`*`<psychnet_casedrop>`*`)`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  [`plot(`*`<psychnet_casedrop>`*`)`](https://pak.dynasite.org/psychnets/reference/net_casedrop_reliability.md)
  : Assess the stability of edge weights when cases are dropped
- [`net_split_reliability()`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
  [`print(`*`<psychnet_reliability>`*`)`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
  [`plot(`*`<psychnet_reliability>`*`)`](https://pak.dynasite.org/psychnets/reference/net_split_reliability.md)
  : Assess the split-half reliability of a network's edge weights

## Comparing networks

- [`net_compare()`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  [`summary(`*`<psychnet_nct>`*`)`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  [`as.data.frame(`*`<psychnet_nct>`*`)`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  [`print(`*`<psychnet_nct>`*`)`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  [`plot(`*`<psychnet_nct>`*`)`](https://pak.dynasite.org/psychnets/reference/net_compare.md)
  : Compare the networks of two groups with a permutation test

## Certificates

Optimality and structural checks that a fitted network is the solution
of its estimator.

- [`certificate()`](https://pak.dynasite.org/psychnets/reference/certificate.md)
  : Check that a fitted network solves its estimator
- [`glasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glasso_kkt.md)
  : Check that a precision matrix solves the graphical lasso
- [`glm_lasso_kkt()`](https://pak.dynasite.org/psychnets/reference/glm_lasso_kkt.md)
  : Check that coefficients solve a lasso-penalized regression
- [`ggm_support_kkt()`](https://pak.dynasite.org/psychnets/reference/ggm_support_kkt.md)
  : Check that a precision matrix is the maximum-likelihood network on a
  graph
- [`lmg_certificate()`](https://pak.dynasite.org/psychnets/reference/lmg_certificate.md)
  : Check the decomposition of a relative-importance network
- [`tmfg_certificate()`](https://pak.dynasite.org/psychnets/reference/tmfg_certificate.md)
  : Check that a network has the structure of a TMFG

## Result objects

- [`plot(`*`<psychnet>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
  [`` `$`( ``*`<psychnet>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
  [`print(`*`<psychnet>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
  [`as.data.frame(`*`<psychnet>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
  [`summary(`*`<psychnet>`*`)`](https://pak.dynasite.org/psychnets/reference/psychnet-class.md)
  : Fitted psychometric network

## Crosswalk

Arguments of qgraph, bootnet, IsingFit and mgm mapped to psychnets.

- [`net_crosswalk()`](https://pak.dynasite.org/psychnets/reference/net_crosswalk.md)
  : Map the arguments of established network packages to psychnets

## Data

- [`SRL_GPT`](https://pak.dynasite.org/psychnets/reference/SRL.md)
  [`SRL_Gemini`](https://pak.dynasite.org/psychnets/reference/SRL.md)
  [`SRL_Claude`](https://pak.dynasite.org/psychnets/reference/SRL.md)
  [`SRL_Mistral`](https://pak.dynasite.org/psychnets/reference/SRL.md)
  [`SRL_LLaMa`](https://pak.dynasite.org/psychnets/reference/SRL.md) :
  Self-regulated learning scores generated by large language models
- [`srl_panel`](https://pak.dynasite.org/psychnets/reference/srl_panel.md)
  : Simulated weekly panel of self-regulated learning
