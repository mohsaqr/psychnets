# net_compare() and its accessors: the network comparison test of two groups (van Borkulo et al. 2022)
# against NetworkComparisonTest::NCT(), run live on the same data after the same set.seed(). Both draw the
# permutations with the same calls to sample(), so the permuted splits are identical and the p-values are
# comparable one to one (the seed is the shared supply of draws). Four configurations of simulated groups
# (n = 120 each, five variables): independent and paired designs, abs = FALSE; and weighted = FALSE, where NCT's
# glasso precision (thr = 1e-4) leaves a residue that counts as an edge in one permuted split (a known
# difference, with the split's networks compared).

configurations <- data.frame(seed = c(1, 2, 3), paired = c(FALSE, TRUE, FALSE), abs = c(TRUE, TRUE, FALSE),
                             weighted = c(TRUE, TRUE, TRUE))
# 100 permutations (raised 2026-11-23 from 15): a p-value is resolved to 1/101, so the parity of the p-values
# rests on a meaningful number of permuted splits; it costs about 25 seconds per configuration.
permutations <- 100L

#' NCT's edge table keyed as psychnets names an edge: `from` the node earlier in the matrix
nct_edges <- function(reference, nodes) {
  table <- reference$einv.pvals
  first <- match(as.character(table$Var1), nodes) < match(as.character(table$Var2), nodes)
  data.frame(from = ifelse(first, as.character(table$Var1), as.character(table$Var2)),
             to = ifelse(first, as.character(table$Var2), as.character(table$Var1)),
             p_value = table[["p-value"]])
}

invisible(lapply(seq_len(nrow(configurations)), \(k) {
  configuration <- configurations[k, ]
  label <- sprintf("configuration %d (paired = %s, abs = %s, weighted = %s)", k, configuration$paired,
                   configuration$abs, configuration$weighted)
  # With abs = FALSE the global statistic is the sum of signed weights, which psychnets names global expected
  # influence; NCT keeps the name glstrinv (global strength).
  global <- if (configuration$abs) "global strength" else "global expected influence"
  test_that(sprintf("the network comparison test equals NetworkComparisonTest::NCT(), %s", label), {
    skip_without_reference("NetworkComparisonTest")
    groups <- ref_nct_groups(configuration$seed)
    nodes <- colnames(groups$a)
    set.seed(configuration$seed)
    ours <- net_compare(groups$a, groups$b, iter = permutations, gamma = 0.5, paired = configuration$paired,
                        abs = configuration$abs, weighted = configuration$weighted)
    set.seed(configuration$seed)
    reference <- NetworkComparisonTest::NCT(groups$a, groups$b, gamma = 0.5, it = permutations,
                                            paired = configuration$paired, abs = configuration$abs,
                                            weighted = configuration$weighted, test.edges = TRUE, edges = "all",
                                            progressbar = FALSE, verbose = FALSE)
    # The same permuted rows: a p-value is a count of permutations over 101 and must be identical.
    expect_reference(summary(ours), data.frame(test = c("network structure", global),
                                               p_value = c(reference$nwinv.pval, reference$glstrinv.pval)),
                     reference = "NetworkComparisonTest", covers = c("net_compare()", "summary.psychnet_nct()"),
                     quantity = sprintf("p-value, %s", label), tolerance = 1e-12, key = "test")
    expect_reference(as.data.frame(ours), nct_edges(reference, nodes), reference = "NetworkComparisonTest",
                     covers = c("net_compare()", "as.data.frame.psychnet_nct()"),
                     quantity = sprintf("edge p-value, %s", label), tolerance = 1e-12, key = c("from", "to"))
    # NCT estimates with qgraph's EBICglasso, whose glasso stops at thr = 1e-4: weights, strengths and the
    # observed and permuted statistics agree to 1e-3 (the package's own test).
    weights_1 <- reference$nw1
    weights_2 <- reference$nw2
    dimnames(weights_1) <- dimnames(weights_2) <- list(nodes, nodes)
    expect_reference(as.data.frame(ours), ref_edge_table(weight_1 = weights_1, weight_2 = weights_2),
                     reference = "NetworkComparisonTest", covers = c("net_compare()", "as.data.frame.psychnet_nct()"),
                     quantity = sprintf("edge weight per group, %s", label), tolerance = 1e-3, key = c("from", "to"))
    expect_reference(summary(ours), data.frame(test = c("network structure", global),
                                               observed = c(reference$nwinv.real, reference$glstrinv.real)),
                     reference = "NetworkComparisonTest", covers = "summary.psychnet_nct()",
                     quantity = sprintf("observed statistic, %s", label), tolerance = 1e-3, key = "test")
    expect_reference(ours$strength$by_group, reference$glstrinv.sep, reference = "NetworkComparisonTest",
                     covers = c("net_compare()", "$.psychnet_nct()"), quantity = sprintf("global strength per group, %s", label),
                     tolerance = 1e-3, note = "Reached with `$strength$by_group`: no accessor (API gap).")
    expect_reference(ours$structure$perm, reference$nwinv.perm, reference = "NetworkComparisonTest",
                     covers = c("net_compare()", "$.psychnet_nct()"),
                     quantity = sprintf("permuted structure statistic, %s", label), tolerance = 1e-3)
    expect_reference(ours$strength$perm, reference$glstrinv.perm, reference = "NetworkComparisonTest",
                     covers = c("net_compare()", "$.psychnet_nct()"),
                     quantity = sprintf("permuted global strength statistic, %s", label), tolerance = 1e-3)
  })
}))

test_that("weighted = FALSE: NCT's observed statistics agree; one permuted split differs by a glasso residue", {
  skip_without_reference("NetworkComparisonTest")
  groups <- ref_nct_groups(4)
  nodes <- colnames(groups$a)
  label <- "configuration 4 (paired = FALSE, abs = TRUE, weighted = FALSE)"
  set.seed(4)
  ours <- net_compare(groups$a, groups$b, iter = permutations, gamma = 0.5, weighted = FALSE)
  set.seed(4)
  reference <- NetworkComparisonTest::NCT(groups$a, groups$b, gamma = 0.5, it = permutations, weighted = FALSE,
                                          test.edges = TRUE, edges = "all", progressbar = FALSE, verbose = FALSE)
  weights_1 <- reference$nw1
  weights_2 <- reference$nw2
  dimnames(weights_1) <- dimnames(weights_2) <- list(nodes, nodes)
  # glasso's thr = 1e-4 in NCT, as above: 1e-3.
  expect_reference(as.data.frame(ours), ref_edge_table(weight_1 = weights_1, weight_2 = weights_2),
                   reference = "NetworkComparisonTest", covers = c("net_compare()", "as.data.frame.psychnet_nct()"),
                   quantity = sprintf("edge weight per group, %s", label), tolerance = 1e-3, key = c("from", "to"))
  expect_reference(summary(ours), data.frame(test = c("network structure", "global strength"),
                                             observed = c(reference$nwinv.real, reference$glstrinv.real)),
                   reference = "NetworkComparisonTest", covers = "summary.psychnet_nct()",
                   quantity = sprintf("observed statistic, %s", label), tolerance = 1e-3, key = "test")
  # The split that differs, drawn as NCT draws it (one sample() of the first group's rows per permutation after
  # set.seed(4); permutation 67), estimated by NCT's own estimator and by psychnets: the networks agree to
  # glasso's precision, but one edge that is exactly 0 in psychnets is left at a residue of order 1e-7 by NCT.
  pooled <- rbind(groups$a, groups$b)
  set.seed(4)
  split <- lapply(seq_len(67L), \(i) sample(seq_len(nrow(pooled)), nrow(groups$a), replace = FALSE))[[67L]]
  nct_network <- NetworkComparisonTest:::NCT_estimator_GGM(pooled[split, ], gamma = 0.5)
  dimnames(nct_network) <- list(nodes, nodes)
  expect_reference(as.data.frame(ebic_glasso(pooled[split, ], gamma = 0.5), include_zero = TRUE),
                   ref_edges(nct_network), reference = "NetworkComparisonTest", covers = "net_compare()",
                   quantity = "edge weight, first group of permuted split 67 (seed 4)", tolerance = 1e-3,
                   key = c("from", "to"),
                   note = "NCT's estimator is NetworkComparisonTest:::NCT_estimator_GGM (qgraph::EBICglasso).")
  residue <- paste("NCT's EBICglasso stops glasso at thr = 1e-4. In permuted split 67 of 100 (seed 4) it leaves an",
                   "edge at a residue of order 1e-7 (-1.58e-7 measured) that is exactly 0 in psychnets (compared",
                   "above); with weighted = FALSE that residue is an edge, so NCT's permuted statistics are 0 where",
                   "psychnets' are 1, and the p-values that count that split move by 1/101.")
  expect_known_difference(ours$structure$perm, reference$nwinv.perm, reference = "NetworkComparisonTest",
                          covers = c("net_compare()", "$.psychnet_nct()"),
                          quantity = sprintf("permuted structure statistic, %s", label), tolerance = 1e-3,
                          reason = residue)
  expect_known_difference(ours$strength$perm, reference$glstrinv.perm, reference = "NetworkComparisonTest",
                          covers = c("net_compare()", "$.psychnet_nct()"),
                          quantity = sprintf("permuted global strength statistic, %s", label), tolerance = 1e-3,
                          reason = residue)
  expect_known_difference(summary(ours), data.frame(test = c("network structure", "global strength"),
                                                    p_value = c(reference$nwinv.pval, reference$glstrinv.pval)),
                          reference = "NetworkComparisonTest", covers = c("net_compare()", "summary.psychnet_nct()"),
                          quantity = sprintf("p-value, %s", label), tolerance = 1e-12, key = "test", reason = residue)
  expect_known_difference(as.data.frame(ours), nct_edges(reference, nodes), reference = "NetworkComparisonTest",
                          covers = c("net_compare()", "as.data.frame.psychnet_nct()"),
                          quantity = sprintf("edge p-value, %s", label), tolerance = 1e-12, key = c("from", "to"),
                          reason = residue)
})

test_that("Holm-adjusted edge p-values equal NCT's p.adjust.methods = \"holm\"", {
  skip_without_reference("NetworkComparisonTest")
  groups <- ref_nct_groups(1)
  set.seed(1)
  ours <- net_compare(groups$a, groups$b, iter = permutations, gamma = 0.5, p_adjust = "holm")
  set.seed(1)
  reference <- NetworkComparisonTest::NCT(groups$a, groups$b, gamma = 0.5, it = permutations, test.edges = TRUE,
                                          edges = "all", p.adjust.methods = "holm", progressbar = FALSE, verbose = FALSE)
  expect_reference(as.data.frame(ours), nct_edges(reference, colnames(groups$a)), reference = "NetworkComparisonTest",
                   covers = "net_compare()", quantity = "Holm-adjusted edge p-value", tolerance = 1e-12,
                   key = c("from", "to"))
})
