# net_centralities(), net_clustering(), net_edge_betweenness(), net_bridge() and net_smallworld() on a
# ten-node signed weighted network (the EBIC graphical lasso network of the bfi items, data/network.tsv).
# References: qgraph::centrality_auto() and clustcoef_auto(), igraph, networktools::bridge() and
# expectedInf(). Shortest paths run over 1 / |weight|.

network <- ref_read_matrix(reference_data("data/network.tsv"))
p <- nrow(network)
nodes <- rownames(network)
graph_of <- \(m) igraph::graph_from_adjacency_matrix(abs(m), mode = "undirected", weighted = TRUE, diag = FALSE)

test_that("strength and expected influence equal qgraph::centrality_auto()", {
  skip_without_reference("qgraph")
  central <- qgraph::centrality_auto(network)$node.centrality
  expected <- data.frame(node = rownames(central), strength = central$Strength,
                         expected_influence = central$ExpectedInfluence)
  # Sums of the same weights.
  expect_reference(net_centralities(network, measures = c("strength", "expected_influence")), expected,
                   reference = "qgraph", covers = "net_centralities()", quantity = "centrality", tolerance = 1e-12,
                   key = "node")
})

test_that("closeness and betweenness: qgraph's convention differs; after the stated rescaling they agree", {
  skip_without_reference("qgraph")
  central <- qgraph::centrality_auto(network)$node.centrality
  ours <- net_centralities(network, measures = c("closeness", "betweenness"))
  raw <- data.frame(node = rownames(central), closeness = central$Closeness, betweenness = central$Betweenness)
  expect_known_difference(ours, raw, reference = "qgraph", covers = "net_centralities()",
                          quantity = "closeness and betweenness, qgraph's scale", tolerance = 1e-10, key = "node",
                          reason = paste("net_centralities() normalizes as igraph does (normalized = TRUE): closeness is",
                                         "(p - 1) / sum of distances and betweenness is divided by (p - 1)(p - 2) / 2;",
                                         "qgraph reports 1 / sum of distances and a count of shortest paths."),
                          agrees_with = "qgraph closeness x (p - 1), betweenness / ((p - 1)(p - 2) / 2)")
  rescaled <- data.frame(node = rownames(central), closeness = central$Closeness * (p - 1),
                         betweenness = central$Betweenness / ((p - 1) * (p - 2) / 2))
  expect_reference(ours, rescaled, reference = "qgraph", covers = "net_centralities()",
                   quantity = "closeness and betweenness, rescaled as documented", tolerance = 1e-10, key = "node")
})

test_that("strength, closeness and betweenness equal igraph", {
  skip_without_reference("igraph")
  graph <- graph_of(network)
  lengths <- 1 / igraph::E(graph)$weight
  expected <- data.frame(node = igraph::V(graph)$name, strength = igraph::strength(graph),
                         closeness = igraph::closeness(graph, weights = lengths, normalized = TRUE),
                         betweenness = igraph::betweenness(graph, weights = lengths, normalized = TRUE))
  expect_reference(net_centralities(network, measures = c("strength", "closeness", "betweenness")), expected,
                   reference = "igraph", covers = "net_centralities()", quantity = "centrality", tolerance = 1e-10,
                   key = "node")
})

test_that("two-step expected influence equals networktools::expectedInf()", {
  skip_without_reference("networktools")
  reference <- networktools::expectedInf(network, step = 2)
  expect_reference(net_centralities(network, measures = "expected_influence_2step"),
                   data.frame(node = nodes, expected_influence_2step = as.numeric(reference$step2)),
                   reference = "networktools", covers = "net_centralities()", quantity = "two-step expected influence",
                   tolerance = 1e-10, key = "node")
})

test_that("local clustering coefficients equal qgraph::clustcoef_auto()", {
  skip_without_reference("qgraph")
  coefficients <- qgraph::clustcoef_auto(network)
  # Closed-form ratios of sums of products of the same weights.
  expect_reference(net_clustering(network), data.frame(node = rownames(coefficients), coefficients),
                   reference = "qgraph", covers = "net_clustering()", quantity = "clustering coefficient",
                   tolerance = 1e-10, key = "node")
})

test_that("edge betweenness equals igraph::edge_betweenness()", {
  skip_without_reference("igraph")
  graph <- graph_of(network)
  ends <- igraph::ends(graph, igraph::E(graph))
  lengths <- 1 / igraph::E(graph)$weight
  # psychnets names an edge by its upper-triangle orientation: the node earlier in the matrix is `from`.
  first <- match(ends[, 1L], nodes) < match(ends[, 2L], nodes)
  expected <- data.frame(from = ifelse(first, ends[, 1L], ends[, 2L]), to = ifelse(first, ends[, 2L], ends[, 1L]),
                         edge_betweenness = igraph::edge_betweenness(graph, weights = lengths, directed = FALSE))
  expect_reference(net_edge_betweenness(network), expected, reference = "igraph", covers = "net_edge_betweenness()",
                   quantity = "edge betweenness, lengths 1 / |weight|", tolerance = 1e-10, key = c("from", "to"))
  unweighted <- (abs(network) > 0) * 1
  graph_u <- igraph::graph_from_adjacency_matrix(unweighted, mode = "undirected", diag = FALSE)
  ends_u <- igraph::ends(graph_u, igraph::E(graph_u))
  first_u <- match(ends_u[, 1L], nodes) < match(ends_u[, 2L], nodes)
  expect_reference(net_edge_betweenness(unweighted),
                   data.frame(from = ifelse(first_u, ends_u[, 1L], ends_u[, 2L]),
                              to = ifelse(first_u, ends_u[, 2L], ends_u[, 1L]),
                              edge_betweenness = igraph::edge_betweenness(graph_u, directed = FALSE)),
                   reference = "igraph", covers = "net_edge_betweenness()", quantity = "edge betweenness, unweighted",
                   tolerance = 1e-10, key = c("from", "to"))
})

test_that("bridge centralities equal networktools::bridge()", {
  skip_without_reference("networktools")
  communities <- c(rep(1, 5L), rep(2, 5L))
  as_table <- \(reference) data.frame(
    node = nodes, bridge_strength = reference[["Bridge Strength"]], bridge_betweenness = reference[["Bridge Betweenness"]],
    bridge_closeness = reference[["Bridge Closeness"]], bridge_ei1 = reference[["Bridge Expected Influence (1-step)"]],
    bridge_ei2 = reference[["Bridge Expected Influence (2-step)"]])
  expect_reference(net_bridge(network, communities = communities),
                   as_table(networktools::bridge(network, communities = communities)), reference = "networktools",
                   covers = "net_bridge()", quantity = "bridge centrality, agreeableness and conscientiousness",
                   tolerance = 1e-10, key = "node")
  expect_reference(net_bridge(network, communities = communities, normalize = TRUE),
                   as_table(networktools::bridge(network, communities = communities, normalize = TRUE)),
                   reference = "networktools", covers = "net_bridge()", quantity = "bridge centrality, normalized",
                   tolerance = 1e-10, key = "node")
})

test_that("the observed transitivity and path length of the small-world index equal qgraph and igraph", {
  skip_without_reference(c("qgraph", "igraph"))
  ours <- net_smallworld(network, n_rand = 20L, seed = 1L)
  observed <- qgraph::smallworldness(network, B = 20)
  observed_terms <- unlist(ours[c("transitivity", "aspl")])
  expect_reference(observed_terms, c(transitivity = observed[["trans_target"]], aspl = observed[["averagelength_target"]]),
                   reference = "qgraph", covers = "net_smallworld()", quantity = "observed transitivity and path length",
                   tolerance = 1e-12)
  graph <- igraph::graph_from_adjacency_matrix((abs(network) > 0) * 1, mode = "undirected", diag = FALSE)
  expect_reference(observed_terms, c(transitivity = igraph::transitivity(graph, type = "global"),
                                     aspl = igraph::mean_distance(graph)),
                   reference = "igraph", covers = "net_smallworld()", quantity = "observed transitivity and path length",
                   tolerance = 1e-12)
})
