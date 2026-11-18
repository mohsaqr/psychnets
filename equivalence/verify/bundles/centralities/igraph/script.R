library(igraph)
network <- as.matrix(utils::read.delim("../data/network.tsv", row.names = 1))
graph <- graph_from_adjacency_matrix(abs(network), mode = "undirected", weighted = TRUE)
lengths <- 1 / E(graph)$weight
# expected influence is the signed row sum; igraph has no function of that name
result <- data.frame(node = V(graph)$name, strength = strength(graph),
                     expected_influence = rowSums(network),
                     closeness = closeness(graph, weights = lengths, normalized = TRUE),
                     betweenness = betweenness(graph, weights = lengths, normalized = TRUE))
result
