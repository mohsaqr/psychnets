library(qgraph)
network <- as.matrix(utils::read.delim("../data/network.tsv", row.names = 1))
central <- centrality_auto(network)$node.centrality
p <- nrow(network)
result <- data.frame(node = rownames(central), strength = central$Strength,
                     expected_influence = central$ExpectedInfluence,
                     closeness = central$Closeness * (p - 1),
                     betweenness = central$Betweenness / ((p - 1) * (p - 2) / 2))
result
