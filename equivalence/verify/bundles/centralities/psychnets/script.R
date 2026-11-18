library(psychnets)
network <- as.matrix(utils::read.delim("../data/network.tsv", row.names = 1))
result <- net_centralities(network, measures = c("strength", "expected_influence", "closeness", "betweenness"))
result
