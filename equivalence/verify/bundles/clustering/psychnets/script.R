library(psychnets)
network <- as.matrix(utils::read.delim("../data/network.tsv", row.names = 1))
result <- net_clustering(network)
result
