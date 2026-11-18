library(qgraph)
network <- as.matrix(utils::read.delim("../data/network.tsv", row.names = 1))
coefficients <- clustcoef_auto(network)
result <- data.frame(node = rownames(coefficients), coefficients)
result
