library(qgraph)
items <- utils::read.csv("../data/bfi10.csv")
edge_names <- function(m) {
  upper <- upper.tri(m)
  paste(rownames(m)[row(m)[upper]], "--", colnames(m)[col(m)[upper]])
}
fit <- EBICglasso(cor(items), n = nrow(items), gamma = 0.5, nlambda = 100, lambda.min.ratio = 0.01,
                  penalize.diagonal = FALSE, refit = FALSE, threshold = FALSE,
                  returnAllResults = TRUE, thr = 1e-10, maxit = 1e5)
network <- fit$optnet
result <- rbind(
  data.frame(quantity = paste("weight:", edge_names(network)), kind = "partial correlation",
             value = network[upper.tri(network)]),
  data.frame(quantity = "selected lambda", kind = "penalty", value = fit$lambda[which.min(fit$ebic)]))
result
