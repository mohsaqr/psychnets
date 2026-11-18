library(qgraph)
items <- utils::read.csv("../data/bfi10.csv")
edge_names <- function(m) {
  upper <- upper.tri(m)
  paste(rownames(m)[row(m)[upper]], "--", colnames(m)[col(m)[upper]])
}
correlations <- cor_auto(items, verbose = TRUE)
result <- data.frame(quantity = edge_names(correlations), kind = "polychoric correlation",
                     value = correlations[upper.tri(correlations)])
result
