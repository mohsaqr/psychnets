srl <- utils::read.csv("../data/srl.csv")
edge_names <- function(m) {
  upper <- upper.tri(m)
  paste(rownames(m)[row(m)[upper]], "--", colnames(m)[col(m)[upper]])
}
correlations <- cor(srl, method = "pearson")
partials <- -cov2cor(solve(correlations))
result <- rbind(
  data.frame(quantity = paste("cor:", edge_names(correlations)), kind = "correlation",
             value = correlations[upper.tri(correlations)]),
  data.frame(quantity = paste("pcor:", edge_names(partials)), kind = "partial correlation",
             value = partials[upper.tri(partials)]))
result
