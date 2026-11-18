answers <- as.matrix(utils::read.csv("../data/bfi6_binary.csv"))
edge_names <- function(m) {
  upper <- upper.tri(m)
  paste(rownames(m)[row(m)[upper]], "--", colnames(m)[col(m)[upper]])
}
items <- colnames(answers)
# one column per regressed item: its intercept, then its coefficient on every item (0 on itself)
nodewise <- sapply(items, function(item) {
  others <- setdiff(items, item)
  fit <- glm(answers[, item] ~ answers[, others], family = binomial,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  replace(setNames(numeric(length(items) + 1), c("threshold", items)), c("threshold", others), coef(fit))
})
nodewise
weights <- (nodewise[items, ] + t(nodewise[items, ])) / 2
result <- rbind(
  data.frame(quantity = paste("weight:", edge_names(weights)), kind = "edge weight", value = weights[upper.tri(weights)]),
  data.frame(quantity = paste("threshold:", items), kind = "threshold", value = nodewise["threshold", ]))
result
