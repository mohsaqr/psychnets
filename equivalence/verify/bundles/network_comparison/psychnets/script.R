library(psychnets)
group_a <- utils::read.csv("../data/group_a.csv")
group_b <- utils::read.csv("../data/group_b.csv")
set.seed(1)
fit <- net_compare(data1 = group_a, data2 = group_b, iter = 15L, gamma = 0.5, paired = FALSE,
                   abs = TRUE, weighted = TRUE, p_adjust = "none")
fit
tests <- summary(fit)
edges <- as.data.frame(fit)
result <- rbind(
  with(tests, data.frame(quantity = paste("p:", test), kind = "p value", value = p_value)),
  with(edges, data.frame(quantity = paste("p: edge", from, "--", to), kind = "p value", value = p_value)),
  with(edges, data.frame(quantity = paste("weight, group 1:", from, "--", to), kind = "edge weight", value = weight_1)),
  with(edges, data.frame(quantity = paste("weight, group 2:", from, "--", to), kind = "edge weight", value = weight_2)),
  data.frame(quantity = paste("global strength, group", seq_len(2L)), kind = "global strength", value = fit$strength$by_group))
result
