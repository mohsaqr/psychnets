library(psychnets)
answers <- utils::read.csv("../data/bfi6_binary.csv")
fit <- ising_sampler(data = answers, rule = "AND", alpha = NULL)
fit
result <- rbind(
  with(as.data.frame(fit), data.frame(quantity = paste("weight:", from, "--", to), kind = "edge weight", value = weight)),
  data.frame(quantity = paste("threshold:", names(fit$thresholds)), kind = "threshold", value = unname(fit$thresholds)))
result
