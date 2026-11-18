library(psychnets)
srl <- utils::read.csv("../data/srl.csv")
correlations <- as.data.frame(cor_network(data = srl, cor_method = "pearson"))
partials <- as.data.frame(pcor_network(data = srl, cor_method = "pearson"))
correlations
partials
result <- rbind(
  with(correlations, data.frame(quantity = paste("cor:", from, "--", to), kind = "correlation", value = weight)),
  with(partials, data.frame(quantity = paste("pcor:", from, "--", to), kind = "partial correlation", value = weight)))
