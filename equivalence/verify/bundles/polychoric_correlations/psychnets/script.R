library(psychnets)
items <- utils::read.csv("../data/bfi10.csv")
cor_auto(data = items)
result <- with(as.data.frame(cor_network(data = items, cor_method = "auto")),
               data.frame(quantity = paste(from, "--", to), kind = "polychoric correlation", value = weight))
result
