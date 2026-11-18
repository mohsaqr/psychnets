library(psychnets)
items <- utils::read.csv("../data/bfi10.csv")
fit <- ebic_glasso(data = items, gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,
                   threshold = 0, cor_method = "pearson", penalize_diagonal = FALSE,
                   native = TRUE, refit = TRUE)
fit
edges <- as.data.frame(fit, include_zero = TRUE)
result <- rbind(
  with(edges, data.frame(quantity = paste("weight:", from, "--", to), kind = "partial correlation", value = weight)),
  data.frame(quantity = "selected lambda", kind = "penalty", value = fit$lambda))
result
