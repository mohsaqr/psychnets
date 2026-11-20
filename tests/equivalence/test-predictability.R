# net_predict() and node_predictability() for Gaussian networks: R-squared of each node given the others,
# 1 - 1 / (Theta_jj S_jj) (Haslbeck and Waldorp 2018). For the unregularized partial-correlation network
# this is exactly the R-squared of the linear regression of the node on all other nodes, which stats::lm()
# computes; for a network that keeps only some edges, the R-squared of the regression on the node's
# neighbours under the model's precision is the same formula. Ising predictability is compared with the
# classification accuracy of the fitted logistic regressions written out in base R.

items <- utils::read.csv(reference_data("data/bfi10.csv"))
nodes <- names(items)

test_that("predictability of the partial-correlation network equals stats::lm()'s R-squared", {
  r_squared <- vapply(nodes, \(node) summary(stats::lm(stats::reformulate(setdiff(nodes, node), node),
                                                        data = items))$r.squared, numeric(1))
  fit <- pcor_network(data = items)
  expect_reference(net_predict(fit), data.frame(node = nodes, predictability = unname(r_squared)), reference = "stats",
                   covers = "net_predict()", quantity = "R-squared of each node on all others", tolerance = 1e-10,
                   key = "node")
  expect_reference(node_predictability(fit), r_squared, reference = "stats", covers = "node_predictability()",
                   quantity = "R-squared of each node on all others, as a vector", tolerance = 1e-10)
})

test_that("predictability of an Ising model equals the accuracy of its nodewise logistic regressions", {
  binary <- utils::read.csv(reference_data("data/bfi6_binary.csv"))
  fit <- ising_sampler(binary, rule = "AND", alpha = NULL)
  # Each node's own logistic regression on all others by stats::glm(); predict 1 when the fitted probability
  # exceeds 1/2. CC is the share predicted correctly; nCC = (CC - m) / (1 - m), m the share of the more
  # frequent category (Haslbeck and Waldorp 2018).
  accuracy <- vapply(names(binary), \(node) {
    model <- stats::glm(stats::reformulate(setdiff(names(binary), node), node), family = stats::binomial(),
                        data = binary, control = stats::glm.control(epsilon = 1e-12, maxit = 100))
    mean((stats::fitted(model) > 0.5) == binary[[node]])
  }, numeric(1))
  majority <- vapply(binary, \(x) max(mean(x), 1 - mean(x)), numeric(1))
  expected <- data.frame(node = names(binary), predictability = (accuracy - majority) / (1 - majority),
                         accuracy = accuracy)
  # A share of 2,632 predictions: equal unless a fitted probability sits within solver noise of 1/2.
  expect_reference(net_predict(fit, data = binary), expected,
                   reference = "independent: classification by each node's logistic regression, stats::glm()",
                   covers = "net_predict()", quantity = "nCC and CC of the unregularized Ising model",
                   tolerance = 1e-12, key = "node")
})
