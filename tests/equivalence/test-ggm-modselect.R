# ggm_modselect() and ggm_support_kkt(): the unregularized Gaussian graphical model selected by the EBIC over
# the graphical lasso path, with and without stepwise search, against qgraph::ggmModSelect(); the
# support-restricted KKT residual against the maximum-likelihood conditions written out in base R.

items <- utils::read.csv(reference_data("data/bfi10.csv"))
correlations <- stats::cor(items)
n_rows <- nrow(items)

test_that("the selected unregularized network equals qgraph::ggmModSelect()", {
  skip_without_reference("qgraph")
  # qgraph fits each candidate graph by maximum likelihood with glasso at rho = 0 (thr = 1e-4 inside
  # ggmModSelect, not reachable from its arguments); psychnets solves the same constrained ML exactly.
  # 1e-4 allows for glasso's stopping rule.
  stepwise <- qgraph::ggmModSelect(correlations, n_rows, gamma = 0, stepwise = TRUE, verbose = FALSE)
  expect_reference(as.data.frame(ggm_modselect(data = items, gamma = 0, stepwise = TRUE), include_zero = TRUE),
                   ref_edges(stepwise$graph), reference = "qgraph", covers = "ggm_modselect()",
                   quantity = "partial correlation, stepwise, gamma = 0", tolerance = 1e-4, key = c("from", "to"))
  no_steps <- qgraph::ggmModSelect(correlations, n_rows, gamma = 0.5, stepwise = FALSE, verbose = FALSE)
  expect_reference(as.data.frame(ggm_modselect(data = items, gamma = 0.5, stepwise = FALSE), include_zero = TRUE),
                   ref_edges(no_steps$graph), reference = "qgraph", covers = "ggm_modselect()",
                   quantity = "partial correlation, no stepwise search, gamma = 0.5", tolerance = 1e-4,
                   key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(items, method = "ggm"), include_zero = TRUE), ref_edges(stepwise$graph),
                   reference = "qgraph", covers = "psychnet()", quantity = "partial correlation, psychnet(method = \"ggm\")",
                   tolerance = 1e-4, key = c("from", "to"))
})

test_that("ggm_support_kkt() equals the maximum-likelihood conditions of a GGM on a graph", {
  fit <- ggm_modselect(data = items, gamma = 0, stepwise = TRUE)
  theta <- fit$precision
  support <- fit$support
  # Dempster (1972): W = theta^-1 equals S on the diagonal and on every edge of the graph, and theta is zero
  # off the graph. The residual is the largest violation.
  residual <- \(theta, S, support) {
    W <- solve(theta)
    off <- upper.tri(theta)
    max(abs(diag(W) - diag(S)), abs((W - S)[off][support[off]]), abs(theta[off][!support[off]]))
  }
  expect_reference(ggm_support_kkt(theta, correlations, support), residual(theta, correlations, support),
                   reference = "independent: Dempster's ML conditions for a covariance selection model, base R",
                   covers = c("ggm_support_kkt()", "$.psychnet()"), quantity = "support KKT residual at the selected model",
                   tolerance = 1e-12, note = "Both values are near machine precision; absolute 1e-12 compares them.")
  # The same residual at a precision that is NOT the optimum: the unconstrained inverse restricted to a graph
  # without its weakest edge.
  perturbed <- solve(correlations)
  weakest <- which(abs(perturbed) == min(abs(perturbed[upper.tri(perturbed)])), arr.ind = TRUE)
  graph <- matrix(TRUE, nrow(perturbed), ncol(perturbed))
  graph[weakest] <- FALSE
  graph[weakest[, 2:1, drop = FALSE]] <- FALSE   # an edge is a property of the unordered pair
  expect_reference(ggm_support_kkt(perturbed, correlations, graph), residual(perturbed, correlations, graph),
                   reference = "independent: Dempster's ML conditions for a covariance selection model, base R",
                   covers = "ggm_support_kkt()", quantity = "support KKT residual off the optimum", tolerance = 1e-10,
                   measure = "relative")
  expect_reference(certificate(fit)$certificate, residual(theta, correlations, support),
                   reference = "independent: Dempster's ML conditions for a covariance selection model, base R",
                   covers = "certificate()", quantity = "certificate of ggm_modselect()", tolerance = 1e-12)
})
