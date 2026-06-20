# Unregularized Ising: the certificate is the maximum-likelihood score residual
# (glm_lasso_kkt at lambda = 0).

make_binary <- function(n, seed) {
  set.seed(seed)
  z <- matrix(stats::rnorm(n * 2), n, 2)
  x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(n * 4), n)
  b <- (x > 0) * 1
  colnames(b) <- paste0("V", 1:4)
  b
}

test_that("ising_sampler certifies the MLE score residual", {
  b <- make_binary(600, 1)
  fit <- ising_sampler(b)
  expect_s3_class(fit, "psychnet")
  expect_lt(fit$kkt, 1e-6)
  expect_equal(fit$graph, t(fit$graph))
  expect_true(all(diag(fit$graph) == 0))
})

test_that("Wald pruning removes edges and recovers the two clusters", {
  b <- make_binary(800, 2)
  full <- ising_sampler(b)
  pruned <- ising_sampler(b, alpha = 0.01, adjust = "BH")
  expect_lte(pruned$n_edges, full$n_edges)
  # the within-cluster edges V1-V2 and V3-V4 should survive
  expect_true(abs(pruned$graph["V1", "V2"]) > 0)
  expect_true(abs(pruned$graph["V3", "V4"]) > 0)
})

test_that("OR keeps at least as many edges as AND", {
  b <- make_binary(600, 3)
  and <- ising_sampler(b, rule = "AND")
  or  <- ising_sampler(b, rule = "OR")
  expect_gte(or$n_edges, and$n_edges)
})

test_that("ising_sampler dispatches via estimate_network", {
  b <- make_binary(500, 4)
  expect_equal(estimate_network(b, "IsingSampler")$method, "IsingSampler")
  expect_equal(estimate_network(b, "isingsampler")$method, "IsingSampler")
})
