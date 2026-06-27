# LMG/Shapley shares satisfy the efficiency identity: incoming shares per node
# sum to that node's full-model R^2 (lmg_certificate).

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("relimp shares decompose each node's R^2 exactly", {
  for (S in list(ar1(5, 0.4), ar1(6, 0.5))) {
    fit <- relimp_network(cor_matrix = S)
    expect_s3_class(fit, "psychnet")
    expect_true(fit$directed)
    expect_lt(lmg_certificate(fit), 1e-8)
    expect_equal(unname(colSums(fit$graph)), unname(fit$r2), tolerance = 1e-8)
  }
})

test_that("the importance network is directed and non-negative", {
  S <- ar1(5, 0.5)
  fit <- relimp_network(cor_matrix = S)
  expect_false(isTRUE(all.equal(fit$graph, t(fit$graph))))  # asymmetric
  expect_true(all(fit$graph >= -1e-12))                     # LMG shares >= 0
  ed <- as.data.frame(fit)
  expect_named(ed, c("from", "to", "weight"))
})

test_that("relimp refuses too many nodes", {
  S <- diag(25)
  expect_error(relimp_network(cor_matrix = S, max_nodes = 21L), "refuses")
})

test_that("relimp dispatches via psychnet", {
  set.seed(1)
  X <- matrix(stats::rnorm(300 * 5), 300, 5) %*% chol(ar1(5, 0.4))
  expect_equal(psychnet(X, "relimp")$method, "relimp")
})
