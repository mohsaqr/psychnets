# TMFG has no convex objective: its certificate is structural (3(p-2) edges,
# chordal, connected). LoGo is the chordal Gaussian MRF on the TMFG, certified
# by ggm_support_kkt (precision reproduces S on the TMFG support).

ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

test_that("TMFG has exactly 3(p-2) edges and passes its structural certificate", {
  set.seed(1)
  for (p in c(6, 8, 10)) {
    X <- matrix(stats::rnorm(400 * p), 400, p) %*% chol(ar1(p, 0.5))
    fit <- tmfg_network(X)
    expect_s3_class(fit, "psychnet")
    expect_equal(fit$n_edges, 3L * (p - 2L))
    expect_equal(tmfg_certificate(fit), 0)            # edge count + chordal + connected
  }
})

test_that("the chordality test rejects a 4-cycle", {
  c4 <- matrix(FALSE, 4, 4)
  e <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  for (k in seq_len(nrow(e))) c4[e[k, 1], e[k, 2]] <- c4[e[k, 2], e[k, 1]] <- TRUE
  expect_false(psychnet:::.is_chordal(c4))
  c4[1, 3] <- c4[3, 1] <- TRUE                        # add a chord
  expect_true(psychnet:::.is_chordal(c4))
})

test_that("LoGo reproduces S on the TMFG support and is positive definite", {
  set.seed(2)
  X <- matrix(stats::rnorm(400 * 8), 400, 8) %*% chol(ar1(8, 0.5))
  fit <- logo_network(X)
  expect_lt(fit$kkt, 1e-8)
  expect_equal(ggm_support_kkt(fit$precision, fit$cor_matrix, fit$support),
               fit$kkt)
  expect_true(all(eigen(fit$precision, only.values = TRUE)$values > 0))
  expect_equal(fit$graph, t(fit$graph))
})

test_that("TMFG / LoGo dispatch and reject too-small inputs", {
  set.seed(3)
  X <- matrix(stats::rnorm(300 * 6), 300, 6)
  expect_equal(psychnet(X, "TMFG")$method, "TMFG")
  expect_equal(psychnet(X, "LoGo")$method, "LoGo")
  expect_error(tmfg_network(matrix(stats::rnorm(60), 30, 2)), "at least 4")
})
