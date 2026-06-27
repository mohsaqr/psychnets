# Regression tests for the adversarial-review fixes (2026-06-22).

mk <- function(seed, n = 200, p = 5) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(0.4^abs(outer(1:p, 1:p, "-")))
  colnames(X) <- paste0("V", seq_len(p))
  X
}

test_that(".new_psychnet stores a zero-diagonal graph", {
  # a correlation-like matrix with unit diagonal must not leak the diagonal
  g <- matrix(0.5, 4, 4); diag(g) <- 1
  obj <- psychnet:::.new_psychnet(g, paste0("V", 1:4), "test", FALSE, 10)
  expect_true(all(diag(obj$graph) == 0))
  expect_equal(centrality(obj)$strength, rep(1.5, 4))   # 3 * 0.5, no diagonal
})

test_that(".new_psychnet rejects a label/dimension mismatch", {
  expect_error(psychnet:::.new_psychnet(matrix(0, 3, 3), c("a", "b"), "t", FALSE, 5),
               "labels length")
})

test_that("glasso_kkt is infinite for an indefinite precision matrix", {
  theta <- matrix(c(1, 2, 2, 1), 2, 2)                  # eigenvalues 3, -1
  expect_true(is.infinite(glasso_kkt(theta, diag(2), 0.1)))
  # a genuine optimum still certifies near zero
  S <- 0.5^abs(outer(1:5, 1:5, "-"))
  fit <- ebic_glasso(cor_matrix = S, n = 250)
  expect_lt(fit$kkt, 1e-7)
})

test_that("glm_lasso_kkt penalizes a non-stationary intercept", {
  set.seed(1)
  # zero predictors isolate the intercept: the slope gradient is exactly 0, so
  # the certificate reflects only mean(y - mu), the intercept stationarity.
  X <- matrix(0, 100, 2)
  y <- as.numeric(stats::rnorm(100) > 0)
  good <- glm_lasso_kkt(X, y, stats::qlogis(mean(y)), c(0, 0), 0, "binomial")
  bad  <- glm_lasso_kkt(X, y, stats::qlogis(mean(y)) + 1, c(0, 0), 0, "binomial")
  expect_lt(good, 1e-8)
  expect_gt(bad, 0.1)
})

test_that("an isolated binary node gets a logit-scale threshold", {
  set.seed(3)
  # V1 independent of the rest; V2..V5 a chain
  base <- mk(3, n = 600, p = 4)
  b <- cbind(V0 = as.numeric(stats::rnorm(600) > 0.4), (base > 0) * 1)
  colnames(b) <- paste0("V", 0:4)
  fit <- ising_fit(b)
  # threshold for the (near) isolated node should track qlogis(prevalence),
  # i.e. be negative for a low-prevalence node, never a probability in (0,1)
  expect_true(fit$thresholds[["V0"]] < 0)
})

test_that("centrality_stability CS is invariant to drop_prop ordering", {
  set.seed(7)
  x <- mk(8, n = 250, p = 5)
  a <- centrality_stability(x, drop_prop = c(0.1, 0.3, 0.5), iter = 25)
  b <- centrality_stability(x, drop_prop = c(0.5, 0.1, 0.3), iter = 25)
  expect_equal(sort(a$cs), sort(b$cs))
})

test_that("bootstrap keeps all directed edges for relimp", {
  set.seed(9)
  x <- mk(9, n = 150, p = 4)
  bs <- bootstrap_network(x, method = "relimp", n_boot = 20, cores = 1)
  expect_equal(nrow(bs$edges), 4 * 3)                   # p*(p-1) directed edges
})

test_that("paired nct requires equal group sizes", {
  expect_error(nct(mk(1, 10, 3), mk(2, 12, 3), iter = 5, paired = TRUE),
               "equal size")
})

test_that("correlation significance does not fabricate edges when df <= 0", {
  # n - 2 - (p-2) = n - p <= 0 for the partial test: nothing can be significant
  set.seed(2)
  x <- mk(2, n = 6, p = 6)
  fit <- pcor_network(x, alpha = 0.05)
  expect_equal(fit$n_edges, 0L)
})

test_that("centrality rejects a non-square bare matrix", {
  expect_error(centrality(matrix(1, 2, 3)), "square")
})

test_that("relimp rejects an indefinite cor_matrix", {
  S <- matrix(c(1, .9, .9, .9, 1, -.9, .9, -.9, 1), 3, 3)
  expect_error(relimp_network(cor_matrix = S), "positive-semidefinite")
})

# ---- round 2 ---------------------------------------------------------------

test_that("the empty graph is returned for zero-association input, not an error", {
  for (fn in list(ebic_glasso, huge_network)) {
    fit <- fn(cor_matrix = diag(5), n = 200)
    expect_equal(fit$n_edges, 0L)
    expect_lt(fit$kkt, 1e-8)
  }
  expect_equal(ggm_modselect(cor_matrix = diag(5), n = 200)$n_edges, 0L)
  # the empty-graph EBIC equals the textbook n*p on a correlation matrix
  expect_equal(ebic_glasso(cor_matrix = diag(4), n = 100)$ebic, 400)
})

test_that("glasso_kkt rejects asymmetric and singular precision matrices", {
  expect_true(is.infinite(glasso_kkt(matrix(c(1, .5, 0, 1), 2, 2), diag(2), .1)))
  expect_true(is.infinite(ggm_support_kkt(matrix(c(1, .5, 0, 1), 2, 2), diag(2),
                                          matrix(TRUE, 2, 2))))
})

test_that("relimp accepts a covariance-scale matrix via cov2cor", {
  S <- matrix(c(4, 1, 1, 1), 2, 2)               # var 4 and 1, cov 1
  fit <- relimp_network(cor_matrix = S)
  expect_true(all(fit$r2 >= 0 & fit$r2 <= 1))    # standardized R^2, not 1.0+
  expect_lt(lmg_certificate(fit), 1e-8)
})

test_that("nct tolerates a constant column and rejects mismatched columns", {
  set.seed(5)
  a <- cbind(x = stats::rnorm(60), y = stats::rnorm(60), z = 1)  # z constant
  b <- cbind(x = stats::rnorm(60), y = stats::rnorm(60), z = 1)
  expect_s3_class(nct(a, b, iter = 5), "psychnet_nct")           # no crash
  d <- cbind(y = stats::rnorm(60), x = stats::rnorm(60), z = 1)  # reordered
  expect_error(nct(a, d, iter = 5), "same columns")
})

test_that("predictability rejects data missing network-node columns", {
  set.seed(6)
  b <- (matrix(stats::rnorm(400 * 3), 400, 3) > 0) * 1L
  colnames(b) <- c("A", "B", "C")
  fit <- ising_fit(b)
  wrong <- b; colnames(wrong) <- c("A", "B", "Z")
  expect_error(predictability(fit, data = wrong), "missing columns")
})

test_that("logo errors on a singular clique block instead of returning garbage", {
  expect_error(logo_network(cor_matrix = matrix(1, 5, 5), n = 100), "singular")
})

# ---- argument tidy pass ----------------------------------------------------

test_that("cor / pcor accept a precomputed cor_matrix like the other GGMs", {
  set.seed(8)
  x <- mk(8, n = 300, p = 5)
  S <- stats::cor(x)
  expect_equal(cor_network(cor_matrix = S)$graph, cor_network(x)$graph)
  expect_equal(pcor_network(cor_matrix = S)$graph, pcor_network(x)$graph)
  # n is only required for significance testing
  expect_error(cor_network(cor_matrix = S, alpha = 0.05), "`n` is required")
})

test_that("cor_method is the correlation-type argument and reaches the front door", {
  set.seed(9)
  x <- mk(9, n = 300, p = 5)
  # spearman selected through psychnet() forwards to the verb
  viafd  <- psychnet(x, "pcor", cor_method = "spearman")$graph
  direct <- pcor_network(x, cor_method = "spearman")$graph
  expect_equal(viafd, direct)
  # and it actually changes the result vs pearson
  expect_false(isTRUE(all.equal(direct, pcor_network(x, cor_method = "pearson")$graph)))
  # the old overloaded `method=` name is gone from the verb
  expect_false("method" %in% names(formals(cor_network)))
  expect_true("cor_method" %in% names(formals(tmfg_network)))
})

test_that("psychnet does not override ising/mgm native gamma (0.25)", {
  set.seed(10)
  b <- (mk(10, n = 600, p = 5) > 0) * 1L
  colnames(b) <- paste0("V", 1:5)
  expect_equal(psychnet(b, "ising")$graph, ising_fit(b)$graph)        # 0.25
  expect_equal(psychnet(b, "ising", gamma = 0.5)$graph,
               ising_fit(b, gamma = 0.5)$graph)                               # override
})

# ---- lean netobject-compatible output shape --------------------------------

test_that("a fitted network has the lean netobject shape and cograph class", {
  fit <- ebic_glasso(cor_matrix = 0.4^abs(outer(1:5, 1:5, "-")), n = 250)
  expect_s3_class(fit, "cograph_network")
  # canonical (str-visible) fields, netobject-aligned
  expect_true(all(c("weights", "nodes", "edges", "directed", "method", "n") %in%
                    names(fit)))
  expect_true(is.matrix(fit$weights) && all(diag(fit$weights) == 0))
  expect_named(fit$nodes, c("id", "label", "name"))
  expect_named(fit$edges, c("from", "to", "weight"))
  expect_identical(as.data.frame(fit), fit$edges)
  # derivable counts are no longer stored as fields...
  expect_false(any(c("graph", "n_nodes", "n_edges", "n_obs") %in% names(fit)))
  # ...but the legacy accessors still resolve via the alias
  expect_identical(fit$graph, fit$weights)
  expect_equal(fit$n_nodes, 5L)
  expect_equal(fit$n_edges, nrow(fit$edges))
  expect_equal(fit$n_obs, fit$n)
})

# ---- closed reference gaps -------------------------------------------------

test_that("mgm_fit gains an AND/OR rule (OR keeps at least as many edges)", {
  set.seed(11)
  f <- matrix(stats::rnorm(2000 * 3), 2000, 3) %*% chol(0.5^abs(outer(1:3, 1:3, "-")))
  d <- data.frame(g1 = f[, 1], b1 = (f[, 2] > 0) * 1, b2 = (f[, 3] > 0) * 1)
  expect_gte(mgm_fit(d, rule = "OR")$n_edges, mgm_fit(d, rule = "AND")$n_edges)
})

test_that("min_sum drops low sum-score rows for the Ising fits", {
  set.seed(12)
  b <- (mk(12, n = 500, p = 5) > 0) * 1L
  colnames(b) <- paste0("V", 1:5)
  expect_s3_class(ising_fit(b, min_sum = 2), "psychnet")
  expect_s3_class(ising_sampler(b, min_sum = 2), "psychnet")
  expect_error(ising_fit(b, min_sum = 99), "removed every row")
})

test_that("observation weights are honoured and stay self-certified", {
  set.seed(13)
  z <- matrix(stats::rnorm(800 * 2), 800, 2)
  x <- cbind(z[, 1], z[, 1], z[, 2], z[, 2]) + matrix(stats::rnorm(800 * 4), 800)
  b <- (x > 0) * 1L; colnames(b) <- paste0("V", 1:4)
  w <- stats::runif(nrow(b), 0.5, 2)
  # unit weights reproduce the unweighted fit exactly
  expect_equal(ising_fit(b)$graph, ising_fit(b, weights = rep(1, nrow(b)))$graph)
  # non-trivial weights change the fit but the KKT certificate still holds
  fw <- ising_fit(b, weights = w)
  expect_lt(fw$kkt, 1e-6)
  expect_false(isTRUE(all.equal(fw$graph, ising_fit(b)$graph)))
  # mgm and the unregularized sampler too
  f <- data.frame(g1 = z[, 1], b1 = b[, 1], b2 = b[, 3])
  expect_lt(mgm_fit(f, weights = stats::runif(800, 0.5, 2))$kkt, 1e-6)
  expect_lt(ising_sampler(b, weights = w)$kkt, 1e-6)
  # bad weights are rejected
  expect_error(ising_fit(b, weights = rep(1, 3)), "one value per row")
  expect_error(mgm_fit(f, weights = rep(-1, 800)), "non-negative")
})
