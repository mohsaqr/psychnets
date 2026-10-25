# Regressions from the complete function review. Each check targets a public
# statistical or data-handling contract rather than duplicating a solver.

test_that("Kendall edge tests use the Kendall null distribution, including ties", {
  set.seed(3)
  x <- rnorm(50); y <- .4 * x + rnorm(50); z <- rnorm(50)
  d <- data.frame(x, y, z)
  fit <- cor_network(d, cor_method = "kendall", alpha = .05)
  expected <- cor.test(x, y, method = "kendall", exact = FALSE)$p.value
  expect_equal(fit$p_values["x", "y"], expected, tolerance = 1e-12)
  expect_lt(expected, .05)
  expect_true(fit$weights["x", "y"] != 0)

  d$x <- round(d$x); d$y[1:3] <- NA_real_
  for (na_method in c("pairwise", "listwise")) {
    fit <- cor_network(d, cor_method = "kendall", alpha = .05,
                       na_method = na_method)
    used <- if (na_method == "listwise") d[complete.cases(d), ] else d
    expected <- cor.test(used$x, used$y, method = "kendall",
                         exact = FALSE)$p.value
    expect_equal(fit$p_values["x", "y"], expected, tolerance = 1e-12)
  }
  expect_error(cor_network(cor_matrix = cor(d, use = "pairwise.complete.obs"),
                           n = 50, cor_method = "kendall", alpha = .05),
               "requires raw")
  expect_error(pcor_network(d, cor_method = "kendall", alpha = .05),
               "partial Kendall")
})

test_that("polyserial correlations are invariant to ordinal code spacing", {
  set.seed(8)
  z <- matrix(rnorm(2000), 1000, 2) %*%
    chol(matrix(c(1, .5, .5, 1), 2))
  d <- data.frame(x = z[, 1],
                  y = as.integer(cut(z[, 2], c(-Inf, -.5, .5, Inf))))
  expected <- cor_auto(d)
  d$y <- 10 * d$y
  expect_equal(cor_auto(d), expected, tolerance = 1e-12)
  d$y[d$y == 20] <- 25
  expect_equal(cor_auto(d), expected, tolerance = 1e-12)
})

test_that("PSD validation uses correlation units for every covariance scale", {
  bad <- matrix(c(1, 1.5, 1.5, 1), 2)
  for (scale in c(1e-12, 1, 1e12)) {
    expect_error(pcor_network(cor_matrix = bad * scale, n = 100), "positive")
    expect_error(relimp_network(cor_matrix = bad * scale), "positive")
  }
  good <- matrix(c(1, .4, .4, 1), 2)
  for (scale in c(1e-12, 1e12)) {
    expect_equal(cor_network(cor_matrix = good * scale)$weights,
                  cor_network(cor_matrix = good)$weights, tolerance = 1e-12)
    expect_equal(relimp_network(cor_matrix = good * scale)$r2,
                  relimp_network(cor_matrix = good)$r2, tolerance = 1e-12)
  }
})

test_that("MGM automatic predictability uses prepared data on its raw scale", {
  set.seed(7)
  f <- rnorm(300)
  d <- data.frame(a = 100 + 15 * (f + rnorm(300)),
                  b = 20 + 3 * (f + rnorm(300)), c = rnorm(300))
  fit <- mgm_fit(d, types = rep("g", 3), threshold = "none", nlambda = 20)
  expected <- pmin(pmax(net_predict(fit, d)$predictability, 0), 1)
  expect_equal(fit$nodes$predictability, expected, tolerance = 1e-12)
  expect_gt(fit$nodes$predictability[1], .2)
  expect_identical(fit$nodewise$source_names, names(d))
})

test_that("MGM records the training categorical levels before numeric encoding", {
  set.seed(9)
  d <- data.frame(a = sample(1:2, 100, TRUE),
                  b = sample(c(-1, 1), 100, TRUE), c = rnorm(100))
  fit <- mgm_fit(d, types = c("c", "c", "g"), nlambda = 10)
  expect_identical(fit$nodewise$category_levels,
                   list(c("1", "2"), c("-1", "1"), NULL))
  skip_if_not_installed("glmnet")
  d$a <- factor(d$a, levels = c("unused", "1", "2"))
  fit <- mgm_fit(d, types = c("c", "c", "g"), native = FALSE, nlambda = 10)
  expect_identical(fit$nodewise$category_levels[[1]], c("1", "2"))
  expect_identical(unname(fit$levels), c(2L, 2L, 1L))
})

test_that("Ising Wald pruning does not certify altered coefficients as the MLE", {
  set.seed(5)
  d <- matrix(rbinom(400 * 4, 1, .5), 400, 4)
  colnames(d) <- letters[1:4]
  fit <- ising_sampler(d, alpha = .01)
  expect_lt(fit$fit_kkt, 1e-6)
  expect_true(is.na(fit$kkt))
  expect_true(is.na(certificate(fit)$certified))
  expect_identical(fit$certificate_target, "pre_pruning_fit")
  full <- ising_sampler(d)
  expect_equal(full$kkt, full$fit_kkt)
  expect_true(certificate(full)$certified)
})

test_that("two-node networks work with the external nodewise engine", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  a <- rbinom(100, 1, .5)
  b <- cbind(a, ifelse(runif(100) < .9, a, 1 - a))
  x <- rnorm(100)
  g <- cbind(x, .8 * x + rnorm(100, sd = .4))
  for (fit in list(ising_fit(b, native = FALSE),
                    mgm_fit(g, native = FALSE, types = c("g", "g")),
                    mgm_fit(b, native = FALSE, types = c("c", "c")))) {
    expect_identical(dim(fit$weights), c(2L, 2L))
    expect_true(all(is.finite(fit$weights)))
    expect_gt(fit$weights[1, 2], 0)
    expect_lt(fit$kkt, 1e-5)
    expect_identical(dim(fit$nodewise$beta_std), c(2L, 2L))
  }
  expect_error(mgm_fit(g, types = c("g", "g"), moderators = 2),
               "at least three variables")
})

test_that("TMFG certificate rejects nonplanar graphs satisfying the old checks", {
  fit <- tmfg_network(cor_matrix = diag(6))
  expect_equal(tmfg_certificate(fit), 0)
  # Triangle plus three common neighbors: contains K3,3, yet is chordal,
  # connected and has exactly 3(p - 2) edges.
  bad <- matrix(FALSE, 6, 6)
  bad[1:3, 1:3] <- TRUE; diag(bad) <- FALSE
  bad[1:3, 4:6] <- bad[4:6, 1:3] <- TRUE
  fit$adjacency <- bad
  expect_gt(tmfg_certificate(fit), 0)
  # K5 with one vertex joined to two clique vertices is another nonplanar
  # chordal graph with the same edge count.
  bad <- matrix(FALSE, 6, 6)
  bad[1:5, 1:5] <- TRUE; diag(bad) <- FALSE
  bad[6, 1:2] <- bad[1:2, 6] <- TRUE
  fit$adjacency <- bad
  expect_gt(tmfg_certificate(fit), 0)
})
