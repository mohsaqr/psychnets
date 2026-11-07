# Regression tests for the Codex review of the 0.7.0 candidate (2026-11-07).
# C1..C9 follow the review's numbering.

test_that("C1 a binary node constant after min_sum is degenerate, not certified", {
  D <- rbind(c(0, 0), c(1, 1), c(1, 1))
  expect_warning(f <- ising_fit(D, min_sum = 2),
                 class = "psychnet_degenerate_sample")
  expect_true(is.na(certificate(f)$certificate))
})

test_that("C2 the weight of a row removed by min_sum does not change the fit", {
  set.seed(81)
  Z <- matrix(stats::rnorm(1200), 300, 4)
  Z[, 2] <- Z[, 1] + 0.7 * Z[, 2]
  B <- rbind(0, (Z > 0) * 1)
  colnames(B) <- paste0("b", 1:4)
  fits <- lapply(c(1, 100), function(v)
    ising_fit(B, weights = c(v, rep(1, 300)), min_sum = 1))
  expect_identical(fits[[1]]$weights, fits[[2]]$weights)
  s1 <- ising_sampler(B, weights = c(1, rep(1, 300)), min_sum = 1)
  s2 <- ising_sampler(B, weights = c(100, rep(1, 300)), min_sum = 1)
  expect_identical(s1$weights, s2$weights)
})

test_that("C3 Ising prediction uses the training coding of factor and logical nodes", {
  b <- as.matrix(expand.grid(a = 0:1, b = 0:1, c = 0:1))
  b <- b[rep(1:8, c(30, 2, 2, 30, 30, 2, 2, 30)), ]
  d <- as.data.frame(b)
  d$a <- factor(ifelse(d$a == 0, "yes", "no"), levels = c("yes", "no"))
  f <- ising_fit(d)
  dr <- d
  dr$a <- factor(dr$a, levels = c("no", "yes"))       # same values, other order
  expect_identical(net_predict(f, d), net_predict(f, dr))
  L <- as.data.frame(b == 1)
  expect_identical(net_predict(ising_fit(L), L),
                   net_predict(ising_fit(b), b))
})

test_that("C4 the constant-vector test is scale-free (CS is unit-invariant)", {
  set.seed(2)
  x <- matrix(stats::rnorm(90), 30, 3)
  run <- function(s) {
    set.seed(1)
    net_stability(x, method = "cor", measures = "custom",
                  centrality_fn = function(W) list(custom = s * (1:3)),
                  drop_prop = 0.2, iter = 3)
  }
  expect_no_warning(small <- run(1e-15))
  expect_identical(as.data.frame(run(1), what = "cs"),
                   as.data.frame(small, what = "cs"))
  expect_true(psychnets:::.psn_is_constant(rep(0, 4)))
  expect_false(psychnets:::.psn_is_constant(1e-15 * (1:3)))
})

test_that("C5 glm_lasso_kkt refuses negative weights and invalid scalars", {
  X <- matrix(c(-1, 0, 1), 3, 1)
  y <- c(0, 1, 0)
  expect_error(glm_lasso_kkt(X, y, b0 = 3, beta = 0, lambda = 0,
                             weights = c(1, -3, 1)),
               class = "psychnet_bad_input")
  expect_error(glm_lasso_kkt(X, y, b0 = 3, beta = 0, lambda = -1),
               class = "psychnet_bad_input")
  expect_error(glm_lasso_kkt(X, c(0, 2, 0), b0 = 0, beta = 0, lambda = 0,
                             family = "binomial"),
               class = "psychnet_bad_input")
})

test_that("C6 character communities keep their alphabetical order", {
  set.seed(99)
  x <- as.data.frame(matrix(stats::rnorm(300), 100, 3))
  expect_identical(names(net_aggregate(x, c("z", "a", "z"))), c("a", "z"))
  expect_identical(names(net_aggregate(x, c(10, 2, 10))), c("2", "10"))
})

test_that("C7 an extreme but finite fit is certified; separation is not", {
  B <- as.matrix(expand.grid(a = 0:1, b = 0:1, c = 0:1))
  w <- ifelse(B[, 1] == B[, 2], 1, 1e-10)
  expect_no_warning(f <- ising_sampler(B, weights = w))
  expect_true(certificate(f)$certified)
  set.seed(1)
  a <- stats::rbinom(200, 1, 0.5)
  X <- cbind(a, b = a, c = stats::rbinom(200, 1, 0.5), d = stats::rbinom(200, 1, 0.5))
  X[1, "b"] <- 1 - X[1, "b"]                           # quasi-complete separation
  expect_warning(g <- ising_sampler(X), class = "psychnet_separation")
  expect_true(is.na(certificate(g)$certificate))
})

test_that("C8 bad choices and flags are classed errors", {
  expect_error(ebic_glasso(diag(3), native = NA), class = "psychnet_bad_input")
  expect_error(cor_network(diag(3), cor_method = "typo"),
               class = "psychnet_bad_input")
  expect_error(ising_fit(dichotomize(SRL_GPT), rule = "XOR"),
               class = "psychnet_bad_input")
})

test_that("C9 a column with community NA is left out before preprocessing", {
  d <- data.frame(a = 1:10, b = 11:20, id = 21:30)
  agg <- net_aggregate(d, c("A", "A", NA))
  expect_identical(names(agg), "A")
  expect_identical(agg, net_aggregate(d[c("a", "b")], c("A", "A")))
  d$id[3] <- NA                                        # an NA in the ID column
  # named form; the NA in the left-out column must not drop a row
  expect_no_warning(agg2 <- net_aggregate(d, c(a = "A", b = "A", id = NA)))
  expect_identical(nrow(agg2), 10L)
})

test_that("net_compare keeps the unadjusted edge p-values next to the adjusted ones", {
  set.seed(1)
  a <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.5^abs(outer(1:5, 1:5, "-")))
  b <- matrix(stats::rnorm(150 * 5), 150, 5) %*% chol(0.2^abs(outer(1:5, 1:5, "-")))
  colnames(a) <- colnames(b) <- paste0("V", 1:5)
  set.seed(2); raw <- as.data.frame(net_compare(a, b, iter = 50))
  expect_identical(raw$p_raw, raw$p_value)              # no adjustment: equal
  for (m in c("holm", "BH")) {
    set.seed(2); fit <- net_compare(a, b, iter = 50, p_adjust = m)
    ed <- as.data.frame(fit)
    expect_identical(ed$p_raw, raw$p_value)             # same permutations
    expect_equal(ed$p_value, stats::p.adjust(ed$p_raw, method = m),
                 tolerance = 1e-15)
    expect_identical(ed$significant, ed$p_value < 0.05) # decided on adjusted p
  }
  set.seed(2); hol <- net_compare(a, b, iter = 50, p_adjust = "holm")
  set.seed(2); none <- net_compare(a, b, iter = 50)
  expect_identical(summary(hol), summary(none))         # M and S never adjusted
})
