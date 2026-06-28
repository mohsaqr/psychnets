# dichotomize: Likert/continuous -> 0/1 for Ising networks.

test_that("dichotomize returns a 0/1 matrix with the same shape and names", {
  b <- dichotomize(SRL_GPT, method = "median")
  expect_true(is.matrix(b))
  expect_equal(dim(b), dim(as.matrix(SRL_GPT)))
  expect_identical(colnames(b), names(SRL_GPT))
  expect_true(all(b %in% c(0L, 1L)))
})

test_that("each method splits as documented", {
  set.seed(1)
  x <- matrix(c(1, 2, 2, 3, 5, 5), ncol = 1)
  expect_equal(as.vector(dichotomize(x, "median")), as.integer(x >= stats::median(x)))
  expect_equal(as.vector(dichotomize(x, "mean")),   as.integer(x > mean(x)))
  # rank gives a balanced split (half ones), robust to ties
  rk <- dichotomize(matrix(rep(1:2, each = 50), ncol = 1), "rank")
  expect_equal(sum(rk), 50L)
})

test_that("dichotomized Likert data feeds ising_fit", {
  b <- dichotomize(SRL_GPT, method = "rank")
  fit <- ising_fit(b)
  expect_s3_class(fit, "psychnet")
  expect_lt(fit$kkt, 1e-5)
})

test_that("dichotomize validates input", {
  expect_error(dichotomize("a"), "numeric")
  expect_error(dichotomize(SRL_GPT, method = "tertile"))   # bad method
})
