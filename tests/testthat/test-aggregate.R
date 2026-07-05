# net_aggregate(): collapse communities into super-nodes. Score methods return a
# reduced data.frame; association methods return a macro psychnet.

agg_data <- function(seed = 1, n = 200) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::rnorm(n * 5), n, 5))
  names(d) <- paste0("V", 1:5)
  d
}
comm5 <- c(1, 1, 2, 2, 2)

test_that("score methods return a reduced data.frame, one column per community", {
  d <- agg_data()
  for (m in c("mean", "median", "sum", "pca", "factor", "loadings")) {
    r <- net_aggregate(d, comm5, method = m)
    expect_s3_class(r, "data.frame")
    expect_equal(dim(r), c(200L, 2L))
    expect_named(r, c("1", "2"))
  }
})

test_that("mean composite equals rowMeans of z-scored members", {
  d <- agg_data()
  Z <- scale(as.matrix(d))
  manual <- rowMeans(Z[, 1:2])
  ours <- net_aggregate(d, comm5, method = "mean")[["1"]]
  expect_equal(ours, manual, ignore_attr = TRUE)
})

test_that("a singleton community passes its z-scored item through", {
  d <- agg_data()
  r <- net_aggregate(d, c(1, 2, 3, 3, 3), method = "mean")
  expect_equal(r[["1"]], as.numeric(scale(d$V1)), ignore_attr = TRUE)
})

test_that("pca composite is sign-aligned with its members", {
  d <- agg_data()
  pc <- net_aggregate(d, comm5, method = "pca")[["1"]]
  expect_gt(mean(stats::cor(pc, scale(as.matrix(d))[, 1:2])), 0)
})

test_that("association methods return a symmetric macro psychnet", {
  d <- agg_data()
  for (m in c("average", "rv", "canonical")) {
    r <- net_aggregate(d, comm5, method = m)
    expect_s3_class(r, "psychnet")
    expect_equal(nrow(r$nodes), 2L)
    expect_equal(r$weights, t(r$weights))
    expect_equal(diag(r$weights), c(0, 0), ignore_attr = TRUE)
  }
})

test_that("canonical macro edge equals stats::cancor; rv is in [0,1]", {
  d <- agg_data(); m <- as.matrix(d)
  ca <- net_aggregate(d, comm5, method = "canonical")$weights[1, 2]
  expect_equal(ca, stats::cancor(m[, 1:2], m[, 3:5])$cor[1])
  rv <- net_aggregate(d, comm5, method = "rv")$weights[1, 2]
  expect_true(rv >= 0 && rv <= 1)
})

test_that("reduced data flows back into psychnet() for a macro network", {
  red <- net_aggregate(agg_data(), comm5, method = "mean")
  macro <- psychnet(red, method = "pcor")
  expect_s3_class(macro, "psychnet")
  expect_equal(nrow(macro$nodes), 2L)
})

test_that("net_aggregate validates community length and alignment", {
  d <- agg_data()
  expect_error(net_aggregate(d, c(1, 1, 2), method = "mean"), "one entry per item")
  named <- stats::setNames(comm5, names(d))
  expect_equal(net_aggregate(d, named)[["1"]], net_aggregate(d, comm5)[["1"]])
})
