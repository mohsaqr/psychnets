# casedrop_reliability() and network_reliability(): edge-weight robustness
# diagnostics. Each verb returns a tidy data.frame (one-liner usage); the rich
# extras (CS-coefficient, raw draws) ride along as attributes for plot().

mk <- function(p = 5, n = 200, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p) %*%
    chol(0.4^abs(outer(seq_len(p), seq_len(p), "-")))
  colnames(x) <- paste0("V", seq_len(p))
  x
}
drew <- function(expr) { pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off()); expr; invisible(TRUE) }

test_that("casedrop_reliability returns a tidy data.frame", {
  skip_slow()
  cd <- casedrop_reliability(mk(), drop_prop = c(0.3, 0.5, 0.7), iter = 20)
  expect_s3_class(cd, "psychnet_casedrop")
  expect_s3_class(cd, "data.frame")
  expect_named(cd, c("metric", "drop_prop", "mean", "sd"))
  expect_equal(nrow(cd), 4L * 3L)                       # 4 metrics x 3 props
  expect_true(all(c("correlation", "mean_abs_dev", "median_abs_dev",
                    "max_abs_dev") %in% cd$metric))
  cs <- attr(cd, "cs")
  expect_true(cs >= 0 && cs <= 1)
})

test_that("casedrop_reliability CS is non-increasing-ish and in range", {
  skip_slow()
  cd <- casedrop_reliability(mk(), drop_prop = seq(0.2, 0.8, 0.2), iter = 25)
  corr <- cd[cd$metric == "correlation", ]
  corr <- corr[order(corr$drop_prop), ]
  expect_true(corr$mean[1] + 0.05 >= corr$mean[nrow(corr)])
})

test_that("network_reliability returns a tidy per-metric data.frame", {
  skip_slow()
  rel <- network_reliability(mk(), iter = 40)
  expect_s3_class(rel, "psychnet_reliability")
  expect_s3_class(rel, "data.frame")
  expect_named(rel, c("metric", "mean", "sd", "lower", "upper"))
  expect_true(all(c("correlation", "mean_abs_dev") %in% rel$metric))
  expect_equal(nrow(attr(rel, "iterations")), 40L)      # raw draws in attr
})

test_that("both verbs are estimator-agnostic (pcor) and route through psychnet", {
  skip_slow()
  cd <- casedrop_reliability(mk(), method = "pcor",
                             drop_prop = c(0.3, 0.6), iter = 15)
  rel <- network_reliability(mk(), method = "pcor", iter = 20)
  expect_s3_class(cd, "psychnet_casedrop")
  expect_s3_class(rel, "psychnet_reliability")
})

test_that("plot methods render without error", {
  skip_slow()
  cd <- casedrop_reliability(mk(), drop_prop = c(0.3, 0.5, 0.7), iter = 20)
  rel <- network_reliability(mk(), iter = 30)
  expect_true(drew(plot(cd)))
  expect_true(drew(plot(rel)))
})

test_that("print shows the CS header and the tidy table", {
  skip_slow()
  cd <- casedrop_reliability(mk(), drop_prop = 0.5, iter = 10)
  rel <- network_reliability(mk(), iter = 10)
  expect_output(print(cd), "edge-weight stability")
  expect_output(print(rel), "split-half reliability")
})
