# Plot methods are drawn for their side effect; we assert they run without error
# on a null graphics device and that the result classes carry a plot method.

mk_data <- function(p = 6, n = 200, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p) %*%
    chol(0.4^abs(outer(seq_len(p), seq_len(p), "-")))
  colnames(x) <- paste0("V", seq_len(p))
  x
}

# Every plot below renders to a throwaway PDF so nothing reaches the screen.
draws_clean <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off())
  expr
  invisible(TRUE)
}

test_that("net_centralities carries a plotting class but stays a data.frame", {
  ct <- net_centralities(ebic_glasso(cor_matrix = cor(mk_data()), n = 200))
  expect_s3_class(ct, "psychnet_centrality")
  expect_s3_class(ct, "data.frame")
  expect_true(is.data.frame(ct))
})

test_that("difference_test carries a plotting class and the observed attribute", {
  skip_slow()
  bs <- net_boot(mk_data(), n_boot = 60, cores = 1)
  dt <- difference_test(bs, type = "strength")
  expect_s3_class(dt, "psychnet_difference")
  expect_s3_class(dt, "data.frame")
  expect_identical(attr(dt, "diff_type"), "strength")
  expect_named(attr(dt, "observed"), bs$node_labels)
})

test_that("plot.psychnet_centrality renders both types without error", {
  ct <- net_centralities(ebic_glasso(cor_matrix = cor(mk_data()), n = 200))
  expect_true(draws_clean(plot(ct)))
  expect_true(draws_clean(plot(ct, type = "line")))
})

test_that("plot.psychnet_bootstrap renders every type without error", {
  skip_slow()
  bs <- net_boot(mk_data(), n_boot = 60, cores = 1, predictability = TRUE)
  expect_true(draws_clean(plot(bs)))
  expect_true(draws_clean(plot(bs, type = "centrality")))
  expect_true(draws_clean(plot(bs, type = "edge_diff")))
  expect_true(draws_clean(plot(bs, type = "centrality_diff", measure = "strength")))
  expect_true(draws_clean(plot(bs, type = "predictability")))
})

test_that("plot.psychnet_bootstrap errors for predictability when absent", {
  skip_slow()
  bs <- net_boot(mk_data(), n_boot = 40, cores = 1)
  pdf(file = tempfile(fileext = ".pdf")); on.exit(dev.off())
  expect_error(plot(bs, type = "predictability"), "predictability")
})

test_that("plot.psychnet_difference renders box and forest styles", {
  skip_slow()
  bs <- net_boot(mk_data(), n_boot = 60, cores = 1)
  expect_true(draws_clean(plot(difference_test(bs, type = "edge"))))
  expect_true(draws_clean(plot(difference_test(bs, type = "strength"))))
  expect_true(draws_clean(plot(difference_test(bs, type = "edge"), style = "forest")))
  expect_true(draws_clean(plot(difference_test(bs, type = "strength"), style = "forest")))
})

test_that("plot.psychnet_stability renders without error", {
  skip_slow()
  st <- net_stability(mk_data(), drop_prop = c(0.3, 0.5, 0.7), iter = 15)
  expect_true(draws_clean(plot(st)))
})

test_that("plot.psychnet_nct renders every type without error", {
  skip_slow()
  cmp <- net_compare(mk_data(seed = 1),
                     mk_data(seed = 2) + 0.2, iter = 60)
  expect_true(draws_clean(plot(cmp)))
  expect_true(draws_clean(plot(cmp, type = "structure")))
  expect_true(draws_clean(plot(cmp, type = "edges")))
})

test_that("plot.psychnet errors cleanly when cograph is unavailable", {
  fit <- ebic_glasso(cor_matrix = cor(mk_data()), n = 200)
  if (!requireNamespace("cograph", quietly = TRUE)) {
    pdf(file = tempfile(fileext = ".pdf")); on.exit(dev.off())
    expect_error(plot(fit), "cograph")
  } else {
    expect_true(draws_clean(plot(fit)))
  }
})
