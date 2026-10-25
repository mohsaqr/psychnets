test_that("undirected rewiring explores the same degree null after relabelling", {
  A <- matrix(0, 6, 6)
  A[1:3, 4:6] <- A[4:6, 1:3] <- 1
  # The 70 labelled simple 3-regular graphs on six vertices consist of
  # 10 K3,3 graphs and 60 triangular prisms (transitivity 0 and 1/3).
  expected <- 2 / 7
  for (o in list(1:6, c(1, 4, 2, 5, 3, 6))) {
    set.seed(1)
    draws <- replicate(150, psychnets:::.psn_rewire(A[o, o], 90L),
                        simplify = FALSE)
    expect_true(all(vapply(draws, function(g)
      isTRUE(all.equal(g, t(g))) && all(diag(g) == 0) &&
        all(rowSums(g) == 3), logical(1))))
    transitivity <- vapply(draws, psychnets:::.psn_transitivity, numeric(1))
    expect_lt(abs(mean(transitivity) - expected), 0.07)
  }
})

test_that("numeric binary MGM prediction preserves training category codes", {
  set.seed(9)
  a <- rbinom(300, 1, 0.5)
  b <- ifelse(runif(300) < 0.9, a, 1 - a)
  raw <- data.frame(a = a + 1, b = b + 1, c = rbinom(300, 1, 0.5) + 1)
  coded <- raw - 1
  for (native in c(TRUE, FALSE)) {
    if (!native && !requireNamespace("glmnet", quietly = TRUE)) next
    fit <- mgm_fit(raw, types = rep("c", 3), threshold = "none",
                   nlambda = 20, native = native)
    pred <- net_predict(fit, raw)
    expected <- psychnets:::.nodewise_predict(fit, coded, encoded = TRUE)
    expect_equal(pred, expected)
    expect_gt(pred$accuracy[1], 0.8)
    expect_equal(fit$nodes$predictability,
                  pmax(0, pmin(1, pred$predictability)))
    factor_data <- raw
    factor_data$a <- factor(raw$a, levels = c(2, 1))
    expect_equal(net_predict(fit, factor_data), pred)
    expect_equal(net_predict(fit, as.matrix(raw)), pred)
    bad <- raw
    bad$a[1] <- 3
    expect_error(net_predict(fit, bad), "unseen category")
  }
})

test_that("nodewise prediction selects source columns independently of labels", {
  set.seed(11)
  d <- data.frame(a = rbinom(200, 1, 0.5), b = rbinom(200, 1, 0.5),
                   c = rbinom(200, 1, 0.5))
  for (estimator in list(ising_fit, ising_sampler)) {
    fit <- estimator(d, labels = c("A", "B", "C"))
    pred <- net_predict(fit, d[c("c", "a", "b")])
    expect_identical(fit$nodes$name, names(d))
    expect_identical(pred$node, c("A", "B", "C"))
    expect_equal(fit$nodes$predictability,
                  pmax(0, pmin(1, pred$predictability)))
    labelled <- d
    names(labelled) <- c("A", "B", "C")
    expect_equal(net_predict(fit, labelled), pred)
    expect_equal(net_predict(fit, unname(as.matrix(d))), pred)
  }
})

test_that("signed NCT plots identify global expected influence", {
  x <- structure(list(
    strength = list(perm = c(0, 0.1, 0.2), observed = 0.1, p_value = 0.5),
    params = list(abs = FALSE)), class = "psychnet_nct")
  hist <- graphics::hist
  axis_label <- NULL
  testthat::local_mocked_bindings(hist = function(..., xlab) {
    axis_label <<- xlab
    hist(..., xlab = xlab)
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(x, type = "strength")
  expect_match(axis_label, "global expected influence")
})
