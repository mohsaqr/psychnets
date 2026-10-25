# Regression checks for the function-by-function framework review.

framework_data <- function(n = 60L, seed = 11L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(x = x, y = x + stats::rnorm(n, sd = 0.2), z = stats::rnorm(n))
}

test_that("custom resampling labels preserve networks when ids are ignored", {
  d <- data.frame(id = paste0("s", 1:60), framework_data())
  set.seed(4)
  plain <- net_boot(d, method = "pcor", n_boot = 3, cores = 1)
  set.seed(4)
  named <- net_boot(d, method = "pcor", n_boot = 3, cores = 1,
                    labels = c("A", "B", "C"))
  expect_equal(unname(named$observed$weights), unname(plain$observed$weights))
  expect_equal(named$edge_boot, plain$edge_boot)
  expect_identical(named$observed$nodes$label, c("A", "B", "C"))
  expect_identical(named$observed$nodes$name, c("x", "y", "z"))
  expect_gt(max(abs(named$observed$weights)), 0.5)
  selected <- net_boot(d, method = "pcor", n_boot = 2, cores = 1,
    vars = c("x", "y"), labels = c("A", "B"))
  expect_equal(unname(selected$observed$weights),
                 unname(pcor_network(d[c("x", "y")])$weights))
  expect_identical(selected$observed$nodes$name, c("x", "y"))
  # A node becoming constant in a draw keeps its label and an isolated row.
  d$z <- 0
  fit <- psychnets:::.psn_refit_network(d, "pcor", c("A", "B", "C"))
  expect_identical(dim(fit$weights), c(3L, 3L))
  expect_true(all(fit$weights["C", ] == 0))
})

test_that("two-node reliability returns deviations and unavailable correlation", {
  d <- framework_data()[1:2]
  split <- net_split_reliability(d, method = "pcor", iter = 2)
  drop <- net_casedrop_reliability(d, method = "pcor", iter = 2,
                                  drop_prop = 0.2)
  expect_true(is.nan(split$mean[split$metric == "correlation"]))
  expect_true(all(is.finite(split$mean[split$metric != "correlation"])))
  expect_true(is.nan(drop$mean[drop$metric == "correlation"]))
  expect_true(all(is.finite(drop$mean[drop$metric != "correlation"])))
  expect_identical(attr(drop, "cs"), 0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(split), split)
  expect_identical(plot(drop), drop)
})

test_that("observed nodewise bootstrap fits retain original prediction columns", {
  set.seed(1)
  d <- data.frame(a = stats::rbinom(60, 1, 0.5),
    b = stats::rbinom(60, 1, 0.5), c = stats::rbinom(60, 1, 0.5))
  fit <- net_boot(d, method = "ising", n_boot = 2, cores = 1,
    labels = c("A", "B", "C"), nlambda = 10)
  expect_identical(fit$observed$nodewise$source_names, names(d))
  expect_s3_class(net_predict(fit$observed, d), "data.frame")
  expect_equal(net_predict(fit$observed, d)$predictability,
                 fit$observed$nodes$predictability)
})

test_that("undefined case-drop correlations count against the CS certainty", {
  skip_on_cran()
  set.seed(1)
  x <- stats::rnorm(100)
  d <- data.frame(x = x, y = 0.8 * x + stats::rnorm(100, sd = 0.6),
                   z = stats::rnorm(100))
  set.seed(3)
  cd <- net_casedrop_reliability(d, drop_prop = 0.9, iter = 40)
  corr <- attr(cd, "correlations")
  expect_true(anyNA(corr))
  expect_identical(attr(cd, "n_success"), 40L)
  expect_lt(sum(corr >= 0.7, na.rm = TRUE) / 40, 0.95)
  expect_identical(attr(cd, "cs"), 0)
  set.seed(3)
  st <- net_stability(d, drop_prop = 0.9, iter = 40)
  expect_true(all(st$table$prop_above < 0.95))
  expect_true(all(st$cs == 0))
  expect_identical(st$n_success, 40L)
})

test_that("composite groups and actor identifiers have collision-free identities", {
  d <- data.frame(g1 = rep(c("a-b", "a", "other"), each = 20),
    g2 = rep(c("c", "b-c", "z"), each = 20), framework_data())
  g <- psychnet(d, group = c("g1", "g2"), method = "pcor")
  expect_length(g, 3L)
  expect_false(anyDuplicated(names(g)) > 0L)
  expect_identical(unname(vapply(attr(g, "subsets"), nrow, integer(1))),
                   rep(20L, 3))
  ev <- data.frame(a = c("a.b", "a", "a.b", "a"),
    b = c("c", "b.c", "c", "b.c"), Action = c("x", "x", "y", "y"))
  f <- event_frequencies(ev, actor = c("a", "b"))
  expect_identical(nrow(f), 2L)
  expect_identical(length(unique(f$actor)), 2L)
  expect_true(all(f$x == 1 & f$y == 1))
})

test_that("grouped resampling forwards estimator dots and explicit overrides", {
  d <- data.frame(group = rep(c("A", "B"), each = 30), framework_data())
  g <- psychnet(d, group = "group", method = "cor")
  sub <- attr(g, "subsets")$A
  set.seed(9)
  grouped <- net_boot(g, n_boot = 3, cores = 1, cor_method = "spearman")
  set.seed(9)
  direct <- net_boot(sub, method = "cor", n_boot = 3, cores = 1,
                     cor_method = "spearman")
  expect_equal(grouped$A$observed$weights, direct$observed$weights)
  expect_equal(grouped$A$edge_boot, direct$edge_boot)
  set.seed(9)
  override <- net_boot(g, n_boot = 3, cores = 1, cor_method = "spearman",
    estimator_args = list(cor_method = "pearson"))
  expect_equal(override$A$observed$weights, g$A$weights)
  set.seed(9)
  sg <- net_split_reliability(g, iter = 2,
    estimator_args = list(cor_method = "spearman"))
  set.seed(9)
  sd <- net_split_reliability(sub, method = "cor", iter = 2,
                              estimator_args = list(cor_method = "spearman"))
  expect_equal(attr(sg$A, "iterations"), attr(sd, "iterations"))
})

test_that("dependent event groups are excluded from independent-row resampling", {
  d <- data.frame(id = rep(1:40, each = 3),
    group = rep(c("A", "B"), each = 60), framework_data(120))
  g <- psychnet(d, id = "id", group = "group", method = "pcor")
  expect_identical(g$A$n, 40L)
  expect_identical(attr(g, "source"), "eventdata")
  for (verb in list(net_boot, net_stability, net_casedrop_reliability,
                    net_split_reliability)) {
    expect_error(verb(g), "cross-sectional data only")
  }
  expect_error(net_compare(g, iter = 2), "cross-sectional data only")
  # Reducing to one independent occasion per person remains supported.
  d$id <- seq_len(nrow(d))
  independent <- psychnet(d, id = "id", group = "group", method = "pcor")
  expect_identical(attr(independent, "source"), "data")
  expect_s3_class(net_boot(independent, n_boot = 2, cores = 1),
                   "psychnet_bootstrap_group")
})

test_that("event feature missingness follows pairwise and listwise handling", {
  d <- data.frame(id = 1:60, framework_data())
  d$y[1] <- NA_real_
  fit <- psychnet(d, id = "id", method = "pcor", na_method = "pairwise")
  direct <- pcor_network(d[c("x", "y", "z")], na_method = "pairwise")
  expect_identical(fit$nodes$label, c("x", "y", "z"))
  expect_equal(fit$weights, direct$weights)
  expect_equal(fit$n, direct$n)
  # Missing observations do not poison the actor's other within observations.
  d$id <- rep(1:20, each = 3)
  dec <- psychnets:::.within_between_decompose(as.matrix(d[-1]), d$id)
  expect_true(is.na(dec$within[1, "y"]))
  expect_true(all(is.finite(dec$within[2:3, "y"])))
  expect_equal(unname(dec$between[1, "y"]), mean(d$y[2:3]))
  expect_s3_class(psychnet(d, id = "id", method = "pcor"), "psychnet")
  lw <- psychnet(d, id = "id", method = "pcor", na_method = "listwise")
  cc <- psychnet(d[stats::complete.cases(d), ], id = "id", method = "pcor")
  expect_equal(lw$weights, cc$weights)
  expect_equal(lw$n, cc$n)
})

test_that("groups of within/between networks print and summarize both levels", {
  d <- data.frame(id = rep(1:40, each = 3),
    group = rep(c("A", "B"), each = 60), framework_data(120))
  g <- psychnet(d, id = "id", group = "group", standardize = FALSE,
                method = "pcor")
  expect_output(print(g), "within")
  tab <- summary(g)
  expect_identical(nrow(tab), 4L)
  expect_identical(tab$group, c("A", "A", "B", "B"))
  expect_identical(tab$level, c("within", "between", "within", "between"))
  expect_identical(tab$nodes, rep(3L, 4))
})
