# Review 2026-11-07, resampling verbs (net_boot, difference_test,
# net_stability, net_casedrop_reliability, net_split_reliability).
# Finding ids R1-R10 refer to the review.

rs_data <- function(n = 120L, p = 4L, seed = 1L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p) %*%
    chol(0.4^abs(outer(seq_len(p), seq_len(p), "-")))
  colnames(x) <- paste0("V", seq_len(p))
  x
}

# A centrality_fn that errors whenever the V1-V2 weight exceeds `cut`: a
# failure OUTSIDE the estimator, deterministic per draw.
failing_centrality <- function(cut) {
  function(g) {
    if (g[1L, 2L] > cut) stop("boom: V1-V2 above the cut")
    list(custom = rowSums(g))
  }
}

# Record every warning class/message, muffled, while evaluating `expr`.
capture_warnings <- function(expr) {
  seen <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen[[length(seen) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

# ---- R1: a failing draw is caught per draw, identically for any core count --

test_that("R1 net_boot failures are per draw and identical for cores 1 and 2", {
  skip_on_cran()
  skip_on_os("windows")
  x <- rs_data()
  cut <- pcor_network(x)$weights[1L, 2L] + 1e-12
  run <- function(cores) {
    set.seed(42)
    capture_warnings(net_boot(x, method = "pcor", n_boot = 20, cores = cores,
                              measures = c("strength", "custom"),
                              centrality_fn = failing_centrality(cut)))
  }
  serial <- run(1L)
  parallel_run <- run(2L)
  s <- serial$value; pr <- parallel_run$value
  expect_gt(s$n_failed, 0L)
  expect_lt(s$n_failed, 20L)
  expect_identical(pr$n_failed, s$n_failed)
  expect_identical(pr$edge_boot, s$edge_boot)
  expect_identical(pr$centrality_boot, s$centrality_boot)
  expect_identical(pr$edges, s$edges)
  for (w in list(serial$warnings[[1L]], parallel_run$warnings[[1L]])) {
    expect_s3_class(w, "psychnet_resample_failed")
    expect_match(conditionMessage(w), "boom", fixed = TRUE)
  }
  # failed draws are NA rows, not dropped rows
  expect_identical(nrow(s$edge_boot), 20L)
  expect_identical(sum(stats::complete.cases(s$edge_boot)), 20L - s$n_failed)
})

test_that("R1 an observed centrality error is raised, not counted as a draw", {
  x <- rs_data()
  expect_error(net_boot(x, method = "pcor", n_boot = 5, cores = 1,
                        measures = "custom",
                        centrality_fn = function(g) stop("boom always")),
               "boom always")
})

test_that("R1 net_boot with every draw failing reports psychnet_resample_failed", {
  x <- rs_data()
  obs_w <- pcor_network(x)$weights[1L, 2L]
  # passes on the observed network only (its exact V1-V2 weight)
  only_observed <- function(g) {
    if (!isTRUE(all.equal(g[1L, 2L], obs_w, tolerance = 1e-12)))
      stop("boom draw")
    list(custom = rowSums(g))
  }
  set.seed(3)
  err <- tryCatch(net_boot(x, method = "pcor", n_boot = 4, cores = 1,
                           measures = "custom", centrality_fn = only_observed),
                  psychnet_resample_failed = function(e) e)
  expect_s3_class(err, "psychnet_resample_failed")
  expect_match(conditionMessage(err), "boom draw", fixed = TRUE)
})

# ---- R2: observed / full-sample failures carry the estimator's message -----

test_that("R2 observed and full-sample fit errors propagate the estimator message", {
  x <- rs_data()
  expect_error(net_boot(x, method = "pcor", n_boot = 2, cores = 1, foo = 1),
               class = "psychnet_estimation_failed")
  expect_error(net_boot(x, method = "pcor", n_boot = 2, cores = 1, foo = 1),
               "unused argument", fixed = TRUE)
  # net_stability has no `cores`: the estimator's complaint reaches the user
  expect_error(net_stability(x, method = "pcor", iter = 2, cores = 1),
               class = "psychnet_estimation_failed")
  expect_error(net_stability(x, method = "pcor", iter = 2, cores = 1),
               "cores", fixed = TRUE)
  expect_error(net_casedrop_reliability(x, method = "pcor", iter = 2,
                                        drop_prop = 0.5, foo = 1),
               class = "psychnet_estimation_failed")
  expect_error(net_split_reliability(x, method = "pcor", iter = 2, foo = 1),
               "unused argument", fixed = TRUE)
  expect_error(net_boot(x, n_boot = 2, cores = 1, gamma = "a"),
               class = "psychnet_estimation_failed")
  # an estimator's own psychnet class survives the re-raise
  e <- tryCatch(net_boot(x, method = "pcor", n_boot = 2, cores = 1,
                         labels = c("a", "a", "b", "c")),
                error = function(e) e)
  expect_s3_class(e, "psychnet_estimation_failed")
  expect_s3_class(e, "psychnet_bad_input")
})

test_that("R2 resample failures in net_stability keep the first error message", {
  x <- rs_data()
  cut <- pcor_network(x)$weights[1L, 2L] + 1e-12
  set.seed(5)
  res <- capture_warnings(net_stability(x, method = "pcor",
    measures = c("strength", "custom"),
    centrality_fn = failing_centrality(cut), drop_prop = c(0.3, 0.5),
    iter = 10))
  classes <- vapply(res$warnings, function(w) class(w)[1L], character(1))
  expect_true("psychnet_resample_failed" %in% classes)
  w <- res$warnings[[match("psychnet_resample_failed", classes)]]
  expect_match(conditionMessage(w), "boom", fixed = TRUE)
  st <- res$value
  expect_gt(st$n_failed, 0L)
  # a failed subset counts as not stable
  cors <- as.data.frame(st, what = "correlations")
  expect_identical(sum(is.na(cors$correlation[cors$measure == "strength"])),
                   st$n_failed)
})

# ---- R3: split-half reliability reports undefined splits -------------------

test_that("R3 an empty half against a non-empty one counts as 0, not undefined", {
  # 3 of these 10 splits leave one half-network empty: that is a disagreement
  # (cor0 = 0, as in bootnet), not a split to drop -- dropping it inflated the
  # mean correlation.
  set.seed(1)
  x <- matrix(stats::rnorm(40 * 3), 40, 3)
  colnames(x) <- c("a", "b", "c")
  x[, 2] <- x[, 2] + 0.7 * x[, 1]
  set.seed(10)
  expect_no_warning(r <- net_split_reliability(x, method = "glasso", iter = 10))
  cor_row <- r[r$metric == "correlation", ]
  draws <- attr(r, "iterations")$correlation
  expect_identical(cor_row$n_valid, 10L)
  expect_identical(sum(draws == 0), 3L)
  expect_equal(cor_row$mean, mean(draws), tolerance = 1e-12)
})

test_that("R3 split reliability counts splits where both halves are constant", {
  # independent items: both half-networks are often empty -> the correlation is
  # undefined there (both agree, but no correlation exists); counted and warned
  set.seed(3)
  x <- matrix(stats::rnorm(40 * 3), 40, 3)
  colnames(x) <- c("a", "b", "c")
  set.seed(4)
  expect_warning(r <- net_split_reliability(x, method = "glasso", iter = 10),
                 class = "psychnet_undefined_statistic")
  cor_row <- r[r$metric == "correlation", ]
  draws <- attr(r, "iterations")$correlation
  expect_lt(cor_row$n_valid, 10L)
  expect_identical(cor_row$n_valid, sum(is.finite(draws)))
  expect_true(all(r$n_valid[r$metric != "correlation"] == 10L))
})

test_that("R3 a metric defined in no split is NA, not NaN", {
  d <- rs_data(p = 2L)
  set.seed(2)
  expect_warning(r <- net_split_reliability(d, method = "pcor", iter = 4),
                 class = "psychnet_undefined_statistic")
  cor_row <- r[r$metric == "correlation", ]
  expect_true(is.na(cor_row$mean))
  expect_false(is.nan(cor_row$mean))
  expect_identical(cor_row$n_valid, 0L)
  # no undefined metric -> no warning
  set.seed(2)
  expect_no_warning(net_split_reliability(rs_data(), method = "pcor", iter = 3))
})

# ---- R4: the bootnet CS rule ------------------------------------------------

test_that("R4 .psn_cor0 follows bootnet's cor0", {
  cor0 <- psychnets:::.psn_cor0
  expect_identical(cor0(c(1, 1, 1), c(1, 2, 3)), 0)
  expect_identical(cor0(c(1, 2, 3), c(5, 5, 5)), 0)
  expect_identical(cor0(1, 2), NA_real_)
  expect_identical(cor0(c(1, NA), c(2, 3)), NA_real_)
  expect_identical(cor0(c(1, Inf, 2), c(1, 2, 3)), NA_real_)
  expect_equal(cor0(c(1, 2, 4), c(2, 1, 5)), stats::cor(c(1, 2, 4), c(2, 1, 5)),
               tolerance = 1e-15)
  # split-half mode: undefined only when BOTH sides are constant; an empty
  # half against a non-empty one is a disagreement (0), as in bootnet's cor0
  expect_identical(cor0(c(1, 1, 1), c(2, 2, 2), constant = "na"), NA_real_)
  expect_identical(cor0(c(0, 0, 0), c(1, 2, 3), constant = "na"), 0)
  # float-noise constancy counts as constant
  expect_identical(cor0(rep(0.1, 7) + c(0, 1e-17, 0, 0, 0, 0, 0), 1:7), 0)
})

test_that("R4 CS rule is strict on threshold and certainty (hand-computed)", {
  cs_rule <- psychnets:::.psn_cs_coefficient
  # 20 subsets x 3 proportions
  corr <- matrix(0.9, 20, 3)
  corr[1, 2] <- 0.7               # exactly at threshold: NOT stable
  corr[1:2, 3] <- 0.5             # 18/20 = 0.90 stable
  expect_identical(cs_rule(corr, c(0.1, 0.2, 0.3), 0.7, 0.95), 0.1)
  # 19/20 = 0.95 exactly is not > 0.95
  corr2 <- matrix(0.9, 20, 1); corr2[1, 1] <- 0.1
  expect_identical(cs_rule(corr2, 0.4, 0.7, 0.95), 0)
  expect_identical(cs_rule(corr2, 0.4, 0.7, 0.9), 0.4)
  # an undefined correlation fails (bootnet would drop it with na.rm)
  corr3 <- matrix(0.9, 20, 1); corr3[1, 1] <- NA
  expect_identical(cs_rule(corr3, 0.4, 0.7, 0.95), 0)
  # order of drop_prop does not matter
  expect_identical(cs_rule(corr[, 3:1], c(0.3, 0.2, 0.1), 0.7, 0.95), 0.1)
})

test_that("R4 net_stability CS equals the strict rule on its stored subsets", {
  skip_on_cran()
  x <- rs_data(n = 100L, p = 5L)
  set.seed(7)
  st <- net_stability(x, method = "pcor", drop_prop = c(0.2, 0.4, 0.6, 0.8),
                      iter = 20)
  expect_identical(st$cor_method, "pearson")
  cors <- as.data.frame(st, what = "correlations")
  expect_identical(nrow(cors), 2L * 4L * 20L)
  for (m in names(st$cs)) {
    by_prop <- split(cors$correlation[cors$measure == m],
                     cors$drop_prop[cors$measure == m])
    share <- vapply(by_prop, function(v) mean(!is.na(v) & v > 0.7), numeric(1))
    hand <- if (any(share > 0.95)) max(as.numeric(names(share))[share > 0.95]) else 0
    expect_identical(st$cs[[m]], hand)
    tab <- as.data.frame(st)
    expect_equal(tab$prop_above[tab$measure == m], unname(share),
                 tolerance = 1e-15)
  }
  # the subsets are Pearson correlations of the subset vs full centralities
  set.seed(7)
  idx <- psychnets:::.psn_casedrop_indices(100L,
    psychnets:::.psn_keep_sizes(c(0.2, 0.4, 0.6, 0.8), 100L), 20L)
  full <- net_centralities(pcor_network(x))
  first <- net_centralities(pcor_network(x[idx[[1L]], ]))
  expect_equal(cors$correlation[cors$measure == "strength"][1L],
               stats::cor(full$strength, first$strength), tolerance = 1e-12)
  # spearman stays available
  set.seed(7)
  sp <- net_stability(x, method = "pcor", drop_prop = 0.2, iter = 5,
                      cor_method = "spearman")
  expect_identical(sp$cor_method, "spearman")
})


test_that("R4 casedrop CS uses pearson and the strict rule", {
  skip_on_cran()
  x <- rs_data(n = 100L, p = 5L)
  set.seed(8)
  cd <- net_casedrop_reliability(x, method = "pcor",
                                 drop_prop = c(0.2, 0.5, 0.8), iter = 20)
  expect_identical(attr(cd, "cor_method"), "pearson")
  corr <- attr(cd, "correlations")
  share <- colMeans(!is.na(corr) & corr > 0.7)
  hand <- if (any(share > 0.95)) max(c(0.2, 0.5, 0.8)[share > 0.95]) else 0
  expect_identical(attr(cd, "cs"), hand)
  expect_named(cd, c("metric", "drop_prop", "actual_drop", "mean", "sd",
                     "n_valid"))
})

test_that("R4 threshold/certainty of 1 are rejected (strict rule cannot pass)", {
  x <- rs_data()
  expect_error(net_stability(x, method = "pcor", certainty = 1),
               class = "psychnet_bad_input")
  expect_error(net_stability(x, method = "pcor", threshold = 1),
               class = "psychnet_bad_input")
  expect_error(net_casedrop_reliability(x, method = "pcor", certainty = 1),
               class = "psychnet_bad_input")
  expect_error(net_stability(x, method = "pcor", cor_method = "kendall"))
})

# ---- R5: tidy accessors -------------------------------------------------------

test_that("R5 as.data.frame.psychnet_stability returns the table and the CS", {
  x <- rs_data()
  set.seed(1)
  st <- net_stability(x, method = "pcor", drop_prop = c(0.3, 0.6), iter = 4)
  tab <- as.data.frame(st)
  expect_identical(tab, st$table)
  expect_named(tab, c("measure", "drop_prop", "actual_drop", "mean_cor",
                      "sd_cor", "prop_above", "n_valid"))
  cs <- as.data.frame(st, what = "cs")
  expect_identical(cs, data.frame(measure = names(st$cs), cs = unname(st$cs),
                                  stringsAsFactors = FALSE))
  expect_error(as.data.frame(st, what = "nope"))
})

test_that("R5 as.data.frame.psychnet_bootstrap serves every table via `what`", {
  x <- rs_data()
  set.seed(2)
  bs <- net_boot(x, method = "glasso", n_boot = 6, cores = 1,
                 threshold = TRUE, predictability = TRUE)
  expect_identical(as.data.frame(bs), bs$edges)
  cent <- as.data.frame(bs, what = "centrality")
  expect_identical(nrow(cent), 4L)
  expect_named(cent, c("node", "strength", "strength_lower", "strength_upper",
                       "expected_influence", "expected_influence_lower",
                       "expected_influence_upper"))
  pred <- as.data.frame(bs, what = "predictability")
  expect_named(pred, c("node", "value", "lower", "upper"))
  th <- as.data.frame(bs, what = "thresholded")
  expect_named(th, c("from", "to", "observed", "weight", "retained"))
  expect_identical(th$retained, bs$edges$significant)
  expect_true(all(th$weight[!th$retained] == 0))
  expect_equal(th$weight[th$retained], th$observed[th$retained],
               tolerance = 0)
  sig <- as.data.frame(bs, significant = TRUE)
  expect_true(all(sig$significant))
  # not requested -> classed error
  set.seed(2)
  plain <- net_boot(x, method = "pcor", n_boot = 3, cores = 1)
  expect_error(as.data.frame(plain, what = "predictability"),
               class = "psychnet_bad_input")
  expect_error(as.data.frame(plain, what = "thresholded"),
               class = "psychnet_bad_input")
  expect_error(as.data.frame(plain, what = "centrality", significant = TRUE),
               class = "psychnet_bad_input")
})

# ---- R6: ci validation --------------------------------------------------------

test_that("R6 ci must be a single number in (0, 1)", {
  x <- rs_data()
  expect_error(net_boot(x, method = "pcor", n_boot = 3, cores = 1,
                        ci = c(0.9, 0.95)), class = "psychnet_bad_input")
  expect_error(net_boot(x, method = "pcor", n_boot = 3, cores = 1, ci = 1),
               class = "psychnet_bad_input")
  set.seed(1)
  bs <- net_boot(x, method = "pcor", n_boot = 5, cores = 1)
  expect_error(difference_test(bs, ci = c(0.9, 0.95)),
               class = "psychnet_bad_input")
  expect_error(difference_test(bs, ci = "0.9"), class = "psychnet_bad_input")
  expect_s3_class(difference_test(bs, ci = 0.9), "psychnet_difference")
})

# ---- R7 / R8: non-negative directed estimators --------------------------------

test_that("R7 relimp threshold / significance raise psychnet_uninformative_test", {
  skip_on_cran()
  x <- rs_data(p = 3L)
  set.seed(4)
  expect_warning(bs <- net_boot(x, method = "relimp", n_boot = 5, cores = 1,
                                threshold = TRUE),
                 class = "psychnet_uninformative_test")
  expect_warning(as.data.frame(bs, significant = TRUE),
                 class = "psychnet_uninformative_test")
  # glasso is not affected
  set.seed(4)
  expect_no_warning(g <- net_boot(x, method = "pcor", n_boot = 3, cores = 1,
                                  threshold = TRUE))
  expect_no_warning(as.data.frame(g, significant = TRUE))
})

test_that("R8 directed edges are labelled A->B in net_boot and difference_test", {
  skip_on_cran()
  x <- rs_data(p = 3L)
  set.seed(4)
  bs <- net_boot(x, method = "relimp", n_boot = 4, cores = 1)
  expect_true(all(grepl("->", bs$edge_labels, fixed = TRUE)))
  expect_false(anyDuplicated(bs$edge_labels) > 0L)
  dt <- difference_test(bs, type = "edge")
  expect_true(all(grepl("->", c(dt$item1, dt$item2), fixed = TRUE)))
  set.seed(4)
  und <- net_boot(x, method = "pcor", n_boot = 3, cores = 1)
  expect_true(all(grepl("--", und$edge_labels, fixed = TRUE)))
})

# ---- R10: explicit checks -----------------------------------------------------

test_that("R10 cores and difference_test type are validated with classes", {
  rc <- psychnets:::.resolve_cores
  expect_error(rc("two"), class = "psychnet_bad_input")
  expect_error(rc(0), class = "psychnet_bad_input")
  expect_error(rc(1.5), class = "psychnet_bad_input")
  expect_error(rc(c(1, 2)), class = "psychnet_bad_input")
  expect_identical(rc(1), 1L)
  x <- rs_data()
  set.seed(1)
  bs <- net_boot(x, method = "pcor", n_boot = 3, cores = 1)
  expect_error(difference_test(bs, type = c("edge", "strength")),
               class = "psychnet_bad_input")
  expect_error(difference_test(bs, type = "betweenness"),
               class = "psychnet_bad_input")
  expect_error(difference_test(list(), type = "edge"),
               class = "psychnet_bad_input")
  expect_error(net_boot(x, method = "pcor", n_boot = 0, cores = 1),
               class = "psychnet_bad_input")
  expect_error(net_stability(x, method = "pcor", drop_prop = c(0.2, NA)),
               class = "psychnet_bad_input")
  expect_error(net_stability(x, method = "pcor", iter = "5"),
               class = "psychnet_bad_input")
})

test_that("R10 a constant full-sample centrality warns and gives CS 0", {
  d <- rs_data(p = 2L)
  set.seed(3)
  res <- capture_warnings(net_stability(d, method = "pcor",
                                        drop_prop = c(0.2, 0.4), iter = 5))
  classes <- vapply(res$warnings, function(w) class(w)[1L], character(1))
  expect_identical(sum(classes == "psychnet_constant_centrality"), 1L)
  expect_true(all(res$value$cs == 0))
  cors <- as.data.frame(res$value, what = "correlations")
  expect_true(all(cors$correlation == 0))
})

test_that("actual_drop reports the proportion dropped after rounding", {
  x <- rs_data(n = 30L)
  set.seed(1)
  st <- net_stability(x, method = "pcor", drop_prop = 0.25, iter = 2)
  tab <- as.data.frame(st)
  expect_equal(unique(tab$actual_drop), 1 - 22 / 30, tolerance = 1e-15)
  expect_identical(unique(tab$drop_prop), 0.25)
})

test_that("pre-drawn case-drop subsets reproduce the former nested-loop order", {
  set.seed(99)
  idx <- psychnets:::.psn_casedrop_indices(50L, c(40L, 25L), 3L)
  set.seed(99)
  loop <- unlist(lapply(c(40L, 25L), function(k)
    lapply(seq_len(3L), function(i) sample.int(50L, k, replace = FALSE))),
    recursive = FALSE)
  expect_identical(idx, loop)
})

test_that("R3/R4 casedrop warns when the edge correlation is undefined", {
  d <- rs_data(p = 2L)
  set.seed(1)
  expect_warning(cd <- net_casedrop_reliability(d, method = "pcor", iter = 3,
                                                drop_prop = 0.3),
                 class = "psychnet_undefined_statistic")
  expect_identical(attr(cd, "cs"), 0)
  expect_identical(cd$n_valid[cd$metric == "correlation"], 0L)
  # an empty full-sample network: every correlation is 0 (cor0), CS 0
  set.seed(2)
  z <- matrix(stats::rnorm(60 * 3), 60, 3)
  expect_warning(cz <- net_casedrop_reliability(z, method = "glasso", iter = 3,
                                                drop_prop = 0.3),
                 class = "psychnet_constant_centrality")
  expect_identical(attr(cz, "cs"), 0)
})

test_that("column-drop notices are raised once by the observed fit, not per draw", {
  set.seed(5)
  d <- data.frame(g1 = stats::rnorm(80), g2 = stats::rnorm(80),
                  g3 = stats::rnorm(80),
                  f = factor(sample(c("a", "b", "c"), 80, TRUE)))
  n_msg <- 0L
  withCallingHandlers(
    bs <- net_boot(d, method = "pcor", n_boot = 4, cores = 1),
    psychnet_dropped_columns = function(m) {
      n_msg <<- n_msg + 1L
      invokeRestart("muffleMessage")
    })
  expect_identical(n_msg, 1L)
  expect_identical(bs$node_labels, c("g1", "g2", "g3"))
})
