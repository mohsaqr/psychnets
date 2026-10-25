# net_compare(): the network comparison test. Generators/configs live in
# helper-nct.R; fixtures/make-nct.R freezes the baseline and NCT reference.

nct_fit <- function(i, ...) {
  cfg <- nct_configs[i, ]
  g <- nct_groups(cfg$seed)
  set.seed(cfg$seed)
  net_compare(g$a, g$b, iter = nct_iter, paired = cfg$paired, abs = cfg$abs,
              weighted = cfg$weighted, ...)
}

test_that("complete-data output is frozen (input handling moves no number)", {
  skip_slow()
  base <- readRDS(test_path("fixtures", "nct-baseline.rds"))
  invisible(lapply(seq_len(nrow(nct_configs)), function(i)
    expect_identical(nct_numbers(nct_fit(i)), base[[i]], label = sprintf("config %d", i))))
})

test_that("permutation p-values equal NetworkComparisonTest::NCT under the same seed", {
  skip_slow()
  ref <- readRDS(test_path("fixtures", "nct-reference.rds"))
  invisible(lapply(seq_len(nrow(nct_configs)), function(i) {
    fit <- nct_fit(i)
    lbl <- sprintf("config %d", i)
    # same seed -> same permuted rows; the only difference is the glasso
    # tolerance (qgraph thr = 1e-4), so p-values agree exactly and the
    # permutation statistics to independent-solver precision
    expect_equal(fit$strength$p_value, ref[[i]]$strength_p, label = lbl)
    expect_equal(fit$structure$p_value, ref[[i]]$structure_p, label = lbl)
    expect_equal(fit$E$p_value, ref[[i]]$E_p, label = lbl)
    expect_lt(max(abs(fit$strength$perm - ref[[i]]$strength_perm)), 1e-3)
    expect_lt(max(abs(fit$structure$perm - ref[[i]]$structure_perm)), 1e-3)
    expect_lt(max(abs(fit$strength$by_group - ref[[i]]$strength_by_group)), 1e-3)
    expect_lt(max(abs(unname(fit$nw1) - ref[[i]]$nw1)), 1e-3)
  }))
})

test_that("a missing cell no longer zeroes the variable's edges (pairwise)", {
  skip_slow()
  g <- nct_groups(1)
  a_na <- g$a; a_na[3, 2] <- NA
  set.seed(1); clean <- net_compare(g$a, g$b, iter = 5)
  set.seed(1); fit <- net_compare(a_na, g$b, iter = 5)
  expect_gt(sum(abs(fit$nw1[2, ])), 0)
  expect_lt(max(abs(fit$nw1 - clean$nw1)), 0.05)
  expect_identical(fit$params$na_method, "pairwise")
})

test_that("listwise equals running on the complete rows", {
  skip_slow()
  g <- nct_groups(2)
  a_na <- g$a; a_na[c(3, 9), 2] <- NA
  set.seed(4); lw <- net_compare(a_na, g$b, iter = 5, na_method = "listwise")
  set.seed(4); cc <- net_compare(a_na[stats::complete.cases(a_na), ], g$b, iter = 5)
  expect_identical(nct_numbers(lw), nct_numbers(cc))
  expect_identical(lw$n, c(nrow(g$a) - 2L, nrow(g$b)))
})

test_that("paired listwise drops the whole pair", {
  skip_slow()
  g <- nct_groups(3)
  a_na <- g$a; a_na[5, 1] <- NA
  set.seed(2)
  fit <- net_compare(a_na, g$b, iter = 3, paired = TRUE, na_method = "listwise")
  expect_identical(fit$n, rep(nrow(g$a) - 1L, 2L))
})

test_that("a constant variable is an isolated node, with a warning when observed", {
  skip_slow()
  g <- nct_groups(1)
  g$a[, 5] <- 0
  set.seed(1)
  expect_warning(fit <- net_compare(g$a, g$b, iter = 5),
                 class = "psychnet_constant_variable")
  expect_identical(unname(fit$nw1[5, ]), rep(0, 5))
  expect_gt(sum(abs(fit$nw2[5, ])), 0)
})

test_that("a variable constant in some permuted splits keeps the test defined", {
  skip_slow()
  g <- nct_groups(1)
  g$a[, 5] <- 0; g$a[1:2, 5] <- 1
  g$b[, 5] <- 0; g$b[1, 5] <- 1
  set.seed(1)
  expect_no_warning(fit <- net_compare(g$a, g$b, iter = 30))
  expect_true(all(is.finite(c(fit$strength$perm, fit$structure$perm, fit$E$perm))))
})

test_that("a never co-observed pair is refused, not zero-filled", {
  g <- nct_groups(1)
  g$a[1:60, 1] <- NA
  g$a[61:120, 2] <- NA
  expect_error(net_compare(g$a, g$b, iter = 3), class = "psychnet_degenerate_sample")
})

test_that("columns must match in name and order", {
  g <- nct_groups(1)
  expect_error(net_compare(g$a, g$b[, c(3, 1, 5, 2, 4)], iter = 3),
               class = "psychnet_bad_input")
  b2 <- g$b; colnames(b2)[5] <- "W5"
  expect_error(net_compare(g$a, b2, iter = 3), class = "psychnet_bad_input")
  expect_error(net_compare(g$a, g$b[, 1:4], iter = 3), class = "psychnet_bad_input")
})

test_that("non-numeric columns are dropped", {
  skip_slow()
  g <- nct_groups(1)
  da <- data.frame(id = seq_len(nrow(g$a)) |> as.character(), g$a)
  set.seed(1); with_id <- net_compare(da, g$b, iter = 3)
  set.seed(1); plain <- net_compare(g$a, g$b, iter = 3)
  expect_identical(nct_numbers(with_id), nct_numbers(plain))
})

test_that("malformed arguments raise psychnet_bad_input", {
  g <- nct_groups(1)
  bad_iter <- list(0, NA, 2.7, "10", c(5, 6))
  invisible(lapply(bad_iter, function(it)
    expect_error(net_compare(g$a, g$b, iter = it), class = "psychnet_bad_input")))
  expect_error(net_compare(g$a, g$b, paired = NA), class = "psychnet_bad_input")
  expect_error(net_compare(g$a, g$b, gamma = -1), class = "psychnet_bad_input")
  expect_error(net_compare(g$a), class = "psychnet_bad_input")
  expect_error(net_compare(g$a, g$b[-1, ], paired = TRUE), class = "psychnet_bad_input")
  expect_error(net_compare(g$a, g$b, na_method = "mean"))
})

test_that("group mode says when it is not comparing the group's own networks", {
  skip_slow()
  g <- nct_groups(1)
  d <- data.frame(grp = rep(c("a", "b"), each = nrow(g$a)), rbind(g$a, g$b))
  expect_warning(net_compare(psychnet(d, group = "grp", method = "cor"), iter = 3),
                 class = "psychnet_estimator_mismatch")
  expect_warning(net_compare(psychnet(d, group = "grp", threshold = 0.1), iter = 3),
                 class = "psychnet_estimator_mismatch")
  expect_warning(net_compare(psychnet(d, group = "grp", cor_method = "spearman"),
                             iter = 3), class = "psychnet_estimator_mismatch")
})

test_that("group mode equals the direct call", {
  skip_slow()
  g <- nct_groups(1)
  d <- data.frame(grp = rep(c("a", "b"), each = nrow(g$a)), rbind(g$a, g$b))
  set.seed(1); via_group <- net_compare(psychnet(d, group = "grp"), iter = 5)
  set.seed(1); direct <- net_compare(g$a, g$b, iter = 5)
  expect_identical(nct_numbers(via_group), nct_numbers(direct))
})

test_that("summary() gives the tests, as.data.frame() the edges", {
  skip_slow()
  fit <- nct_fit(1)
  ed <- as.data.frame(fit)
  expect_named(ed, c("from", "to", "weight_1", "weight_2", "difference", "p_value",
                     "significant"))
  expect_identical(nrow(ed), 10L)                   # 5 nodes -> 10 pairs
  expect_equal(ed$difference, abs(ed$weight_1 - ed$weight_2))
  expect_true(all(ed$p_value >= 1 / (nct_iter + 1) & ed$p_value <= 1))
  expect_identical(ed$significant, ed$p_value < 0.05)
  gl <- summary(fit)
  expect_identical(gl$statistic, c("M", "S"))
  expect_identical(gl$p_value, c(fit$structure$p_value, fit$strength$p_value))
  expect_identical(gl$test, c("network structure", "global strength"))
  expect_named(gl, c("test", "statistic", "observed", "p_value", "significant"))
  expect_identical(as.data.frame(fit, alpha = 0.5)$significant, ed$p_value < 0.5)
})

test_that("swapping the groups leaves every statistic unchanged", {
  skip_slow()
  g <- nct_groups(2)
  set.seed(9); ab <- net_compare(g$a, g$b, iter = 5)
  set.seed(9); ba <- net_compare(g$b, g$a, iter = 5)
  # observed statistics are symmetric in the two groups
  expect_equal(ab$strength$observed, ba$strength$observed)
  expect_equal(ab$structure$observed, ba$structure$observed)
  expect_equal(ab$E$observed, ba$E$observed)
  expect_identical(ab$nw1, ba$nw2)
})

test_that("print shows the three tests and every edge", {
  skip_slow()
  expect_snapshot(print(nct_fit(1)))
})

test_that("M and S follow the NCT notation; the old fields fail loudly", {
  skip_slow()
  fit <- nct_fit(1)
  # M = largest edge difference (network structure), S = global strength
  expect_identical(fit$structure$observed, max(fit$E$observed))
  ut <- upper.tri(fit$nw1)
  expect_equal(fit$strength$by_group, c(sum(abs(fit$nw1[ut])), sum(abs(fit$nw2[ut]))))
  expect_equal(fit$strength$observed, abs(diff(fit$strength$by_group)))
  expect_error(fit$M, class = "psychnet_renamed")
  expect_error(fit$S, class = "psychnet_renamed")
  ei <- nct_fit(3)                                   # config 3 has abs = FALSE
  expect_identical(summary(ei)$test[2],
                   "global expected influence")
  expect_equal(ei$strength$by_group, c(sum(ei$nw1[ut]), sum(ei$nw2[ut])))
})

# ---- invariance / property tests (no reference package needed) ---------------

test_that("statistics are invariant to rescaling and shifting the variables", {
  skip_slow()
  g <- nct_groups(2)
  scale_shift <- function(x) sweep(sweep(x, 2, c(2, 0.5, 10, 3, 0.1), "*"), 2,
                                   c(5, -3, 100, 0, 1), "+")
  set.seed(5); raw <- net_compare(g$a, g$b, iter = 8)
  set.seed(5); tx <- net_compare(scale_shift(g$a), scale_shift(g$b), iter = 8)
  expect_equal(tx$nw1, raw$nw1, tolerance = 1e-8)
  expect_equal(tx$strength$perm, raw$strength$perm, tolerance = 1e-8)
  expect_equal(tx$structure$perm, raw$structure$perm, tolerance = 1e-8)
  expect_identical(tx$E$p_value, raw$E$p_value)
})

test_that("relabelling the variables permutes the result and nothing else", {
  skip_slow()
  g <- nct_groups(3)
  ord <- c(4, 2, 5, 1, 3)
  set.seed(6); a <- net_compare(g$a, g$b, iter = 8)
  set.seed(6); b <- net_compare(g$a[, ord], g$b[, ord], iter = 8)
  expect_equal(b$nw1, a$nw1[ord, ord], tolerance = 1e-8)
  expect_equal(b$strength$perm, a$strength$perm, tolerance = 1e-8)
  expect_equal(b$structure$perm, a$structure$perm, tolerance = 1e-8)
  key <- function(d) paste(pmin(d$from, d$to), pmax(d$from, d$to))
  ea <- as.data.frame(a); eb <- as.data.frame(b)
  expect_equal(eb$p_value[order(key(eb))], ea$p_value[order(key(ea))])
})

test_that("identical groups give zero statistics and p = 1", {
  skip_slow()
  g <- nct_groups(1)
  set.seed(2); fit <- net_compare(g$a, g$a, iter = 10, paired = TRUE)
  expect_identical(fit$structure$observed, 0)
  expect_identical(fit$strength$observed, 0)
  expect_true(all(fit$E$observed == 0))
  # every permuted split of identical pairs is the same pair of networks
  expect_true(all(c(fit$structure$p_value, fit$strength$p_value, fit$E$p_value) == 1))
})

test_that("p-values lie on the permutation lattice k / (iter + 1)", {
  skip_slow()
  fit <- nct_fit(1)
  p <- c(fit$structure$p_value, fit$strength$p_value, fit$E$p_value)
  k <- p * (nct_iter + 1)
  expect_true(all(abs(k - round(k)) < 1e-9))
  expect_true(all(k >= 1 & k <= nct_iter + 1))
})

test_that("p_adjust equals stats::p.adjust of the unadjusted edge p-values", {
  skip_slow()
  raw <- nct_fit(1)
  invisible(lapply(c("holm", "BH", "bonferroni", "BY"), function(m) {
    adj <- nct_fit(1, p_adjust = m)
    expect_equal(adj$E$p_value, stats::p.adjust(raw$E$p_value, method = m))
    # the global tests are never adjusted
    expect_identical(adj$structure$p_value, raw$structure$p_value)
    expect_identical(adj$strength$p_value, raw$strength$p_value)
  }))
})

test_that("weighted = FALSE compares binary networks", {
  skip_slow()
  fit <- nct_fit(4)                                  # config 4: weighted = FALSE
  ut <- upper.tri(fit$nw1)
  expect_true(all(c(fit$nw1, fit$nw2) %in% c(0, 1)))
  expect_equal(fit$strength$by_group, c(sum(fit$nw1[ut]), sum(fit$nw2[ut])))
  expect_true(all(fit$E$observed %in% c(0, 1)))
})

test_that("results are reproducible under a seed", {
  skip_slow()
  expect_identical(nct_numbers(nct_fit(2)), nct_numbers(nct_fit(2)))
})

test_that("print names the groups and lists the largest edge differences", {
  skip_slow()
  g <- nct_groups(1)
  group_a <- g$a; group_b <- g$b
  set.seed(1)
  fit <- net_compare(group_a, group_b, iter = 10)
  expect_identical(fit$groups, c("group_a", "group_b"))
  all_lines <- capture.output(print(fit))
  expect_true(any(grepl("Edges, largest difference first", all_lines, fixed = TRUE)))
  expect_identical(sum(grepl("^  V[0-9] - V[0-9]", all_lines)), 10L)   # every edge
  out <- capture.output(print(fit, n = 3))
  expect_match(out[1], "group_a vs group_b", fixed = TRUE)
  expect_true(any(grepl("Largest edge differences (3 of 10)", out, fixed = TRUE)))
  # the listed edges are the three largest differences, in order
  ed <- as.data.frame(fit)
  top <- ed[order(-ed$difference), ][1:3, ]
  listed <- out[grep("Edge", out) + 1:3]
  expect_true(all(startsWith(trimws(listed), paste(top$from, "-", top$to))))
  # group mode uses the level names
  d <- data.frame(grp = rep(c("early", "late"), each = nrow(g$a)), rbind(g$a, g$b))
  set.seed(1)
  expect_identical(net_compare(psychnet(d, group = "grp"), iter = 3)$groups,
                   c("early", "late"))
})

test_that("two empty networks give M = S = 0 and p = 1 (no floating residue)", {
  skip_slow()
  # p = 12, n = 60 independent data: both groups and nearly every permuted
  # split select the empty graph; a residue edge used to make M = 5.6e-17 and
  # turn "no difference" into p = 0.02.
  set.seed(30)
  a <- matrix(stats::rnorm(60 * 12), 60, 12); b <- matrix(stats::rnorm(60 * 12), 60, 12)
  colnames(a) <- colnames(b) <- paste0("V", 1:12)
  set.seed(1)
  fit <- net_compare(a, b, iter = 20)
  if (all(fit$nw1 == 0) && all(fit$nw2 == 0)) {
    expect_identical(fit$structure$observed, 0)
    expect_identical(fit$strength$observed, 0)
    expect_identical(c(fit$structure$p_value, fit$strength$p_value), c(1, 1))
  }
  expect_true(all(c(fit$nw1, fit$nw2) == 0 | abs(c(fit$nw1, fit$nw2)) > 1e-10))
})
