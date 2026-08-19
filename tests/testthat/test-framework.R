mk <- function(seed, n = 150, p = 5) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p) %*% chol(0.4^abs(outer(1:p, 1:p, "-")))
  colnames(X) <- paste0("V", seq_len(p))
  X
}

test_that("net_compare returns valid invariants and p-values", {
  skip_slow()
  fit <- net_compare(mk(1), mk(2), iter = 30)
  expect_s3_class(fit, "psychnet_nct")
  expect_true(fit$M$p_value >= 0 && fit$M$p_value <= 1)
  expect_true(fit$S$p_value >= 0 && fit$S$p_value <= 1)
  expect_equal(dim(fit$nw1), c(5L, 5L))
  expect_length(fit$M$perm, 30)
})

test_that("nearest-correlation projection returns a valid correlation matrix", {
  set.seed(9)
  S <- stats::cor(matrix(stats::rnorm(6 * 8), 6, 8))   # n < p -> not PD
  P <- .nearest_pd_cor(S)
  expect_true(all(abs(diag(P) - 1) < 1e-8))
  expect_gt(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("net_boot returns tidy edge and centrality CIs", {
  skip_slow()
  bs <- net_boot(mk(3), n_boot = 40, cores = 1)
  expect_s3_class(bs, "psychnet_bootstrap")
  expect_named(bs$edges,
               c("from", "to", "observed", "mean", "lower", "upper",
                 "prop_nonzero", "significant"))
  expect_true(all(bs$edges$lower <= bs$edges$upper))
  expect_true(all(bs$edges$prop_nonzero >= 0 & bs$edges$prop_nonzero <= 1))
  expect_equal(nrow(bs$edges), 10)              # 5 nodes -> 10 upper-tri edges
  expect_true(all(bs$centrality$strength_lower <= bs$centrality$strength_upper))
  expect_type(bs$edges$significant, "logical")
  # significance must agree with the stored interval
  expect_equal(bs$edges$significant, bs$edges$lower > 0 | bs$edges$upper < 0)
  # raw draws are retained for difference_test()
  expect_equal(dim(bs$edge_boot), c(40L, 10L))
  expect_identical(as.data.frame(bs), bs$edges)
})

test_that("parallel bootstrap is byte-identical to the serial run", {
  skip_slow()
  set.seed(123); bs1 <- net_boot(mk(3), n_boot = 40, cores = 1)
  set.seed(123); bs2 <- net_boot(mk(3), n_boot = 40, cores = 2)
  expect_equal(bs1$edge_boot, bs2$edge_boot)
  expect_equal(bs1$str_boot, bs2$str_boot)
  expect_equal(bs1$edges, bs2$edges)
})

test_that("difference_test returns a tidy pairwise table with sound intervals", {
  skip_slow()
  bs <- net_boot(mk(3), n_boot = 60, cores = 1)
  for (ty in c("edge", "strength", "expected_influence")) {
    dt <- difference_test(bs, type = ty)
    expect_named(dt, c("item1", "item2", "value1", "value2", "obs_diff",
                       "lower", "upper", "p_value", "significant"))
    expect_true(all(dt$lower <= dt$upper))
    expect_true(all(dt$p_value >= 0 & dt$p_value <= 1))
    expect_equal(dt$obs_diff, dt$value1 - dt$value2)
    expect_equal(dt$significant, dt$lower > 0 | dt$upper < 0)
  }
  # 5 edges of centrality -> choose(5,2) = 10 node pairs; 10 edges -> 45 pairs
  expect_equal(nrow(difference_test(bs, type = "strength")), 10)
  expect_equal(nrow(difference_test(bs, type = "edge")), 45)
})

test_that("difference_test edge difference brackets the observed difference", {
  skip_slow()
  bs <- net_boot(mk(7), n_boot = 80, cores = 1)
  dt <- difference_test(bs, type = "strength")
  # the observed difference lies inside its own bootstrap interval, allowing
  # for a small percentile margin
  inside <- dt$obs_diff >= dt$lower - 1e-8 & dt$obs_diff <= dt$upper + 1e-8
  expect_gt(mean(inside), 0.8)
})

test_that("net_stability returns CS-coefficients in [0,1]", {
  skip_slow()
  cs <- net_stability(mk(4), drop_prop = c(0.3, 0.5, 0.7), iter = 15)
  expect_s3_class(cs, "psychnet_stability")
  expect_true(all(cs$cs >= 0 & cs$cs <= 1))
  expect_named(cs$table,
               c("measure", "drop_prop", "mean_cor", "sd_cor", "prop_above"))
  # larger drop proportions never increase stability
  str_tab <- cs$table[cs$table$measure == "strength", ]
  expect_true(str_tab$mean_cor[1] + 1e-9 >= str_tab$mean_cor[nrow(str_tab)])
})

test_that("net_boot enrichments: extra centralities, predictability, threshold, diff_test", {
  skip_slow()
  set.seed(3)
  x <- matrix(stats::rnorm(200 * 5), 200, 5) %*%
    chol(0.4^abs(outer(1:5, 1:5, "-")))
  colnames(x) <- paste0("V", 1:5)
  bs <- net_boot(x, n_boot = 60, cores = 1,
                 measures = c("strength", "expected_influence",
                              "betweenness", "closeness"),
                 predictability = TRUE, threshold = TRUE,
                 diff_test = TRUE, p_adjust = "BH")
  expect_true(all(c("betweenness", "betweenness_lower", "closeness_upper") %in%
                    names(bs$centrality)))
  expect_named(bs$centrality_boot,
               c("strength", "expected_influence", "betweenness", "closeness"))
  expect_s3_class(bs$predictability, "data.frame")
  expect_true(all(bs$predictability$value >= 0 & bs$predictability$value <= 1))
  expect_equal(dim(bs$thresholded), c(5L, 5L))
  expect_equal(dim(bs$edge_diff_p), c(10L, 10L))
  expect_named(bs$centrality_diff_p,
               c("strength", "expected_influence", "betweenness", "closeness"))
  # difference_test exposes a two-sided p-value and works on any measure
  dt <- difference_test(bs, type = "betweenness", p_adjust = "holm")
  expect_true("p_value" %in% names(dt))
  expect_true(all(dt$p_value >= 0 & dt$p_value <= 1))
})

# --- review fixes 2026-08-17: resampling verbs hand each draw to the estimator
# exactly as the full data would be (factor columns and incomplete rows kept)

.count_warnings <- function(expr, class = "warning") {
  n <- 0L
  val <- withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, class)) n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  list(value = val, n = n)
}

.gen_frame <- function(seed, n = 150) {
  set.seed(seed)
  d <- data.frame(g1 = stats::rnorm(n), g2 = stats::rnorm(n),
                  b = stats::rbinom(n, 1, 0.5),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$g2 <- d$g2 + 0.5 * d$g1
  d
}

test_that("a factor node survives every resampling verb for mgm, and is dropped for glasso, like a direct call", {
  skip_if_not_installed("glmnet")
  d <- .gen_frame(5)
  # mgm(native = FALSE) keeps `f` as one categorical node, in the observed fit
  # and in every draw ...
  bs <- net_boot(d, method = "mgm", n_boot = 3, cores = 1, native = FALSE)
  expect_equal(bs$node_labels, c("g1", "g2", "b", "f"))
  expect_true("f" %in% c(bs$edges$from, bs$edges$to))
  st <- net_stability(d, method = "mgm", iter = 2, drop_prop = c(0.2, 0.4),
                      native = FALSE)
  expect_s3_class(st, "psychnet_stability")
  cd <- net_casedrop_reliability(d, method = "mgm", iter = 2, drop_prop = 0.3,
                                 native = FALSE)
  expect_s3_class(cd, "psychnet_casedrop")
  sp <- net_split_reliability(d, method = "mgm", iter = 2, native = FALSE)
  expect_s3_class(sp, "psychnet_reliability")
  # ... while glasso drops it, silently, exactly as ebic_glasso(d) does.
  expect_silent(bs_g <- net_boot(d, method = "glasso", n_boot = 2, cores = 1))
  expect_equal(bs_g$node_labels, c("g1", "g2", "b"))
  expect_equal(ebic_glasso(d)$nodes$label, c("g1", "g2", "b"))
})

test_that("resampling draws from the rows the estimator would use (na_method = 'pairwise')", {
  d <- .gen_frame(7)[, c("g1", "g2", "b")]
  d$g1[1:15] <- NA
  # pairwise (every estimator's default): the observed fit and the draws use all
  # 150 rows, as ebic_glasso(d) itself does; listwise: the estimator drops the
  # incomplete rows inside each fit.
  bs <- net_boot(d, method = "glasso", n_boot = 2, cores = 1)
  expect_identical(bs$observed$weights, ebic_glasso(d)$weights)
  expect_identical(bs$observed$na_method, "pairwise")
  bs_lw <- net_boot(d, method = "glasso", n_boot = 2, cores = 1,
                    na_method = "listwise")
  expect_identical(bs_lw$observed$weights,
                   ebic_glasso(d, na_method = "listwise")$weights)
  expect_false(identical(bs$observed$weights, bs_lw$observed$weights))
  # complete numeric data: identical draws to a numeric-matrix input
  dc <- .gen_frame(7)[, c("g1", "g2", "b")]
  set.seed(1); a <- net_boot(dc, method = "glasso", n_boot = 3, cores = 1)
  set.seed(1); b <- net_boot(as.matrix(dc), method = "glasso", n_boot = 3, cores = 1)
  expect_identical(a$edges, b$edges)
})

test_that("mgm node types are decided once on the full data and pinned across draws", {
  skip_if_not_installed("glmnet")
  # `k` has 11 distinct integer values in full (gaussian by mgm's rule) but a
  # bootstrap draw can hold <= 10 of them; without pinning, that draw would
  # silently re-classify it as categorical.
  set.seed(8); n <- 60
  d <- data.frame(g1 = stats::rnorm(n), g2 = stats::rnorm(n),
                  k = c(0:10, sample(0:10, n - 11, TRUE)))
  d$g2 <- d$g2 + 0.5 * d$g1
  r <- .count_warnings(net_boot(d, method = "mgm", n_boot = 12, cores = 1),
                       "psychnet_type_promotion")
  expect_equal(r$n, 0L)                                # never promoted
  expect_equal(unname(r$value$observed$types), c("g", "g", "g"))
  expect_equal(nrow(r$value$edges), 3L)               # every draw fit p = 3
  # and a genuinely promoted column warns exactly once per verb call
  d$k <- sample(1:3, n, TRUE)
  for (verb in list(
    function() net_boot(d, method = "mgm", n_boot = 4, cores = 1, native = FALSE),
    function() net_stability(d, method = "mgm", iter = 2, drop_prop = c(0.2, 0.4),
                             native = FALSE),
    function() net_casedrop_reliability(d, method = "mgm", iter = 2,
                                        drop_prop = 0.3, native = FALSE),
    function() net_split_reliability(d, method = "mgm", iter = 3, native = FALSE)
  )) expect_equal(.count_warnings(verb(), "psychnet_type_promotion")$n, 1L)
})
