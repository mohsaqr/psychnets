# Regression tests for the 2026-11-07 review of net_compare(), group mode,
# event-data multilevel networks and panel_gvar(). Each block names the
# finding it pins (N* = NCT / groups / multilevel, P* = panel GVAR).

.rv_groups <- function(n = 80, p = 4, seed = 11) {
  set.seed(seed)
  a <- matrix(stats::rnorm(n * p), n, p)
  a <- a + 0.4 * matrix(rowSums(a), n, p)
  b <- matrix(stats::rnorm(n * p), n, p)
  b[, 2] <- b[, 2] + 0.7 * b[, 1]
  b[, 4] <- b[, 4] + 0.5 * b[, 3]
  colnames(a) <- colnames(b) <- paste0("V", seq_len(p))
  data.frame(grp = rep(c("A", "B"), each = n), rbind(a, b))
}

.rv_events <- function(seed = 1, n = 600) {
  set.seed(seed)
  data.frame(
    Actor   = sample(paste0("s", 1:30), n, replace = TRUE),
    Session = sample(1:3, n, replace = TRUE),
    Action  = sample(c("read", "quiz", "note", "watch"), n, replace = TRUE),
    stringsAsFactors = FALSE)
}

# ---- N1: stacked group results surface errors and forward `...` ------------

test_that("N1 as.data.frame() of a result group raises a member's error", {
  registerS3method("as.data.frame", "psn_rv_broken",
                   function(x, ...) stop("cannot convert"))
  bad <- structure(list(A = structure(list(), class = "psn_rv_broken"),
                        B = structure(list(), class = "psn_rv_broken")),
                   class = c("psychnet_rv", "psychnet_result_group"))
  expect_error(as.data.frame(bad), "cannot convert")
})

test_that("N1 as.data.frame() of a result group forwards `...` (e.g. what =)", {
  registerS3method("as.data.frame", "psn_rv_what",
                   function(x, ..., what = "a")
                     data.frame(what = what, stringsAsFactors = FALSE))
  el <- structure(list(), class = "psn_rv_what")
  grp <- structure(list(A = el, B = el),
                   class = c("psychnet_rv", "psychnet_result_group"))
  out <- as.data.frame(grp, what = "b")
  expect_identical(out$group, c("A", "B"))
  expect_identical(out$what, c("b", "b"))
})

test_that("N1 stacked group stability tables are not silently empty", {
  skip_slow()
  # Needs as.data.frame.psychnet_stability (R/stability.R).
  nets <- psychnet(.rv_groups(), group = "grp")
  set.seed(3)
  st <- net_stability(nets, iter = 5, drop_prop = c(0.2, 0.5))
  tab <- as.data.frame(st)
  expect_gt(nrow(tab), 0L)
  expect_setequal(unique(tab$group), c("A", "B"))
  cs <- as.data.frame(st, what = "cs")
  expect_named(cs, c("group", "measure", "cs"))
})

# ---- N2: infinite values are refused, not silently isolated -----------------

test_that("N2 net_compare() refuses infinite values", {
  d <- .rv_groups()
  a <- as.matrix(d[d$grp == "A", -1L])
  b <- as.matrix(d[d$grp == "B", -1L])
  a[5, 2] <- Inf
  expect_error(net_compare(a, b, iter = 3), class = "psychnet_bad_input")
  b[1, 1] <- -Inf
  expect_error(net_compare(b, a, iter = 3), class = "psychnet_bad_input")
})

test_that("N2 a constant variable is still an isolated node", {
  d <- .rv_groups()
  a <- as.matrix(d[d$grp == "A", -1L])
  b <- as.matrix(d[d$grp == "B", -1L])
  a[, 3] <- 1
  set.seed(1)
  expect_warning(fit <- net_compare(a, b, iter = 3),
                 class = "psychnet_constant_variable")
  expect_true(all(fit$nw1[3, ] == 0))
})

test_that("N2 .nct_cor() refuses a non-finite correlation of varying variables", {
  # finite values whose cross-products overflow: cor() returns NaN
  x <- cbind(V1 = c(1, -1, 3, 2) * 1.7e308, V2 = c(1, 2, 3, 5),
             V3 = c(2, 1, 4, 3))
  expect_error(psychnets:::.nct_cor(x), class = "psychnet_bad_input")
})

# ---- N3: explicit gamma / na_method win over the group's settings -----------

test_that("N3 net_compare() on a group honours an explicit gamma and warns", {
  skip_slow()
  nets <- psychnet(.rv_groups(), group = "grp", gamma = 0.25)
  set.seed(2)
  expect_no_warning(own <- net_compare(nets, iter = 3))
  expect_identical(own$params$gamma, 0.25)
  set.seed(2)
  expect_warning(over <- net_compare(nets, iter = 3, gamma = 0),
                 class = "psychnet_estimator_mismatch")
  expect_identical(over$params$gamma, 0)
  # supplying the group's own value is not a mismatch
  set.seed(2)
  expect_no_warning(same <- net_compare(nets, iter = 3, gamma = 0.25))
  expect_identical(same$params$gamma, 0.25)
})

test_that("N3 net_compare() on a group honours an explicit na_method and warns", {
  skip_slow()
  d <- .rv_groups()
  d$V2[c(3, 90)] <- NA
  nets <- psychnet(d, group = "grp")
  set.seed(2)
  expect_no_warning(own <- net_compare(nets, iter = 3))
  expect_identical(own$params$na_method, "pairwise")
  set.seed(2)
  expect_warning(lw <- net_compare(nets, iter = 3, na_method = "listwise"),
                 class = "psychnet_estimator_mismatch")
  expect_identical(lw$params$na_method, "listwise")
  expect_identical(lw$n, c(79L, 79L))
})

# ---- N4: missing actor / session identifiers are refused ---------------------

test_that("N4 events with a missing actor or session are refused", {
  ev <- .rv_events()
  ev$Actor[c(2, 7)] <- NA
  expect_error(event_frequencies(ev, session = "Session"),
               class = "psychnet_bad_input")
  expect_error(psychnet(ev, actor = "Actor", session = "Session"),
               class = "psychnet_bad_input")
  ev <- .rv_events()
  ev$Session[4] <- NA
  expect_error(event_frequencies(ev, session = "Session"),
               class = "psychnet_bad_input")
  # complete identifiers: every event is counted
  ev <- .rv_events()
  fr <- event_frequencies(ev, session = "Session")
  expect_identical(sum(as.matrix(fr[c("note", "quiz", "read", "watch")])),
                   as.numeric(nrow(ev)))
})

# ---- N5: per-level sample-size guards and dropped NA groups ------------------

test_that("N5 a group level with fewer than 3 rows is refused by name", {
  d <- .rv_groups()
  d$grp[1:2] <- "tiny"
  err <- expect_error(psychnet(d, group = "grp"), class = "psychnet_bad_input")
  expect_match(conditionMessage(err), "tiny", fixed = TRUE)
  d$grp[1:2] <- c("one", "B")
  expect_error(psychnet(d, group = "grp"), class = "psychnet_bad_input")
})

test_that("N5 a group level with no more rows than variables warns", {
  d <- .rv_groups()
  d$grp[1:4] <- "small"                  # 4 rows, 4 variables
  expect_warning(psychnet(d, group = "grp"), class = "psychnet_small_sample")
})

test_that("N5 rows with a missing group value are dropped with a warning", {
  d <- .rv_groups()
  d$grp[c(1, 100)] <- NA
  expect_warning(nets <- psychnet(d, group = "grp"),
                 class = "psychnet_rows_dropped")
  expect_identical(vapply(attr(nets, "subsets"), nrow, integer(1)),
                   c(A = 79L, B = 79L))
})

# ---- N6 / N7: group labels reach net_compare; group errors are classed -------

test_that("N6 net_compare() on a group uses the group's labels", {
  skip_slow()
  nets <- psychnet(.rv_groups(), group = "grp", labels = c("a", "b", "c", "d"))
  set.seed(1)
  ed <- as.data.frame(net_compare(nets, iter = 3))
  expect_setequal(unique(c(ed$from, ed$to)), c("a", "b", "c", "d"))
})

test_that("N7 group-mode net_compare() errors are classed", {
  nets <- psychnet(.rv_groups(), group = "grp")
  expect_error(net_compare(nets, c("A", "Z")), class = "psychnet_bad_input")
  expect_error(net_compare(nets, "A"), class = "psychnet_bad_input")
  expect_error(net_compare(nets, c("A", "A")), class = "psychnet_bad_input")
  ev <- .rv_events()
  ev$cohort <- rep(c("x", "y"), length.out = nrow(ev))
  ml <- psychnet(ev, actor = "Actor", session = "Session", group = "cohort",
                 standardize = FALSE)
  expect_error(net_compare(ml), class = "psychnet_bad_input")
})

# ---- N8: summary() and plot() for psychnet_multilevel -----------------------

test_that("N8 summary() of a within/between pair is one row per level", {
  ml <- psychnet(.rv_events(), actor = "Actor", session = "Session",
                 standardize = FALSE)
  sm <- summary(ml)
  expect_s3_class(sm, "data.frame")
  expect_named(sm, c("level", "nodes", "edges", "mean_abs_weight", "n"))
  expect_identical(sm$level, c("within", "between"))
  expect_identical(sm$edges, c(nrow(as.data.frame(ml$within)),
                               nrow(as.data.frame(ml$between))))
  expect_identical(sm$n[2L], attr(ml, "n_actors"))
})

test_that("N8 plot() of a within/between pair draws through cograph", {
  skip_if_not_installed("cograph")
  ml <- psychnet(.rv_events(), actor = "Actor", session = "Session",
                 standardize = FALSE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(ml), ml)
})

# ---- N9: the multilevel-group refusal names a verb, not `$` -----------------

test_that("N9 a multilevel group refusal points to an argument, not `$`", {
  ev <- .rv_events()
  ev$cohort <- rep(c("x", "y"), length.out = nrow(ev))
  ml <- psychnet(ev, actor = "Actor", session = "Session", group = "cohort",
                 standardize = FALSE)
  err <- expect_error(net_centralities(ml), class = "psychnet_bad_input")
  expect_false(grepl("$", conditionMessage(err), fixed = TRUE))
  expect_match(conditionMessage(err), "standardize = TRUE", fixed = TRUE)
})

# ---- N10: event_frequencies() column types -----------------------------------

test_that("N10 event_frequencies() returns character labels and double counts", {
  fr <- event_frequencies(.rv_events(), session = "Session")
  expect_type(fr$actor, "character")
  expect_type(fr$session, "character")
  expect_true(all(vapply(fr[c("note", "quiz", "read", "watch")], is.double,
                         logical(1))))
})

# ---- N11: unparseable `time` raises the package's classed error -------------

test_that("N11 an unparseable time column raises psychnet_bad_input", {
  ev <- .rv_events(n = 40)
  ev$Time <- "not a time"
  expect_error(event_frequencies(ev, time = "Time"), class = "psychnet_bad_input")
  ev$Time <- c(rep("2024-01-01 10:00", 39), "garbage")
  expect_error(event_frequencies(ev, time = "Time"), class = "psychnet_bad_input")
  ev$Time <- format(as.POSIXct("2024-01-01", tz = "UTC") + seq_len(40) * 60)
  expect_s3_class(event_frequencies(ev, time = "Time"), "data.frame")
})

# ---- N12: factor level order and empty-string levels ------------------------

test_that("N12 a factor grouping keeps its level order", {
  d <- .rv_groups()
  d$grp <- factor(d$grp, levels = c("B", "A"))
  nets <- psychnet(d, group = "grp")
  expect_identical(names(nets), c("B", "A"))
  # character grouping keeps order of appearance
  d$grp <- as.character(d$grp)
  expect_identical(names(psychnet(d, group = "grp")), c("A", "B"))
})

test_that("N12 an empty-string group level is refused", {
  d <- .rv_groups()
  d$grp[d$grp == "B"] <- ""
  expect_error(psychnet(d, group = "grp"), class = "psychnet_bad_input")
})

test_that("N12 print of a result group lists every group", {
  nets <- psychnet(.rv_groups(), group = "grp")
  out <- capture.output(print(net_centralities(nets)))
  expect_true(any(grepl("--- A ---", out, fixed = TRUE)))
  expect_true(any(grepl("--- B ---", out, fixed = TRUE)))
})

# ---- panel GVAR -------------------------------------------------------------

.rv_pv <- function(p = 3) paste0("V", seq_len(p))

# Three-wave synthetic panels often end on a between-person boundary; those
# diagnostics are tested in test-panelgvar.R and muffled here, by class only.
suppress_panel_fit_warnings <- function(expr) {
  withCallingHandlers(expr,
    psychnet_boundary = function(w) invokeRestart("muffleWarning"),
    psychnet_singular_information = function(w) invokeRestart("muffleWarning"),
    psychnet_approximate_se = function(w) invokeRestart("muffleWarning"))
}

test_that("P1 numeric character waves are fitted in numeric time order", {
  d <- panel_sim(n = 150, p = 3, waves = 11, seed = 4)
  f_num <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave")
  d$wave <- as.character(d$wave)
  expect_no_warning(f_chr <- panel_gvar(d, vars = .rv_pv(), id = "id",
                                        wave = "wave"),
                    class = "psychnet_wave_order")
  expect_identical(attr(f_chr, "waves"), as.character(1:11))
  expect_equal(as.data.frame(f_chr, what = "fit"),
               as.data.frame(f_num, what = "fit"), tolerance = 1e-10)
})

test_that("P1 non-numeric character waves warn with the order used", {
  d <- panel_sim(n = 150, p = 3, waves = 4, seed = 4)
  d$wave <- c("t1", "t2", "t3", "t4")[d$wave]
  expect_warning(panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave"),
                 class = "psychnet_wave_order")
  d$wave <- factor(d$wave, levels = c("t1", "t2", "t3", "t4"))
  expect_no_warning(suppress_panel_fit_warnings(
    panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave")),
    class = "psychnet_wave_order")
  expect_identical(psychnets:::.pgv_wave_levels(factor(c("b", "a"),
                                                       levels = c("b", "a"))),
                   c("b", "a"))
  expect_identical(psychnets:::.pgv_wave_levels(c("10", "2", "1")),
                   c("1", "2", "10"))
})

test_that("P1 a missing wave value is refused", {
  d <- panel_sim(n = 60, p = 3, waves = 3, seed = 4)
  d$wave[5] <- NA
  expect_error(panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave"),
               class = "psychnet_bad_input")
})

test_that("P2 duplicate or missing labels are refused before fitting", {
  d <- panel_sim(n = 60, p = 3, waves = 3, seed = 4)
  expect_error(panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
                          labels = c("a", "a", "b")),
               class = "psychnet_bad_input")
  expect_error(panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
                          labels = c("a", NA, "b")),
               class = "psychnet_bad_input")
})

test_that("P3 Fisher polish halves the step instead of stopping early", {
  d <- panel_sim(n = 25, p = 4, waves = 3, seed = 2)
  fit <- withCallingHandlers(
    panel_gvar(d, vars = .rv_pv(4), id = "id", wave = "wave",
               select = "prune"),
    psychnet_approximate_se = function(w) invokeRestart("muffleWarning"))
  dg <- as.data.frame(fit, what = "diagnostics")
  expect_false(dg$boundary)
  expect_lt(dg$score, 1e-6)
})

test_that("P4 starts agreement is judged on the unit-free objective", {
  st <- data.frame(start = c("a", "b"), objective = c(1, 1 + 1e-5),
                   stringsAsFactors = FALSE)
  # A large log-Jacobian (big units) must not make distinct optima agree.
  tab <- psychnets:::.pgv_starts_table(st, log_jacobian = 1e4)
  expect_identical(tab$best, c(TRUE, FALSE))
  expect_equal(tab$objective, c(1, 1 + 1e-5) + 1e4)
  d <- panel_sim(n = 150, p = 3, waves = 3, seed = 5)
  big <- d
  big[.rv_pv()] <- big[.rv_pv()] * 1e6
  a <- as.data.frame(panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave"),
                     what = "diagnostics")
  b <- as.data.frame(panel_gvar(big, vars = .rv_pv(), id = "id", wave = "wave"),
                     what = "diagnostics")
  expect_identical(a$starts_agree, b$starts_agree)
})

test_that("P5 adjust without pruning warns that it is ignored", {
  d <- panel_sim(n = 150, p = 3, waves = 4, seed = 4)
  expect_warning(f <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
                                 adjust = "BH"),
                 class = "psychnet_ignored_argument")
  # p-values are unadjusted: identical to the default fit
  f0 <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave")
  expect_identical(as.data.frame(f)$p, as.data.frame(f0)$p)
  expect_no_warning(suppress_panel_fit_warnings(
    panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
               select = "prune", adjust = "BH")),
    class = "psychnet_ignored_argument")
})

test_that("P6/P7 an EBIC refit reports its own starts; path and covariance read tidily", {
  skip_slow()
  d <- panel_sim(n = 200, p = 3, waves = 3, seed = 1)
  fit <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
                    select = "ebic", nlambda = 5, refit = TRUE)
  expect_identical(attr(fit, "starts")$start, c("selected", "baseline"))
  expect_true(any(grepl("refitted (unpenalized) optimum",
                        capture.output(print(fit)), fixed = TRUE)))
  path <- as.data.frame(fit, what = "path")
  expect_named(path, c("lambda", "ebic", "edges", "kkt", "final", "selected"))
  expect_identical(sum(path$selected), 1L)
  expect_equal(path$lambda[path$selected], attr(fit, "settings")$lambda)
  pen <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave",
                    select = "ebic", nlambda = 5)
  expect_true(any(grepl("best penalized optimum",
                        capture.output(print(pen)), fixed = TRUE)))
})

test_that("P7 covariance accessor matches the fitted covariances", {
  d <- panel_sim(n = 150, p = 3, waves = 4, seed = 4)
  fit <- panel_gvar(d, vars = .rv_pv(), id = "id", wave = "wave")
  cv <- as.data.frame(fit, what = "covariance")
  expect_named(cv, c("matrix", "row", "col", "value"))
  expect_identical(nrow(cv), 12L)                   # 2 x 6 lower-triangle cells
  inn <- subset(cv, matrix == "innovation")
  sig <- fit$contemporaneous$sigma
  expect_equal(inn$value, sig[lower.tri(sig, diag = TRUE)])
  expect_error(as.data.frame(fit, what = "path"), class = "psychnet_bad_input")
  expect_error(as.data.frame(fit, what = "nope"), class = "psychnet_bad_input")
})

test_that("P8 bad scalar arguments raise psychnet_bad_input", {
  d <- panel_sim(n = 60, p = 3, waves = 3, seed = 4)
  call <- function(...) panel_gvar(d, vars = .rv_pv(), id = "id",
                                   wave = "wave", ...)
  expect_error(call(alpha = 2), class = "psychnet_bad_input")
  expect_error(call(refit = "yes"), class = "psychnet_bad_input")
  expect_error(call(lambda_min_ratio = NA), class = "psychnet_bad_input")
  expect_error(call(select = "lasso"), class = "psychnet_bad_input")
  expect_error(call(standardize = "quantile"), class = "psychnet_bad_input")
  expect_error(call(criterion = "aic"), class = "psychnet_bad_input")
  expect_error(call(adjust = "nope"), class = "psychnet_bad_input")
  expect_error(call(gamma = -1), class = "psychnet_bad_input")
  expect_error(call(stepup = NA), class = "psychnet_bad_input")
})

test_that("P9 duplicated vars or id/wave among vars are refused at input", {
  d <- panel_sim(n = 60, p = 3, waves = 3, seed = 4)
  expect_error(panel_gvar(d, vars = c("V1", "V1", "V2"), id = "id",
                          wave = "wave"), class = "psychnet_bad_input")
  expect_error(panel_gvar(d, vars = c("V1", "V2", "wave"), id = "id",
                          wave = "wave"), class = "psychnet_bad_input")
  expect_error(panel_gvar(d, vars = c("V1", "V2", "id"), id = "id",
                          wave = "wave"), class = "psychnet_bad_input")
  w <- panel_sim_wide(n = 60, p = 3, waves = 3, seed = 4)
  design <- w$design
  design[2, 1] <- design[1, 1]
  expect_error(panel_gvar(w$data, vars = design), class = "psychnet_bad_input")
})
