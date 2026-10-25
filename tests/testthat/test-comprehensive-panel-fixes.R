# Regressions from the complete panel-function review. These cases exercise
# public output and plotting, including panels with no selected edges.

review_panel_empty <- function() {
  k <- 8L
  d <- as.data.frame(rbind(diag(k), -diag(k)) * sqrt(k))
  panel_gvar(d, vars = matrix(names(d), 2L), select = "ebic", nlambda = 4L)
}

test_that("missing person IDs cannot combine unrelated people across waves", {
  x <- panel_sim(n = 100)
  x$id[c(1, 102, 203, 304)] <- NA
  expect_error(panel_gvar(x, c("V1", "V2", "V3"), "id", "wave"),
               "identifier.*missing", class = "psychnet_bad_input")
})

test_that("nlambda must be a finite whole integer before any fitting", {
  bad <- list(2.5, Inf, -Inf, NA_real_, NaN, numeric(), c(2, 3), 1,
              .Machine$integer.max + 1, 2 + 1i, "5")
  for (n in bad) {
    # No usable data: this must fail at nlambda validation, not model fitting.
    expect_error(panel_gvar(NULL, vars = NULL, nlambda = n),
                 "nlambda", class = "psychnet_bad_input")
  }
})

test_that("an empty EBIC panel has readable zero-row edge output", {
  f <- review_panel_empty()
  edges <- as.data.frame(f)
  expect_identical(nrow(edges), 0L)
  expect_type(edges$significant, "logical")
  expect_named(edges, c("network", "from", "to", "weight", "se", "p",
                        "ci_lower", "ci_upper", "significant"))
  sm <- summary(f)
  expect_identical(sm$edges, rep(0L, 3L))
  expect_true(all(is.na(sm$significant)))
  expect_output(print(f), "EBIC-selected")
  all_edges <- as.data.frame(f, include_zero = TRUE)
  expect_identical(nrow(all_edges), 6L)
  expect_true(all(is.na(all_edges$significant)))
})

test_that("perfectly fitting null panels have a finite CFI of one", {
  ft <- as.data.frame(review_panel_empty(), what = "fit")
  expect_lt(abs(ft$chisq), 1e-10)
  expect_lt(abs(ft$baseline_chisq), 1e-10)
  expect_identical(ft$cfi, 1)
  expect_identical(ft$rmsea, 0)
})

test_that("empty EBIC panels plot all three networks successfully", {
  skip_if_not_installed("cograph")
  f <- review_panel_empty()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(f), f)
})

test_that("labels with spaces keep each temporal edge's significance", {
  labels <- c("a b", "c", "a", "b c")
  f <- suppressWarnings(panel_gvar(panel_sim(p = 4),
    c("V1", "V2", "V3", "V4"), "id", "wave", labels = labels))
  edges <- as.data.frame(f)
  drawn <- psychnets:::.pgv_plot_temporal(f$temporal)
  actual <- psychnets:::.pgv_edge_significance(drawn, "temporal", edges)
  expected <- vapply(seq_len(nrow(drawn$edges)), function(i) {
    from <- labels[drawn$edges$from[i]]
    to <- labels[drawn$edges$to[i]]
    edges$significant[edges$network == "temporal" &
                        edges$from == from & edges$to == to]
  }, logical(1))
  expect_identical(actual, expected)
  expect_true(actual[drawn$edges$from == 1L & drawn$edges$to == 2L])
  expect_false(actual[drawn$edges$from == 3L & drawn$edges$to == 4L])
})

test_that("plot abbreviations override defaults and custom styles have an honest key", {
  skip_if_not_installed("cograph")
  f <- panel_gvar(panel_sim(), c("V1", "V2", "V3"), "id", "wave")
  calls <- list()
  keys <- character()
  # Keep splot's true formals so the wrapper must resolve abbreviations as it
  # does for a real plot; inspect the options that actually reach the renderer.
  spy <- function() {
    calls[[length(calls) + 1L]] <<- list(fill = node_fill,
      positive = edge_positive_color, style = edge_style)
    invisible(x)
  }
  formals(spy) <- formals(cograph::splot)
  testthat::local_mocked_bindings(splot = spy, .package = "cograph")
  testthat::local_mocked_bindings(mtext = function(text, ...) {
    keys <<- c(keys, as.character(text))
  }, .package = "graphics")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(f, node_f = "red", edge_pos = "green",
                        edge_sty = 1L), f)
  expect_length(calls, 3L)
  for (call in calls) {
    expect_identical(call$fill, "red")
    expect_identical(call$positive, "green")
    expect_identical(call$style, 1L)
  }
  expect_match(keys, "custom line styles")
  expect_match(keys, "green: positive")
  expect_false(any(grepl("solid: p", keys, fixed = TRUE)))
  expect_error(plot(f, tit = "Custom"), "ambiguous",
               class = "psychnet_bad_input")
})
