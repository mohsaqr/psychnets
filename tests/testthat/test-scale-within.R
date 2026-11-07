# scale_within(): within-person centering and scaling (mlVAR scaleWithin).

v <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")

test_that("scale_within centers and scales each person to mean 0, sd 1", {
  sw <- scale_within(srl_panel, id = "student", vars = v)
  expect_identical(dim(sw), dim(srl_panel))
  expect_identical(sw[c("student", "week")], srl_panel[c("student", "week")])
  means <- vapply(v, function(x) max(abs(tapply(sw[[x]], sw$student, mean))), 1)
  sds <- vapply(v, function(x) max(abs(tapply(sw[[x]], sw$student, stats::sd) - 1)), 1)
  expect_lt(max(means), 1e-12)
  expect_lt(max(sds), 1e-12)
})

test_that("center-only and scale-only follow their definitions", {
  ce <- scale_within(srl_panel, "student", v, scale = FALSE)
  expect_equal(ce$effort, srl_panel$effort - stats::ave(srl_panel$effort, srl_panel$student),
               tolerance = 1e-14)
  sc <- scale_within(srl_panel, "student", v, center = FALSE)
  expect_equal(sc$effort,
               srl_panel$effort / stats::ave(srl_panel$effort, srl_panel$student,
                                             FUN = stats::sd),
               tolerance = 1e-14)
})

test_that("scale_within equals mlVAR's within-person scaling", {
  skip_if_not_installed("mlVAR")
  ave_center <- utils::getFromNamespace("aveCenter", "mlVAR")
  sw <- scale_within(srl_panel, "student", v)
  ref <- lapply(srl_panel[v], function(x)
    stats::ave(x, srl_panel$student, FUN = function(z) ave_center(z, scale = TRUE)))
  expect_equal(as.list(sw[v]), ref, tolerance = 1e-14)
})

test_that("the result is invariant to a per-person shift and positive rescaling", {
  set.seed(4)
  d <- srl_panel
  shift <- stats::setNames(stats::rnorm(100, 0, 20), unique(d$student))
  mult <- stats::setNames(stats::runif(100, 0.5, 3), unique(d$student))
  key <- as.character(d$student)
  d2 <- d
  d2[v] <- lapply(d[v], function(x) x * mult[key] + shift[key])
  expect_equal(scale_within(d2, "student", v)[v], scale_within(d, "student", v)[v],
               tolerance = 1e-10)
})

test_that("two occasions per person warn that deviations are mirror images", {
  two <- subset(srl_panel, week <= 2)
  # (a few students have equal values in both weeks: no variation -> NA)
  expect_warning(
    expect_warning(sw <- scale_within(two, "student", v),
                   class = "psychnet_mirror_deviations"),
    class = "psychnet_constant_within")
  # the two weeks' centered deviations are exact negatives
  expect_equal(subset(sw, week == 1, select = v)$effort,
               -subset(sw, week == 2, select = v)$effort, tolerance = 1e-12)
  expect_no_warning(scale_within(srl_panel, "student", v))
})

test_that("persons without within-person variation get NA with a warning", {
  d <- srl_panel
  d$effort[d$student == d$student[1]] <- 50           # one constant person
  expect_warning(sw <- scale_within(d, "student", v),
                 class = "psychnet_constant_within")
  expect_true(all(is.na(sw$effort[d$student == d$student[1]])))
  expect_false(anyNA(sw$effort[d$student != d$student[1]]))
  one <- rbind(srl_panel, transform(srl_panel[1, ], student = 999L))
  expect_warning(scale_within(one, "student", v, scale = FALSE),
                 class = "psychnet_constant_within")
})

test_that("bad input is a classed error", {
  expect_error(scale_within(as.matrix(srl_panel), "student"), class = "psychnet_bad_input")
  expect_error(scale_within(srl_panel, "nope"), class = "psychnet_bad_input")
  expect_error(scale_within(srl_panel, "student", vars = "nope"),
               class = "psychnet_bad_input")
  d <- srl_panel; d$student[3] <- NA
  expect_error(scale_within(d, "student", v), class = "psychnet_bad_input")
  d <- srl_panel; d$tag <- "x"
  expect_error(scale_within(d, "student", vars = c(v, "tag")),
               class = "psychnet_bad_input")
  expect_error(scale_within(srl_panel, "student", v, scale = NA),
               class = "psychnet_bad_input")
})

test_that("within-scaled occasions feed a paired network comparison", {
  sw <- scale_within(srl_panel, "student", v)
  set.seed(1)
  fit <- net_compare(subset(sw, week == 1, select = v),
                     subset(sw, week == 4, select = v), paired = TRUE, iter = 20)
  expect_s3_class(fit, "psychnet_nct")
  expect_identical(nrow(as.data.frame(fit)), 10L)
})

# ---- psychnet(within =) and person-paired net_compare -----------------------

test_that("psychnet(within =) fits each level on within-person scaled data", {
  nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                   within = "student")
  sw <- scale_within(srl_panel, "student", v)
  expect_identical(nets[["2"]]$weights,
                   ebic_glasso(subset(sw, week == 2, select = v))$weights)
  expect_false("student" %in% nets[["1"]]$nodes$label)
  # without vars: every numeric column except the person and group is a node
  all_num <- psychnet(srl_panel, group = "week", within = "student")
  expect_identical(all_num[["1"]]$nodes$label, v)
})

test_that("a paired comparison of a within group pairs by person, not row order", {
  nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                   within = "student")
  set.seed(1)
  a <- net_compare(nets, c(1, 4), iter = 30)
  set.seed(7)
  shuffled <- srl_panel[sample(nrow(srl_panel)), ]
  nets_s <- psychnet(shuffled, vars = efficacy:anxiety, group = "week",
                     within = "student")
  set.seed(1)
  b <- net_compare(nets_s, c(1, 4), iter = 30)
  # same pairs and same permutation draws; the numbers agree to float noise
  # (the per-person means are summed in a different order after shuffling)
  expect_equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-10)
  expect_equal(summary(a), summary(b), tolerance = 1e-10)
  expect_identical(as.data.frame(a)$p_value, as.data.frame(b)$p_value)
})

test_that("persons seen in one level only are left out with a warning", {
  d <- subset(srl_panel, !(student == "S001" & week == 1))
  nets <- psychnet(d, vars = efficacy:anxiety, group = "week", within = "student")
  set.seed(1)
  expect_warning(cmp <- net_compare(nets, c(1, 4), iter = 10),
                 class = "psychnet_rows_dropped")
  expect_identical(cmp$n, c(99L, 99L))
})

test_that("within = refuses event logs and repeated persons in a paired level", {
  ev <- data.frame(Actor = rep(c("a", "b"), each = 10),
                   Action = rep(c("x", "y"), 10))
  expect_error(psychnet(ev, actor = "Actor", within = "Actor"),
               class = "psychnet_bad_input")
  expect_error(psychnet(srl_panel, vars = efficacy:anxiety, within = "nope"),
               class = "psychnet_bad_input")
  d <- srl_panel
  d$block <- ifelse(d$week <= 2, "early", "late")   # two rows per student per level
  nets <- psychnet(d, vars = efficacy:anxiety, group = "block", within = "student")
  expect_error(net_compare(nets, paired = TRUE, iter = 10),
               class = "psychnet_bad_input")
})

test_that("stacked group tables have plain row names", {
  nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                   within = "student")
  expect_identical(rownames(as.data.frame(nets)),
                   as.character(seq_len(nrow(as.data.frame(nets)))))
})

test_that("a within group is compared paired by person by default", {
  nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                   within = "student")
  set.seed(1); auto <- net_compare(nets, c(1, 4), iter = 10)
  set.seed(1); expl <- net_compare(nets, c("1", "4"), paired = TRUE, iter = 10)
  expect_true(auto$paired)
  expect_identical(as.data.frame(auto), as.data.frame(expl))
  set.seed(1); unpaired <- net_compare(nets, c(1, 4), paired = FALSE, iter = 10)
  expect_false(unpaired$paired)
  # an ordinary group keeps the independent-groups default
  d <- rbind(data.frame(SRL_GPT, src = "GPT"), data.frame(SRL_Claude, src = "Claude"))
  set.seed(1)
  expect_false(net_compare(psychnet(d, group = "src"), iter = 10)$paired)
})

test_that("two networks of one group can be compared directly", {
  nets <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week",
                   within = "student")
  set.seed(1); by_level <- net_compare(nets, c(1, 4), iter = 10)
  set.seed(1); by_fit <- net_compare(nets$`1`, nets$`4`, iter = 10)
  expect_identical(as.data.frame(by_fit), as.data.frame(by_level))
  expect_identical(summary(by_fit), summary(by_level))
  expect_true(by_fit$paired)                         # within group: paired
  expect_identical(by_fit$groups, c("1", "4"))
  # an ordinary two-level group
  d <- rbind(data.frame(SRL_GPT, src = "GPT"), data.frame(SRL_Claude, src = "Claude"))
  g <- psychnet(d, group = "src")
  set.seed(2); a <- net_compare(g$GPT, g$Claude, iter = 10)
  set.seed(2); b <- net_compare(g, iter = 10)
  expect_identical(as.data.frame(a), as.data.frame(b))
  # refusals: a standalone fit, the same level twice, two different calls
  expect_error(net_compare(ebic_glasso(SRL_GPT), ebic_glasso(SRL_Claude), iter = 5),
               class = "psychnet_bad_input")
  expect_error(net_compare(nets$`1`, nets$`1`, iter = 5), class = "psychnet_bad_input")
  other <- psychnet(srl_panel, vars = efficacy:anxiety, group = "week")
  expect_error(net_compare(nets$`1`, other$`4`, iter = 5), class = "psychnet_bad_input")
})
