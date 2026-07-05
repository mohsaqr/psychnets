# Group mode: one network per level of a grouping variable, and the per-group
# dispatch of the framework verbs (Nestimate netobject_group parity).

# A two-group data set: group A has a strong common factor, B a weak one, so the
# two networks differ in edge strength.
.mk_groups <- function(seed = 1, n = 200) {
  set.seed(seed)
  mk <- function(r) {
    f <- stats::rnorm(n)
    X <- vapply(1:5, function(j) r * f + stats::rnorm(n, 0, .5), numeric(n))
    colnames(X) <- paste0("V", 1:5); X
  }
  data.frame(g = rep(c("A", "B"), each = n), rbind(mk(0.9), mk(0.3)))
}

test_that("psychnet(group=) returns a keyed psychnet_group", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  expect_s3_class(g, "psychnet_group")
  expect_s3_class(g, "netobject_group")               # so cograph::splot grids it
  expect_identical(names(g), c("A", "B"))
  expect_s3_class(g$A, "psychnet")
  expect_identical(attr(g, "group_col"), "g")
  # the grouping column is never itself a node
  expect_false("g" %in% g$A$nodes$label)
})

test_that("psychnet_group prints, stacks edges, and summarises by group", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  ed <- as.data.frame(g)
  expect_true(all(c("group", "from", "to", "weight") %in% names(ed)))
  expect_identical(sort(unique(ed$group)), c("A", "B"))
  sm <- summary(g)
  expect_identical(nrow(sm), 2L)
  # A's common factor is stronger, so its mean absolute edge weight is larger
  expect_gt(sm$mean_abs_weight[sm$group == "A"],
            sm$mean_abs_weight[sm$group == "B"])
})

test_that("net_centralities dispatches per group", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  cc <- net_centralities(g)
  expect_s3_class(cc, "psychnet_result_group")
  df <- as.data.frame(cc)
  expect_true(all(c("group", "node", "strength") %in% names(df)))
  expect_identical(nrow(df), 10L)                     # 2 groups x 5 nodes
  # stronger network -> higher mean strength
  expect_gt(mean(df$strength[df$group == "A"]),
            mean(df$strength[df$group == "B"]))
})

test_that("net_predict dispatches per group using each level's data", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  pr <- as.data.frame(net_predict(g))
  expect_true(all(c("group", "node", "predictability") %in% names(pr)))
  expect_identical(nrow(pr), 10L)
})

test_that("net_boot bootstraps each group and stacks CIs", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  b <- net_boot(g, n_boot = 40)
  expect_s3_class(b, "psychnet_bootstrap_group")
  expect_identical(names(b), c("A", "B"))
  expect_s3_class(b$A, "psychnet_bootstrap")
  bd <- as.data.frame(b)
  expect_true(all(c("group", "from", "to", "lower", "upper") %in% names(bd)))
})

test_that("net_stability dispatches per group", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  s <- net_stability(g, iter = 15)
  expect_s3_class(s, "psychnet_stability_group")
  expect_identical(names(s), c("A", "B"))
})

test_that("net_compare compares two group levels", {
  g <- psychnet(.mk_groups(), method = "glasso", group = "g")
  cmp <- net_compare(g, c("A", "B"), iter = 60)
  expect_s3_class(cmp, "psychnet_nct")
  # defaulting to the two levels when there are exactly two
  cmp2 <- net_compare(g, iter = 30)
  expect_s3_class(cmp2, "psychnet_nct")
})

test_that("group mode validates levels and missing columns", {
  expect_error(psychnet(.mk_groups(), group = "nope"),
               "Group column")
  d1 <- data.frame(g = rep("only", 50), V1 = stats::rnorm(50),
                   V2 = stats::rnorm(50))
  expect_error(psychnet(d1, group = "g"), "at least 2 levels")
  expect_error(net_compare(data.frame(stats::rnorm(10))), "`data2` is required")
})

test_that("group mode honors vars for an id feature table", {
  set.seed(8)
  ft <- data.frame(id = rep(1:60, each = 4), g = rep(c("A", "B"), each = 120),
                   V1 = stats::rnorm(240), V2 = stats::rnorm(240),
                   V3 = stats::rnorm(240), V4 = stats::rnorm(240))
  gf <- psychnet(ft, id = "id", group = "g", vars = V1:V3)
  expect_s3_class(gf, "psychnet_group")
  expect_identical(gf$A$nodes$label, c("V1", "V2", "V3"))  # V4 excluded
  expect_false("id" %in% gf$A$nodes$label)                 # id never a node
})

test_that("a raw event log + group is accepted as-is and is bootstrappable", {
  set.seed(9)
  ev <- data.frame(
    Actor  = rep(paste0("s", 1:120), each = 12),
    Cohort = rep(c("X", "Y"), each = 720),
    Action = sample(c("read", "quiz", "note", "watch"), 1440, replace = TRUE))
  g <- psychnet(ev, actor = "Actor", action = "Action", group = "Cohort")
  expect_s3_class(g, "psychnet_group")
  expect_identical(names(g), c("X", "Y"))
  expect_false("Cohort" %in% g$X$nodes$label)        # group col never a node
  # event-data group reduces to a design matrix -> bootstrap / compare just work
  b <- net_boot(g, n_boot = 30)
  expect_s3_class(b, "psychnet_bootstrap_group")
  cmp <- net_compare(g, c("X", "Y"), iter = 40)
  expect_s3_class(cmp, "psychnet_nct")
})

test_that("group mode works through the front door for any method", {
  g_pcor <- psychnet(.mk_groups(), method = "pcor", group = "g")
  expect_s3_class(g_pcor$A, "psychnet")
  expect_identical(attr(g_pcor, "call")$method, "pcor")
})
