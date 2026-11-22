# psychnet() with `group`, `within`, `id` and event logs; the group and multilevel accessors; and
# event_frequencies(). References: qgraph::EBICglasso() per group; the definitions in base R (partial
# correlations per group, within-person centring and scaling, the within/between decomposition of
# Hamaker and Muthen, counts of events); mlVAR's within-person scaling.

with_gender <- utils::read.csv(reference_data("data/bfi10_gender.csv"))
items <- setdiff(names(with_gender), "gender")
partials <- \(x) -stats::cov2cor(solve(stats::cor(x)))
by_gender <- split(with_gender[items], with_gender$gender)

test_that("one partial-correlation network per group equals the definition on each group's rows", {
  groups <- psychnet(with_gender, method = "pcor", group = "gender")
  expected <- do.call(rbind, lapply(names(by_gender), \(level) data.frame(group = level, ref_edges(partials(by_gender[[level]])))))
  expect_reference(as.data.frame(groups, include_zero = TRUE), expected,
                   reference = "independent: -cov2cor(solve(cor(x))) on each group's rows, base R",
                   covers = c("psychnet()", "as.data.frame.psychnet_group()"), quantity = "partial correlation per gender",
                   tolerance = 1e-12, key = c("group", "from", "to"))
  counts <- data.frame(group = names(by_gender), edges = vapply(by_gender, \(x) sum(abs(partials(x)[upper.tri(diag(10))]) > 0),
                                                                 numeric(1)),
                       mean_abs_weight = vapply(by_gender, \(x) mean(abs(partials(x)[upper.tri(diag(10))])), numeric(1)))
  expect_reference(summary(groups), counts, reference = "independent: -cov2cor(solve(cor(x))) on each group's rows, base R",
                   covers = "summary.psychnet_group()", quantity = "edges and mean absolute weight per gender",
                   tolerance = 1e-12, key = "group")
  strength <- do.call(rbind, lapply(names(by_gender), \(level)
    data.frame(group = level, node = items, strength = rowSums(abs(partials(by_gender[[level]]))) - 1)))
  expect_reference(as.data.frame(net_centralities(groups, measures = "strength")), strength,
                   reference = "independent: -cov2cor(solve(cor(x))) on each group's rows, base R",
                   covers = c("as.data.frame.psychnet_result_group()", "net_centralities()"),
                   quantity = "strength per gender", tolerance = 1e-12, key = c("group", "node"),
                   note = "Strength is the row sum of |partial correlation| (the diagonal, 1, is subtracted).")
})

test_that("one EBIC graphical lasso network per group equals qgraph::EBICglasso() on each group's rows", {
  skip_without_reference("qgraph")
  expected <- do.call(rbind, lapply(names(by_gender), \(level) {
    rows <- by_gender[[level]]
    network <- qgraph::EBICglasso(stats::cor(rows), n = nrow(rows), gamma = 0.5, refit = FALSE, threshold = FALSE,
                                  thr = 1e-10, maxit = 1e5, verbose = FALSE)
    data.frame(group = level, ref_edges(network))
  }))
  # Two independent solvers of one convex problem per group: 1e-6, as for ebic_glasso(). The edges are
  # compared both ways: an edge present on one side only fails.
  expect_reference(as.data.frame(psychnet(with_gender, method = "glasso", group = "gender")),
                   subset(expected, weight != 0), reference = "qgraph", covers = c("psychnet()", "as.data.frame.psychnet_group()"),
                   quantity = "partial correlation per gender, EBIC glasso", tolerance = 1e-6,
                   key = c("group", "from", "to"))
})

test_that("as.data.frame() of a group passes include_zero = TRUE on to each network", {
  groups <- psychnet(with_gender, method = "glasso", group = "gender")
  pairs <- length(by_gender) * choose(length(items), 2L)
  expect_reference(nrow(as.data.frame(groups, include_zero = TRUE)), pairs,
                   reference = "independent: two groups times choose(10, 2) node pairs", covers = "as.data.frame.psychnet_group()",
                   quantity = "rows of as.data.frame(include_zero = TRUE), two groups", tolerance = 0,
                   defect = "as.data.frame.psychnet_group() does not pass `...` to each network: include_zero = TRUE is ignored and zero edges are dropped; FINDINGS.md")
})

test_that("within-person scaling equals mlVAR's and the definition", {
  variables <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
  # Each person's mean subtracted and the result divided by the person's standard deviation.
  scaled <- as.data.frame(lapply(srl_panel[variables], \(x)
    stats::ave(x, srl_panel$student, FUN = \(z) (z - mean(z)) / stats::sd(z))))
  expect_reference(as.data.frame(psychnet(srl_panel, vars = efficacy:anxiety, method = "cor", within = "student"),
                                 include_zero = TRUE),
                   ref_edges(stats::cor(scaled)), reference = "independent: within-person centring and scaling, base R",
                   covers = "psychnet()", quantity = "correlation of within-person scaled scores", tolerance = 1e-12,
                   key = c("from", "to"))
  skip_without_reference("mlVAR")
  # mlVAR::mlVAR(scaleWithin = TRUE) scales each person's series with its internal aveCenter(scale = TRUE).
  ave_center <- utils::getFromNamespace("aveCenter", "mlVAR")
  by_mlvar <- as.data.frame(lapply(srl_panel[variables], \(x)
    stats::ave(x, srl_panel$student, FUN = \(z) ave_center(z, scale = TRUE))))
  expect_reference(as.data.frame(psychnet(srl_panel, vars = efficacy:anxiety, method = "cor", within = "student"),
                                 include_zero = TRUE),
                   ref_edges(stats::cor(by_mlvar)), reference = "mlVAR", covers = "psychnet()",
                   quantity = "correlation of within-person scaled scores", tolerance = 1e-12, key = c("from", "to"))
})

test_that("within and between networks equal partial correlations of deviations and of person means", {
  variables <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
  fit <- psychnet(srl_panel, id = "student", vars = efficacy:anxiety, method = "pcor", standardize = FALSE)
  values <- as.matrix(srl_panel[variables])
  means <- rowsum(values, srl_panel$student) / as.vector(table(srl_panel$student))
  deviations <- values - means[as.character(srl_panel$student), ]
  decomposition <- "independent: deviations from and means of each person's rows, partial correlations, base R"
  expected <- rbind(data.frame(level = "within", ref_edges(partials(deviations))),
                    data.frame(level = "between", ref_edges(partials(means))))
  expect_reference(as.data.frame(fit, include_zero = TRUE), expected, reference = decomposition,
                   covers = c("psychnet()", "as.data.frame.psychnet_multilevel()"),
                   quantity = "within- and between-person partial correlation", tolerance = 1e-12,
                   key = c("level", "from", "to"))
  sizes <- data.frame(level = c("within", "between"), n = c(nrow(values) - nrow(means), nrow(means)),
                      mean_abs_weight = c(mean(abs(partials(deviations)[upper.tri(diag(5))])),
                                          mean(abs(partials(means)[upper.tri(diag(5))]))))
  expect_reference(summary(fit), sizes, reference = decomposition, covers = "summary.psychnet_multilevel()",
                   quantity = "sample size and mean absolute weight per level", tolerance = 1e-12, key = "level")
})

test_that("event_frequencies() counts each actor's actions per session", {
  set.seed(4)
  events <- data.frame(Actor = rep(sprintf("s%02d", seq_len(30L)), each = 20L),
                       Session = rep(rep(seq_len(4L), each = 5L), 30L),
                       Action = sample(c("read", "quiz", "note", "watch"), 600L, replace = TRUE))
  counted <- as.data.frame.matrix(table(paste(events$Actor, events$Session, sep = "."), events$Action))
  expected <- data.frame(session = rownames(counted), counted, check.names = FALSE)
  expect_reference(event_frequencies(events, session = "Session"), expected,
                   reference = "independent: a table of actions by actor and session, base R",
                   covers = "event_frequencies()", quantity = "count of each action per session", tolerance = 0,
                   key = "session")
})

test_that("within and between networks equal partial correlations of psych::statsBy()'s rwg and rbg", {
  skip_without_reference(c("psych", "corpcor"))
  variables <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
  fit <- psychnet(srl_panel, id = "student", vars = efficacy:anxiety, method = "pcor", standardize = FALSE)
  # statsBy(): rwg is the pooled within-person correlation, rbg the correlation of the person means.
  by_person <- psych::statsBy(srl_panel[c("student", variables)], group = "student", cors = TRUE)
  named <- \(m) `dimnames<-`(corpcor::cor2pcor(m), list(variables, variables))
  within <- named(by_person$rwg)
  between <- named(by_person$rbg)
  expected <- rbind(data.frame(level = "within", ref_edges(within)), data.frame(level = "between", ref_edges(between)))
  expect_reference(as.data.frame(fit, include_zero = TRUE), expected, reference = "psych",
                   covers = c("psychnet()", "as.data.frame.psychnet_multilevel()"),
                   quantity = "within- and between-person partial correlation, statsBy() rwg / rbg", tolerance = 1e-10,
                   key = c("level", "from", "to"), note = "corpcor::cor2pcor() of psych::statsBy()'s rwg and rbg")
  persons <- nrow(by_person$n)
  sizes <- data.frame(level = c("within", "between"), n = c(sum(by_person$n[, "efficacy"]) - persons, persons),
                      mean_abs_weight = c(mean(abs(within[upper.tri(within)])), mean(abs(between[upper.tri(between)]))))
  expect_reference(summary(fit), sizes, reference = "psych", covers = "summary.psychnet_multilevel()",
                   quantity = "sample size and mean absolute weight per level, statsBy()", tolerance = 1e-10,
                   key = "level")
})

test_that("each group's strength and summary equal qgraph and corpcor on that group's rows", {
  skip_without_reference(c("qgraph", "corpcor"))
  groups <- psychnet(with_gender, method = "pcor", group = "gender")
  per_group <- lapply(by_gender, \(x) `dimnames<-`(corpcor::cor2pcor(stats::cor(x)), list(items, items)))
  strength <- do.call(rbind, lapply(names(per_group), \(level) {
    network <- per_group[[level]]
    diag(network) <- 0
    data.frame(group = level, node = items, strength = unname(qgraph::centrality(network)$OutDegree))
  }))
  expect_reference(as.data.frame(net_centralities(groups, measures = "strength")), strength, reference = "qgraph",
                   covers = c("as.data.frame.psychnet_result_group()", "net_centralities()"),
                   quantity = "strength per gender, qgraph::centrality() of corpcor's partial correlations",
                   tolerance = 1e-10, key = c("group", "node"))
  glasso <- psychnet(with_gender, method = "glasso", group = "gender")
  expected <- do.call(rbind, lapply(names(by_gender), \(level) {
    rows <- by_gender[[level]]
    network <- qgraph::EBICglasso(stats::cor(rows), n = nrow(rows), gamma = 0.5, refit = FALSE, threshold = FALSE,
                                  thr = 1e-10, maxit = 1e5, verbose = FALSE)
    weights <- network[upper.tri(network)]
    data.frame(group = level, edges = sum(weights != 0), mean_abs_weight = mean(abs(weights[weights != 0])))
  }))
  # Two solvers of one convex problem (1e-6, as for ebic_glasso()); the edge count must match exactly.
  expect_reference(summary(glasso), expected, reference = "qgraph", covers = "summary.psychnet_group()",
                   quantity = "edges and mean absolute weight per gender, qgraph::EBICglasso()", tolerance = 1e-6,
                   key = "group")
})
