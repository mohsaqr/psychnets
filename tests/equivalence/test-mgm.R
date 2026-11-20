# mgm_fit(), condition() and net_predict() on a mixed graphical model: Gaussian, binary and multi-level
# categorical nodes, simulated, against mgm::mgm(), mgm::condition() and mgm's predict(). mgm leaves an edge
# between two categorical nodes unsigned, so edge MAGNITUDES are compared. With native = FALSE psychnets fits
# each nodewise lasso with glmnet, as mgm does; the selection, scaling, AND rule and aggregation over a
# categorical node's parameters are each package's own code.

ar1 <- \(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

#' Two Gaussian and three binary nodes driven by a latent AR(1) normal
mixed_data <- function() {
  set.seed(3)
  latent <- matrix(stats::rnorm(2000L * 5L), 2000L, 5L) %*% chol(ar1(5L, 0.5))
  data.frame(g1 = latent[, 1L], g2 = latent[, 2L], b1 = (latent[, 3L] > 0) * 1, b2 = (latent[, 4L] > 0) * 1,
             b3 = (latent[, 5L] > 0) * 1)
}

#' Three Gaussian, one binary and one k-level categorical node on one factor
categorical_data <- function(levels, seed) {
  set.seed(seed)
  factor_score <- stats::rnorm(400L)
  data.frame(g1 = factor_score + stats::rnorm(400L, sd = 0.6), g2 = factor_score + stats::rnorm(400L, sd = 0.6),
             b1 = ((factor_score + stats::rnorm(400L, sd = 0.6)) > 0) * 1L,
             fc = factor(cut(factor_score + stats::rnorm(400L, sd = 0.8), breaks = levels, labels = FALSE)),
             g3 = stats::rnorm(400L))
}

magnitudes <- \(fit) within(as.data.frame(fit, include_zero = TRUE), weight <- abs(weight))
mgm_edges <- \(reference, data) {
  wadj <- reference$pairwise$wadj
  dimnames(wadj) <- list(names(data), names(data))
  ref_edges(wadj)
}

test_that("native = FALSE: edge magnitudes equal mgm::mgm() on Gaussian and binary nodes", {
  skip_without_reference(c("mgm", "glmnet"))
  data <- mixed_data()
  reference <- mgm::mgm(as.matrix(data), type = c("g", "g", "c", "c", "c"), level = c(1, 1, 2, 2, 2),
                        lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW", scale = TRUE,
                        pbar = FALSE, signInfo = FALSE)
  # glmnet's convergence threshold is shared; the package's own test held this to 1e-6.
  expect_reference(magnitudes(mgm_fit(data, gamma = 0.25, native = FALSE)), mgm_edges(reference, data),
                   reference = "mgm", covers = "mgm_fit()", quantity = "edge magnitude, Gaussian and binary nodes",
                   tolerance = 1e-6, key = c("from", "to"))
})

test_that("native = FALSE: edge magnitudes equal mgm::mgm() with 3- and 4-level categorical nodes", {
  skip_without_reference(c("mgm", "glmnet"))
  cases <- expand.grid(levels = c(3L, 4L), seed = c(1L, 2L))
  invisible(lapply(seq_len(nrow(cases)), \(k) {
    data <- categorical_data(cases$levels[[k]], cases$seed[[k]])
    numeric_data <- transform(data, fc = as.integer(fc) - 1L)
    reference <- mgm::mgm(as.matrix(numeric_data), type = c("g", "g", "c", "c", "g"),
                          level = c(1, 1, 2, cases$levels[[k]], 1), lambdaSel = "EBIC", lambdaGam = 0.25,
                          ruleReg = "AND", threshold = "LW", scale = TRUE, pbar = FALSE, signInfo = FALSE)
    # The package's own test held these to 1e-10: both sides run the same glmnet fits.
    expect_reference(magnitudes(mgm_fit(data, gamma = 0.25, native = FALSE, threshold = "LW", rule = "AND")),
                     mgm_edges(reference, data), reference = "mgm", covers = "mgm_fit()",
                     quantity = sprintf("edge magnitude, %d-level categorical node, seed %d", cases$levels[[k]],
                                        cases$seed[[k]]),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("native = TRUE scans its own penalty grid: a known difference from mgm", {
  skip_without_reference("mgm")
  data <- mixed_data()
  reference <- mgm::mgm(as.matrix(data), type = c("g", "g", "c", "c", "c"), level = c(1, 1, 2, 2, 2),
                        lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW", scale = TRUE,
                        pbar = FALSE, signInfo = FALSE)
  expect_known_difference(magnitudes(mgm_fit(data, gamma = 0.25)), mgm_edges(reference, data), reference = "mgm",
                          covers = "mgm_fit()", quantity = "edge magnitude, native solver", tolerance = 1e-6,
                          key = c("from", "to"),
                          reason = paste("The native solver scans its own penalty grid per node (lambda_min_ratio = 0.01);",
                                         "mgm scans glmnet's default path. The EBIC chooses among different candidate",
                                         "penalties, so the selected lasso fits differ in the third decimal."),
                          agrees_with = "native = FALSE (glmnet's path), compared above")
})

test_that("condition() on a moderated model equals mgm::condition()", {
  skip_without_reference(c("mgm", "glmnet"))
  set.seed(2)
  moderator <- stats::rnorm(400L, 10, 2)
  x1 <- stats::rnorm(400L)
  data <- data.frame(x1 = x1, x2 = stats::rnorm(400L), y = x1 * (moderator - 10) / 2 + stats::rnorm(400L),
                     m = moderator)
  fit <- mgm_fit(data, moderators = 4, native = FALSE)
  reference <- mgm::mgm(as.matrix(data), type = rep("g", 4L), level = rep(1, 4L), moderators = 4,
                        lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW",
                        overparameterize = FALSE, scale = TRUE, pbar = FALSE, signInfo = FALSE, warnings = FALSE)
  # condition() takes the moderator on its raw scale; mgm (scale = TRUE) on the z scale.
  invisible(lapply(c(8, 10, 13), \(value) {
    conditioned <- mgm::condition(reference, values = list("4" = (value - mean(moderator)) / stats::sd(moderator)))
    expect_reference(magnitudes(condition(fit, value)), mgm_edges(conditioned, data), reference = "mgm",
                     covers = c("condition()", "mgm_fit()"), quantity = sprintf("edge magnitude at moderator = %g", value),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("node predictability of a mixed graphical model equals mgm's predict()", {
  skip_without_reference(c("mgm", "glmnet"))
  data <- mixed_data()
  fit <- mgm_fit(data, gamma = 0.25, native = FALSE)
  reference <- mgm::mgm(as.matrix(data), type = c("g", "g", "c", "c", "c"), level = c(1, 1, 2, 2, 2),
                        lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW", scale = TRUE,
                        pbar = FALSE, signInfo = FALSE)
  errors <- stats::predict(reference, data = as.matrix(data), errorCon = "R2", errorCat = c("CC", "nCC"),
                           errordecimals = 15)$errors
  expected <- data.frame(node = errors$Variable, predictability = ifelse(is.na(errors$R2), errors$nCC, errors$R2),
                         accuracy = errors$CC)
  # The same model on both sides (compared above to 1e-6), so the predictions agree to that order.
  expect_reference(net_predict(fit, data = data), expected, reference = "mgm", covers = "net_predict()",
                   quantity = "predictability (R2, nCC) and accuracy (CC)", tolerance = 1e-6, key = "node")
  expect_reference(node_predictability(fit, data = data), stats::setNames(expected$predictability, expected$node),
                   reference = "mgm", covers = "node_predictability()", quantity = "predictability as a vector",
                   tolerance = 1e-6)
})

test_that("as.data.frame() of a moderated model gives mgm's aggregated pairwise and three-way weights", {
  # mgm stores the selected interactions in $interactions: indicator[[1]] the pairs and indicator[[2]] the
  # triples, weightsAgg[[1]] and [[2]] their aggregated weights (the mean absolute nodewise estimate). A pair's
  # main effect is its pairwise weight; its moderation effect is the weight of the triple of the pair and the
  # moderator, and is not defined for a pair that holds the moderator.
  skip_without_reference(c("mgm", "glmnet"))
  set.seed(2)
  moderator <- stats::rnorm(400L, 10, 2)
  x1 <- stats::rnorm(400L)
  x2 <- stats::rnorm(400L)
  # x1 - y moderated by m, x2 - y a main effect
  data <- data.frame(x1 = x1, x2 = x2, y = x1 * (moderator - 10) / 2 + 0.5 * x2 + stats::rnorm(400L), m = moderator)
  fit <- mgm_fit(data, moderators = 4, native = FALSE)
  reference <- mgm::mgm(as.matrix(data), type = rep("g", 4L), level = rep(1, 4L), moderators = 4,
                        lambdaSel = "EBIC", lambdaGam = 0.25, ruleReg = "AND", threshold = "LW",
                        overparameterize = FALSE, scale = TRUE, pbar = FALSE, signInfo = FALSE, warnings = FALSE)
  labels <- names(data)
  key_of <- \(i, j) paste(labels[pmin(i, j)], labels[pmax(i, j)])
  pairs <- matrix(reference$interactions$indicator[[1L]], ncol = 2L)
  triples <- matrix(reference$interactions$indicator[[2L]], ncol = 3L)
  main <- stats::setNames(unlist(reference$interactions$weightsAgg[[1L]]), key_of(pairs[, 1L], pairs[, 2L]))
  others <- t(apply(triples, 1L, \(t) sort(setdiff(t, 4L))))
  moderation <- stats::setNames(unlist(reference$interactions$weightsAgg[[2L]]), key_of(others[, 1L], others[, 2L]))
  all_pairs <- utils::combn(4L, 2L)
  keys <- key_of(all_pairs[1L, ], all_pairs[2L, ])
  expected <- data.frame(from = labels[all_pairs[1L, ]], to = labels[all_pairs[2L, ]],
                         main_effect = ifelse(keys %in% names(main), main[keys], 0),
                         interaction = ifelse(all_pairs[2L, ] == 4L, NA, ifelse(keys %in% names(moderation), moderation[keys], 0)))
  expected <- expected[expected$main_effect > 0 | (!is.na(expected$interaction) & expected$interaction > 0), ]
  expect_reference(as.data.frame(fit), expected, reference = "mgm", covers = "as.data.frame.psychnet_moderated()",
                   quantity = "main and moderation effects", tolerance = 1e-10, key = c("from", "to"))
})
