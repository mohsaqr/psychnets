# ising_fit(), ising_sampler(), glm_lasso_kkt() and dichotomize(): Ising models of six binary bfi items
# (A1-A5 and C1, 1 = agree; 2,632 rows). The eLasso model against IsingFit::IsingFit(); the unregularized
# nodewise model against IsingSampler::EstimateIsing(method = "uni") and against stats::glm() (the
# definition); the lasso KKT residual against its conditions written out at a glmnet solution; the split
# of ordinal items against its definition.

binary <- utils::read.csv(reference_data("data/bfi6_binary.csv"))
items <- utils::read.csv(reference_data("data/bfi10.csv"))
nodes <- names(binary)

isingfit <- function(and) IsingFit::IsingFit(binary, gamma = 0.25, AND = and, plot = FALSE, progressbar = FALSE)

test_that("native = FALSE: the eLasso network and thresholds equal IsingFit::IsingFit()", {
  skip_without_reference(c("IsingFit", "glmnet"))
  and <- isingfit(TRUE)
  or <- isingfit(FALSE)
  # Both run glmnet's path for each item and pick it by the EBIC; psychnets calls glmnet itself, IsingFit
  # calls glmnet: the selection, the AND / OR rule and the averaging are the two packages' own code.
  # glmnet's convergence threshold (1e-7) is shared, so only rounding may differ.
  expect_reference(as.data.frame(ising_fit(binary, rule = "AND", native = FALSE), include_zero = TRUE),
                   ref_edges(and$weiadj), reference = "IsingFit", covers = "ising_fit()",
                   quantity = "edge weight, AND rule, native = FALSE", tolerance = 1e-8, key = c("from", "to"),
                   note = "psychnets (native = FALSE) and IsingFit both fit the nodewise lasso with glmnet.")
  expect_reference(as.data.frame(ising_fit(binary, rule = "OR", native = FALSE), include_zero = TRUE),
                   ref_edges(or$weiadj), reference = "IsingFit", covers = "ising_fit()",
                   quantity = "edge weight, OR rule, native = FALSE", tolerance = 1e-8, key = c("from", "to"))
  expect_reference(ising_fit(binary, rule = "AND", native = FALSE)$thresholds, stats::setNames(and$thresholds, nodes),
                   reference = "IsingFit", covers = c("ising_fit()", "$.psychnet()"),
                   quantity = "threshold, AND rule, native = FALSE", tolerance = 1e-8)
})

test_that("native = TRUE: the eLasso network equals IsingFit::IsingFit() on glmnet's path", {
  skip_without_reference("IsingFit")
  and <- isingfit(TRUE)
  # The native solver selects on glmnet's default path (lambda_min_ratio 1e-4 here, glmnet's early stop), as
  # IsingFit does, and solves each penalty to its KKT conditions; IsingFit's glmnet stops at thresh = 1e-7, which
  # leaves errors of order 1e-5 in an edge weight: 1e-4 bounds them.
  expect_reference(as.data.frame(ising_fit(binary, rule = "AND"), include_zero = TRUE), ref_edges(and$weiadj),
                   reference = "IsingFit", covers = "ising_fit()", quantity = "edge weight, AND rule, native solver",
                   tolerance = 1e-4, key = c("from", "to"))
})

test_that("the unregularized nodewise model equals IsingSampler::EstimateIsing(method = \"uni\")", {
  skip_without_reference("IsingSampler")
  reference <- IsingSampler::EstimateIsing(as.matrix(binary), responses = c(0L, 1L), method = "uni")
  graph <- reference$graph
  dimnames(graph) <- list(nodes, nodes)
  # IsingSampler fits each regression with stats::glm() at its default epsilon (1e-8 on the deviance):
  # coefficients of order one agree to about 1e-6.
  expect_reference(as.data.frame(ising_sampler(binary, rule = "AND", alpha = NULL), include_zero = TRUE),
                   ref_edges(graph), reference = "IsingSampler", covers = "ising_sampler()",
                   quantity = "edge weight, no pruning", tolerance = 1e-6, key = c("from", "to"))
  expect_reference(ising_sampler(binary, rule = "AND", alpha = NULL)$thresholds,
                   stats::setNames(reference$thresholds, nodes), reference = "IsingSampler",
                   covers = "ising_sampler()", quantity = "threshold, no pruning", tolerance = 1e-6)
})

test_that("the unregularized nodewise model equals logistic regressions by stats::glm()", {
  answers <- as.matrix(binary)
  # one column per regressed item: its intercept, then its coefficient on every item (0 on itself)
  nodewise <- vapply(nodes, \(item) {
    others <- setdiff(nodes, item)
    fit <- stats::glm(answers[, item] ~ answers[, others], family = stats::binomial(),
                      control = stats::glm.control(epsilon = 1e-12, maxit = 100))
    replace(stats::setNames(numeric(length(nodes) + 1L), c("threshold", nodes)), c("threshold", others), stats::coef(fit))
  }, numeric(length(nodes) + 1L))
  weights <- (nodewise[nodes, ] + t(nodewise[nodes, ])) / 2
  fit <- ising_sampler(binary, rule = "AND", alpha = NULL)
  by_glm <- "independent: nodewise logistic regression by stats::glm(), edges the mean of the two coefficients"
  # Two IRLS solvers of six concave problems: glm() to epsilon = 1e-12, psychnets certified to about 1e-8.
  expect_reference(as.data.frame(fit, include_zero = TRUE), ref_edges(weights), reference = by_glm,
                   covers = c("ising_sampler()", "as.data.frame.psychnet()"), quantity = "edge weight, no pruning",
                   tolerance = 1e-6, key = c("from", "to"))
  expect_reference(fit$thresholds, nodewise["threshold", ], reference = by_glm, covers = "ising_sampler()",
                   quantity = "threshold, no pruning", tolerance = 1e-6)
  expect_reference(as.data.frame(psychnet(binary, method = "ising_sampler"), include_zero = TRUE), ref_edges(weights),
                   reference = by_glm, covers = "psychnet()", quantity = "edge weight, psychnet(method = \"ising_sampler\")",
                   tolerance = 1e-6, key = c("from", "to"))
})

test_that("glm_lasso_kkt() equals the lasso KKT conditions at a glmnet solution", {
  skip_without_reference("glmnet")
  X <- as.matrix(items[c("A2", "A3", "A4", "A5", "C1")])
  y <- items$A1
  fit <- glmnet::glmnet(X, y, family = "gaussian", lambda = 0.05, standardize = FALSE)
  beta <- as.numeric(fit$beta)
  b0 <- as.numeric(fit$a0)
  # Friedman, Hastie and Tibshirani (2010): with the objective RSS / (2n) + lambda |beta|_1, the gradient
  # g = X'(y - mu) / n must equal lambda sign(beta_j) on the active set, |g_j| <= lambda elsewhere, and the
  # intercept's score sum(y - mu) / n must be 0. The residual is the largest violation.
  mu <- b0 + as.numeric(X %*% beta)
  gradient <- as.numeric(crossprod(X, y - mu)) / nrow(X)
  active <- abs(beta) > 1e-8
  residual <- max(abs(sum(y - mu) / nrow(X)), abs(gradient[active] - 0.05 * sign(beta[active])),
                  pmax(abs(gradient[!active]) - 0.05, 0))
  expect_reference(glm_lasso_kkt(X, y, b0 = b0, beta = beta, lambda = 0.05), residual,
                   reference = "independent: the lasso KKT conditions at glmnet's gaussian solution, base R",
                   covers = "glm_lasso_kkt()", quantity = "lasso KKT residual, gaussian", tolerance = 1e-12)
  z <- binary$A1
  logistic <- glmnet::glmnet(as.matrix(binary[-1L]), z, family = "binomial", lambda = 0.02, standardize = FALSE)
  beta_l <- as.numeric(logistic$beta)
  b0_l <- as.numeric(logistic$a0)
  eta <- b0_l + as.numeric(as.matrix(binary[-1L]) %*% beta_l)
  mu_l <- 1 / (1 + exp(-eta))
  gradient_l <- as.numeric(crossprod(as.matrix(binary[-1L]), z - mu_l)) / length(z)
  active_l <- abs(beta_l) > 1e-8
  residual_l <- max(abs(sum(z - mu_l) / length(z)), abs(gradient_l[active_l] - 0.02 * sign(beta_l[active_l])),
                    pmax(abs(gradient_l[!active_l]) - 0.02, 0))
  expect_reference(glm_lasso_kkt(as.matrix(binary[-1L]), z, b0 = b0_l, beta = beta_l, lambda = 0.02, family = "binomial"),
                   residual_l, reference = "independent: the lasso KKT conditions at glmnet's gaussian solution, base R",
                   note = "The same conditions at glmnet's binomial (logistic) solution.",
                   covers = "glm_lasso_kkt()", quantity = "lasso KKT residual, binomial", tolerance = 1e-12)
})

test_that("dichotomize() equals its split rules written out", {
  rule <- "independent: the split rule of each column, base R"
  items <- utils::head(items, 300L)
  at_median <- vapply(items, \(x) as.integer(x >= stats::median(x)), integer(nrow(items)))
  above_mean <- vapply(items, \(x) as.integer(x > mean(x)), integer(nrow(items)))
  above_mid_rank <- vapply(items, \(x) as.integer(rank(x) > (length(x) + 1) / 2), integer(nrow(items)))
  expect_reference(dichotomize(items, method = "median"), at_median, reference = rule, covers = "dichotomize()",
                   quantity = "1 at or above the median", tolerance = 0)
  expect_reference(dichotomize(items, method = "mean"), above_mean, reference = rule, covers = "dichotomize()",
                   quantity = "1 above the mean", tolerance = 0)
  expect_reference(dichotomize(items, method = "rank"), above_mid_rank, reference = rule, covers = "dichotomize()",
                   quantity = "1 above the middle rank", tolerance = 0)
})

test_that("dichotomize() at the median and the mean equals bootnet::binarize()", {
  skip_without_reference("bootnet")
  items <- utils::head(items, 300L)
  expect_reference(dichotomize(items, method = "median"), as.matrix(bootnet::binarize(items, split = "median", verbose = FALSE)),
                   reference = "bootnet", covers = "dichotomize()", quantity = "1 at or above the median, binarize()",
                   tolerance = 0)
  expect_reference(dichotomize(items, method = "mean"), as.matrix(bootnet::binarize(items, split = "mean", verbose = FALSE)),
                   reference = "bootnet", covers = "dichotomize()", quantity = "1 at or above the mean, binarize()",
                   tolerance = 0)
})
