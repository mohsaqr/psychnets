# panel_gvar(): the panel graphical VAR. Formula tests (implied moments, both
# gradient routes), equivalence with psychonetrics::panelgvar() via a stored
# reference fixture, invariance properties, pruning, error/warning classes, and
# the accessors. Data come from helper-panelgvar.R.

vars3 <- c("V1", "V2", "V3")
# Muffle the boundary warning where a test deliberately fits a boundary-prone
# (small-n or saturated) model and checks something else.
quiet_boundary <- function(expr) withCallingHandlers(expr,
  psychnet_boundary = function(w) invokeRestart("muffleWarning"),
  psychnet_singular_information = function(w) invokeRestart("muffleWarning"),
  psychnet_approximate_se = function(w) invokeRestart("muffleWarning"))
# Muffle only the identification warnings (to assert a boundary warning).
quiet_identification <- function(expr) withCallingHandlers(expr,
  psychnet_singular_information = function(w) invokeRestart("muffleWarning"),
  psychnet_approximate_se = function(w) invokeRestart("muffleWarning"))
fit_long <- function(long, ...) {
  panel_gvar(long, vars = setdiff(names(long), c("id", "wave")),
             id = "id", wave = "wave", ...)
}
# Rebuild the internal context the estimator works in, from wide data.
pgv_ctx <- function(w) {
  wd <- psychnets:::.pgv_wide(w$data, w$design, NULL, NULL)
  y <- wd$wide
  means <- colMeans(y)
  lay <- psychnets:::.pgv_layout(length(wd$labels))
  list(lay = lay, waves = length(wd$waves),
       cov = crossprod(sweep(y, 2L, means)) / nrow(y), means = means,
       n = nrow(y), free = rep(TRUE, lay$n_par))
}

test_that("implied moments are the stationary block-Toeplitz structure", {
  w <- panel_sim_wide()
  fit <- panel_gvar(w$data, vars = w$design, select = "none")
  ctx <- pgv_ctx(w)
  st <- psychnets:::.pgv_state(as.data.frame(fit, what = "parameters")$estimate,
                               ctx)
  b <- st$par$beta
  # Stationarity: Sigma0 = Sigma_zeta + B Sigma0 B'.
  expect_equal(st$s0, st$zw$sigma + b %*% st$s0 %*% t(b), tolerance = 1e-12)
  # Block (3, 1) = B^2 Sigma0 + Sigma_between; (1, 3) is its transpose.
  blk <- function(i, j) st$sigma[(i - 1) * 3 + 1:3, (j - 1) * 3 + 1:3]
  expect_equal(blk(3, 1), b %*% b %*% st$s0 + st$zb$sigma, tolerance = 1e-12)
  expect_equal(blk(1, 3), t(blk(3, 1)), tolerance = 1e-12)
  # GGM parameterization: Sigma = D (I - Omega)^-1 D.
  expect_equal(st$zw$sigma,
               diag(st$par$delta_w) %*% solve(diag(3) - st$par$omega_w) %*%
                 diag(st$par$delta_w), tolerance = 1e-12)
})

test_that("reverse-mode gradient equals forward-mode score and finite differences", {
  ctx <- pgv_ctx(panel_sim_wide(n = 200, p = 4, waves = 4, seed = 3))
  theta <- psychnets:::.pgv_start(ctx)
  st <- psychnets:::.pgv_state(theta, ctx)
  g_rev <- psychnets:::.pgv_gradient_full(st, ctx)
  jac <- psychnets:::.pgv_jacobian(st, ctx, seq_len(ctx$lay$n_par))
  g_fwd <- psychnets:::.pgv_score(st, ctx, jac)
  f <- function(th) psychnets:::.pgv_state(th, ctx)$f
  g_num <- vapply(seq_along(theta), function(i) {
    h <- 1e-6
    e <- replace(numeric(length(theta)), i, h)
    (f(theta + e) - f(theta - e)) / (2 * h)
  }, numeric(1))
  expect_equal(g_rev, g_fwd, tolerance = 1e-10)
  expect_lt(max(abs(g_rev - g_num)), 1e-6)
  expect_gt(max(abs(g_rev)), 1e-3)   # a real gradient, not a trivial zero
})

test_that("the fit is certified: ML score ~ 0 at the returned estimate", {
  fit <- fit_long(panel_sim(), select = "none")
  expect_lt(attr(fit, "kkt"), 1e-8)
  pruned <- fit_long(panel_sim(), select = "prune")
  expect_lt(attr(pruned, "kkt"), 1e-8)
})

test_that("matches psychonetrics::panelgvar() (stored reference fixture)", {
  ref <- readRDS(test_path("fixtures", "psychonetrics-panelgvar.rds"))
  cases <- list(default = list(n = 400, p = 3, waves = 4, seed = 1,
                               alpha = 0.01, adjust = "none"),
                large = list(n = 600, p = 5, waves = 5, seed = 8,
                             alpha = 0.01, adjust = "none"),
                bh = list(n = 150, p = 4, waves = 4, seed = 10,
                          alpha = 0.05, adjust = "BH"),
                stepup_small = list(n = 45, p = 4, waves = 4, seed = 21,
                                    alpha = 0.01, adjust = "none"),
                stepup_waves = list(n = 300, p = 4, waves = 5, seed = 12,
                                    alpha = 0.01, adjust = "none"))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    rf <- ref[[nm]]
    w <- panel_sim_wide(n = cs$n, p = cs$p, waves = cs$waves, seed = cs$seed)
    sat <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "none"))
    ps <- as.data.frame(sat, what = "parameters")
    ctx <- pgv_ctx(w)
    f_own <- psychnets:::.pgv_state(ps$estimate, ctx)$f
    f_ref <- psychnets:::.pgv_state(rf$estimate, ctx)$f
    expect_equal(f_ref, rf$objective, tolerance = 1e-10)  # same objective
    expect_lte(f_own, f_ref + 1e-12)
    # The two step-up designs are weakly identified (a boundary at n = 45, a
    # flat ridge over five waves: psychonetrics sits ~1e-5 higher in F with
    # estimates 0.1 away), so only the objective and the selected supports are
    # comparable there; estimate/SE/fit-index agreement is checked on the
    # well-identified designs.
    if (startsWith(nm, "stepup")) {
      pru <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                        adjust = cs$adjust, stepup = FALSE))
      expect_identical(as.data.frame(pru, what = "parameters")$fixed,
                       rf$pruned_fixed, label = nm)
      up <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                       adjust = cs$adjust))
      expect_identical(as.data.frame(up, what = "parameters")$fixed,
                       rf$stepup_fixed, label = nm)
      next
    }
    # Saturated estimates and SEs agree to psychonetrics' optimizer tolerance:
    # its default nlminb stop leaves a gradient of ~1e-3 (measured 6.8e-4 on
    # the default case), which moves estimates by ~1e-3 on the 5 x 5 and
    # n = 150 designs. The decisive checks are the next ones: ours is the
    # deeper optimum of the identical objective, with the identical support.
    expect_lt(max(abs(ps$estimate - rf$estimate)), 5e-3)
    expect_lt(max(abs(ps$se / rf$se - 1)), 5e-3)
    ft <- as.data.frame(sat, what = "fit")
    expect_identical(ft$df, rf$fit[["df"]])
    expect_lt(abs(ft$chisq - rf$fit[["chisq"]]), 1e-3)
    expect_lt(abs(ft$baseline_chisq - rf$fit[["baseline.chisq"]]), 1e-3)
    expect_lt(abs(ft$cfi - rf$fit[["cfi"]]), 1e-5)
    expect_lt(abs(ft$tli - rf$fit[["tli"]]), 1e-5)
    expect_lt(abs(ft$rmsea - rf$fit[["rmsea"]]), 1e-5)
    expect_lt(abs(ft$bic - rf$fit[["bic"]]), 1e-3)
    # Pruning alone (stepup = FALSE) selects psychonetrics::prune()'s support,
    # and the default prune + stepup selects prune() |> stepup()'s support.
    pru <- panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                      adjust = cs$adjust, stepup = FALSE)
    pp <- as.data.frame(pru, what = "parameters")
    expect_identical(pp$fixed, rf$pruned_fixed, label = nm)
    expect_lt(max(abs(pp$estimate - rf$pruned_estimate)), 5e-3)
    up <- panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                     adjust = cs$adjust)
    pu <- as.data.frame(up, what = "parameters")
    expect_identical(pu$fixed, rf$stepup_fixed, label = nm)
    expect_lt(max(abs(pu$estimate - rf$stepup_estimate)), 5e-3)
  }
})

test_that("BH adjustment prunes more than no adjustment on the reference case", {
  long <- panel_sim(n = 150, p = 4, waves = 4, seed = 10)
  none <- fit_long(long, select = "prune", alpha = 0.05)
  bh <- fit_long(long, select = "prune", alpha = 0.05, adjust = "BH")
  expect_gt(nrow(as.data.frame(none)), nrow(as.data.frame(bh)))
})

test_that("long and wide input give the same fit", {
  long <- panel_sim()
  w <- panel_sim_wide()
  a <- fit_long(long)
  b <- panel_gvar(w$data, vars = w$design)
  expect_equal(as.data.frame(a, what = "parameters"),
               as.data.frame(b, what = "parameters"), tolerance = 1e-10)
  # Wave labels are only an order: rescaling them changes nothing.
  long$wave <- long$wave * 10
  c10 <- fit_long(long)
  expect_equal(as.data.frame(c10, what = "parameters")$estimate,
               as.data.frame(a, what = "parameters")$estimate, tolerance = 1e-10)
})

test_that("variable order is a relabelling (permutation invariance)", {
  long <- panel_sim(n = 300, p = 4, waves = 4, seed = 2)
  a <- fit_long(long, select = "none")
  b <- panel_gvar(long, vars = c("V3", "V1", "V4", "V2"), id = "id",
                  wave = "wave", select = "none")
  for (net in c("temporal", "contemporaneous", "between")) {
    wa <- a[[net]]$weights
    wb <- b[[net]]$weights
    expect_equal(wb, wa[rownames(wb), colnames(wb)], tolerance = 1e-7)
  }
  expect_equal(b$temporal$beta, a$temporal$beta[rownames(b$temporal$beta),
                                                colnames(b$temporal$beta)],
               tolerance = 1e-7)
})

test_that("rescaling a variable leaves partial correlations and rescales beta", {
  long <- panel_sim(n = 300, p = 3, waves = 4, seed = 4)
  a <- fit_long(long, select = "none")
  long$V2 <- 3 * long$V2
  b <- fit_long(long, select = "none")
  expect_equal(b$contemporaneous$weights, a$contemporaneous$weights,
               tolerance = 1e-7)
  expect_equal(b$between$weights, a$between$weights, tolerance = 1e-7)
  scale <- c(1, 3, 1)
  expect_equal(b$temporal$beta, a$temporal$beta * outer(scale, 1 / scale),
               tolerance = 1e-7)
})

test_that("pruning fixes non-significant edges to zero; recursive leaves none", {
  long <- panel_sim(n = 200, p = 4, waves = 4, seed = 5)
  pru <- fit_long(long, select = "prune", stepup = FALSE)
  pp <- as.data.frame(pru, what = "parameters")
  fixed <- pp[pp$fixed, ]
  expect_gt(nrow(fixed), 0)
  expect_true(all(fixed$estimate == 0))
  expect_true(all(fixed$matrix %in%
                    c("beta", "omega_zeta_within", "omega_zeta_between")))
  expect_true(all(is.na(fixed$se)))
  rec <- fit_long(long, select = "prune", recursive = TRUE, stepup = FALSE)
  edges <- as.data.frame(rec)
  expect_true(all(edges$p <= 0.01))
  # The temporal network carries the pruned beta, transposed to from -> to.
  expect_equal(pru$temporal$weights, t(pru$temporal$beta) * (1 - diag(4)),
               ignore_attr = TRUE)
})

test_that("input errors are classed", {
  long <- panel_sim(n = 100, p = 3, waves = 4)
  expect_error(fit_long(long[long$wave <= 2, ]), class = "psychnet_bad_input")
  expect_error(fit_long(long, na_method = "pairwise"),
               class = "psychnet_bad_input")
  expect_error(panel_gvar(long, vars = c("V1", "nope"), id = "id",
                          wave = "wave"), class = "psychnet_bad_input")
  expect_error(panel_gvar(long, vars = vars3), class = "psychnet_bad_input")
  expect_error(fit_long(rbind(long, long[1, ])), class = "psychnet_bad_input")
  bad <- long
  bad$V1 <- as.character(bad$V1)
  expect_error(fit_long(bad), class = "psychnet_bad_input")
  expect_error(fit_long(long, labels = c("a", "b")),
               class = "psychnet_bad_input")
  expect_error(fit_long(long[long$id <= 10, ]), class = "psychnet_bad_input")
  expect_error(fit_long(long, alpha = 2))
})

test_that("missing values drop persons listwise with a classed warning", {
  long <- panel_sim(n = 300)
  long$V2[long$id == 7 & long$wave == 3] <- NA
  expect_warning(fit <- fit_long(long), class = "psychnet_listwise_drop")
  expect_identical(attr(fit, "n"), 299L)
  # A person missing a whole wave row is dropped the same way.
  expect_warning(fit_long(long[!(long$id == 9 & long$wave == 2), ]),
                 class = "psychnet_listwise_drop")
})

test_that("a near-singular between covariance raises a boundary warning", {
  # Three waves, saturated: the between partial correlations go to +/-1, the
  # same boundary optimum psychonetrics reaches (with SEs > 10).
  expect_warning(fit <- quiet_identification(
    fit_long(panel_sim(n = 500, p = 3, waves = 3, seed = 7), select = "none")),
    class = "psychnet_boundary")
  dg <- as.data.frame(fit, what = "diagnostics")
  expect_true(dg$boundary)
  expect_match(dg$boundary_part, "between-person covariance near-singular")
})

test_that("accessors return tidy data frames", {
  fit <- fit_long(panel_sim())
  expect_s3_class(fit, c("psychnet_panelgvar", "netobject_group"))
  expect_named(fit, c("temporal", "contemporaneous", "between"))
  expect_true(fit$temporal$directed)
  expect_false(fit$between$directed)
  edges <- as.data.frame(fit)
  expect_named(edges, c("network", "from", "to", "weight", "se", "p",
                        "ci_lower", "ci_upper", "significant"))
  # Default: the full model, every edge reported; significance flagged at alpha.
  expect_identical(nrow(edges), 9L + 3L + 3L)
  expect_identical(edges$significant, !is.na(edges$p) & edges$p < 0.01)
  expect_identical(summary(fit)$significant,
                   vapply(c("temporal", "contemporaneous", "between"), function(nm)
                     sum(edges$significant[edges$network == nm & edges$from != edges$to]),
                     integer(1), USE.NAMES = FALSE))
  expect_setequal(unique(edges$network),
                  c("temporal", "contemporaneous", "between"))
  expect_true(all(edges$ci_lower < edges$weight & edges$weight < edges$ci_upper))
  pruned_fit <- fit_long(panel_sim(), select = "prune")
  expect_lt(nrow(as.data.frame(pruned_fit)), 15L)
  expect_identical(nrow(as.data.frame(pruned_fit, include_zero = TRUE)), 15L)
  params <- as.data.frame(fit, what = "parameters")
  expect_identical(nrow(params), 3L + 3L + 3L + 9L + 3L + 3L)
  ft <- as.data.frame(fit, what = "fit")
  expect_identical(nrow(ft), 1L)
  expect_named(ft, c("logl", "npar", "df", "chisq", "pvalue", "baseline_chisq",
                     "baseline_df", "cfi", "tli", "rmsea", "aic", "bic"))
  sm <- summary(fit)
  expect_identical(sm$network, c("temporal", "contemporaneous", "between"))
  dg <- as.data.frame(fit, what = "diagnostics")
  expect_identical(nrow(dg), 1L)
  expect_named(dg, c("n", "waves", "score", "boundary", "boundary_part",
                     "se_available", "starts_agree", "starts"))
  expect_identical(dg$n, 400L)
  expect_false(dg$boundary)
  expect_true(dg$se_available)
  expect_identical(dg$starts_agree, 3L)
  # Variance split: within + between reproduce the implied total variance.
  vr <- as.data.frame(fit, what = "variance")
  expect_named(vr, c("variable", "between", "within", "share_between"))
  ctx <- pgv_ctx(panel_sim_wide())
  st <- psychnets:::.pgv_state(as.data.frame(fit, what = "parameters")$estimate,
                               ctx)
  expect_equal(vr$between + vr$within, unname(diag(st$sigma))[1:3],
               tolerance = 1e-10)
  expect_true(all(vr$share_between > 0 & vr$share_between < 1))
  expect_snapshot(print(fit), transform = function(x)
    sub("score\\): .*$", "score): <certificate>", x))
})

test_that("plots the three networks, with autoregressions as self-loops", {
  skip_if_not_installed("cograph")
  fit <- fit_long(panel_sim())
  drawn <- psychnets:::.pgv_plot_temporal(fit$temporal)
  loops <- subset(drawn$edges, from == to)
  ar <- diag(fit$temporal$beta)
  expect_identical(nrow(loops), sum(abs(ar) > 1e-12))
  expect_equal(loops$weight, ar[abs(ar) > 1e-12], ignore_attr = TRUE)
  expect_equal(unname(diag(drawn$weights)), unname(ar))
  # The stored network keeps the zero-diagonal psychnet contract.
  expect_true(all(diag(fit$temporal$weights) == 0))
  # Solid/dashed follow the edge table's significance, edge by edge.
  edges <- as.data.frame(fit)
  sig_t <- psychnets:::.pgv_edge_significance(drawn, "temporal", edges)
  lab <- drawn$nodes$label
  expect_identical(sig_t, vapply(seq_len(nrow(drawn$edges)), function(i)
    subset(edges, network == "temporal" & from == lab[drawn$edges$from[i]] &
             to == lab[drawn$edges$to[i]])$significant, logical(1)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit), fit)
  expect_identical(plot(fit, node_fill = "red", layout = "circle",
                        title = "Custom panel", edge_positive_color = "green",
                        edge_label_size = 0.8), fit)
  expect_identical(plot(fit_long(panel_sim(), select = "none"),
                        edge_labels = FALSE), fit_long(panel_sim(), select = "none"))
})

# Regression: raw 0-100-scale data with few persons (the weekly-mean ESM case)
# sent the single-start, unscaled optimizer to F ~ 1e5 against psychonetrics'
# ~98. The fit now runs on standardized variables from three starts.
test_that("raw-scale, small-n data fit as well as their standardized version", {
  long <- panel_sim(n = 45, p = 4, waves = 4, seed = 21)
  vars <- c("V1", "V2", "V3", "V4")
  raw <- long
  raw[vars] <- 50 + 18 * long[vars]
  a <- quiet_boundary(fit_long(long, select = "none"))
  b <- quiet_boundary(fit_long(raw, select = "none"))
  expect_true(all(attr(b, "starts")$best))
  expect_equal(b$contemporaneous$weights, a$contemporaneous$weights,
               tolerance = 1e-6)
  # This small-n fit is a boundary solution: the between partial correlations
  # lie in a nearly flat direction of the likelihood, so the two scales agree
  # there to ~1e-5 (measured 6.3e-6) while the objective is identical.
  expect_equal(b$between$weights, a$between$weights, tolerance = 1e-4)
  expect_equal(b$temporal$beta, a$temporal$beta, tolerance = 1e-6)
  # Fit indices are scale-invariant; the model never ends above its baseline.
  fa <- as.data.frame(a, what = "fit")
  fb <- as.data.frame(b, what = "fit")
  expect_equal(fb$chisq, fa$chisq, tolerance = 1e-6)
  expect_lt(fb$chisq, fb$baseline_chisq)
})

test_that("modification indices match psychonetrics' at the pruned model", {
  ref <- readRDS(test_path("fixtures", "psychonetrics-panelgvar.rds"))
  for (nm in c("default", "bh", "stepup_small")) {
    rf <- ref[[nm]]
    cs <- switch(nm, default = list(400, 3, 4, 1, 0.01, "none"),
                 bh = list(150, 4, 4, 10, 0.05, "BH"),
                 stepup_small = list(45, 4, 4, 21, 0.01, "none"))
    w <- panel_sim_wide(n = cs[[1]], p = cs[[2]], waves = cs[[3]],
                        seed = cs[[4]])
    pru <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "prune",
                                     alpha = cs[[5]], adjust = cs[[6]],
                                     stepup = FALSE))
    pp <- as.data.frame(pru, what = "parameters")
    ctx <- pgv_ctx(w)
    ctx$free <- !pp$fixed
    st <- psychnets:::.pgv_state(pp$estimate, ctx)
    cand <- which(pp$fixed)
    mod <- psychnets:::.pgv_modification(st, ctx, cand)
    # psychonetrics evaluates the score test at its own (looser) optimum, so
    # agreement is to that tolerance; the decision-relevant ranking is exact.
    expect_equal(mod$mi, rf$pruned_mi[cand], tolerance = 2e-2, label = nm)
    expect_identical(order(-mod$mi), order(-rf$pruned_mi[cand]), label = nm)
  }
})

test_that("an MI equals the likelihood-ratio statistic to first order", {
  w <- panel_sim_wide(n = 45, p = 4, waves = 4, seed = 21)
  pru <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "prune",
                                   stepup = FALSE))
  pp <- as.data.frame(pru, what = "parameters")
  ctx <- pgv_ctx(w)
  ctx$free <- !pp$fixed
  st <- psychnets:::.pgv_state(pp$estimate, ctx)
  mod <- psychnets:::.pgv_modification(st, ctx, which(pp$fixed))
  best <- mod$index[which.max(mod$mi)]
  up <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "prune"))
  lr <- as.data.frame(pru, what = "fit")$chisq - as.data.frame(up, what = "fit")$chisq
  # The step-up added exactly this parameter; its MI predicts the LR gain.
  expect_identical(which(xor(pp$fixed, as.data.frame(up, what = "parameters")$fixed)),
                   best)
  expect_equal(max(mod$mi), lr, tolerance = 0.25)
})

test_that("stepup re-adds a parameter one-pass pruning removed", {
  long <- panel_sim(n = 45, p = 4, waves = 4, seed = 21)
  pruned <- quiet_boundary(fit_long(long, select = "prune", stepup = FALSE))
  stepped <- quiet_boundary(fit_long(long, select = "prune"))
  expect_identical(attr(stepped, "settings")$added, 1L)
  expect_identical(nrow(as.data.frame(stepped)), nrow(as.data.frame(pruned)) + 1L)
  # The re-added model fits significantly better and has the lower BIC.
  fp <- as.data.frame(pruned, what = "fit")
  fs <- as.data.frame(stepped, what = "fit")
  expect_lt(fs$bic, fp$bic)
  expect_lt(stats::pchisq(fp$chisq - fs$chisq, 1, lower.tail = FALSE), 0.01)
  # criterion = "none" can only add at least as many parameters.
  none <- quiet_boundary(fit_long(long, select = "prune", criterion = "none"))
  expect_gte(attr(none, "settings")$added, 1L)
  # stepup is skipped when nothing was pruned.
  expect_identical(attr(quiet_boundary(fit_long(long, select = "none")),
                       "settings")$added, 0L)
})

test_that("unidentified parameters get infinite SEs for pruning", {
  set.seed(5)
  a <- crossprod(matrix(rnorm(40), 10, 4))
  h <- matrix(0, 5, 5)
  h[1:4, 1:4] <- a                       # parameter 5 is not identified at all
  se <- psychnets:::.pgv_identified_se(h)
  expect_identical(se[5], Inf)
  expect_equal(se[1:4], sqrt(diag(solve(a))), tolerance = 1e-10)
  # An invertible matrix is inverted exactly.
  expect_equal(psychnets:::.pgv_identified_se(a), sqrt(diag(solve(a))),
               tolerance = 1e-10)
})

test_that("a variance driven to zero is reported as a boundary, by name", {
  # V2 has no stable trait, so its between-person variance heads to zero.
  long <- panel_sim(n = 300, p = 3, waves = 4, seed = 2,
                    trait_scale = c(1, 0, 1))
  expect_warning(fit <- quiet_identification(fit_long(long, select = "none")),
                 class = "psychnet_boundary")
  dg <- as.data.frame(fit, what = "diagnostics")
  expect_true(dg$boundary)
  expect_match(dg$boundary_part, "between-person variance of V2")
})

# ---- select = "ebic": penalized ML with EBIC (psychonetrics PML) -------------

test_that("the EBIC fit satisfies the lasso KKT conditions on the data scale", {
  skip_on_cran()
  w <- panel_sim_wide()
  fit <- panel_gvar(w$data, vars = w$design, select = "ebic")
  pp <- as.data.frame(fit, what = "parameters")
  lambda <- attr(fit, "settings")$lambda
  ctx <- pgv_ctx(w)
  st <- psychnets:::.pgv_state(pp$estimate, ctx)
  g <- psychnets:::.pgv_gradient_full(st, ctx)
  pen <- pp$matrix %in% c("beta", "omega_zeta_within", "omega_zeta_between")
  nonzero <- pen & pp$estimate != 0
  zero <- pen & pp$estimate == 0
  # Recomputed here, independently of the package's certificate function.
  expect_lt(max(abs(g[!pen])), 1e-8)
  expect_lt(max(abs(g[nonzero] + lambda * sign(pp$estimate[nonzero]))), 1e-8)
  expect_true(all(abs(g[zero]) <= lambda + 1e-8))
  expect_gt(sum(zero), 0)                      # the lasso really zeroes edges
  expect_lt(attr(fit, "kkt"), 1e-8)
  # The path: lambda_max empties every network; the choice is the EBIC minimum.
  path <- attr(fit, "settings")$path
  expect_identical(path$edges[1], 0)
  expect_identical(lambda, path$lambda[which.min(path$ebic)])
  expect_true(all(attr(fit, "starts")$best))
})

test_that("select = 'ebic' matches psychonetrics' PML (stored reference)", {
  skip_on_cran()
  ref <- readRDS(test_path("fixtures", "psychonetrics-panelgvar.rds"))$pml
  designs <- list(default = list(400, 3, 4, 1), bh = list(150, 4, 4, 10),
                  stepup_waves = list(300, 4, 5, 12))
  for (nm in names(designs)) {
    d <- designs[[nm]]
    rf <- ref[[nm]]
    w <- panel_sim_wide(n = d[[1]], p = d[[2]], waves = d[[3]], seed = d[[4]])
    fit <- panel_gvar(w$data, vars = w$design, select = "ebic")
    pp <- as.data.frame(fit, what = "parameters")
    pen <- pp$matrix %in% c("beta", "omega_zeta_within", "omega_zeta_between")
    # Same selected lambda (the grids differ only through lambda_max, i.e. the
    # precision of each package's null-model fit), same zero pattern, and
    # estimates to psychonetrics' optimizer tolerance.
    expect_equal(attr(fit, "settings")$lambda, rf$lambda, tolerance = 1e-4,
                 label = nm)
    expect_identical(pen & pp$estimate == 0, pen & abs(rf$estimate) < 1e-5,
                     label = nm)
    expect_lt(max(abs(pp$estimate - rf$estimate)), 5e-3)
    # Ours is at least as good an optimum of the same penalized objective.
    ctx <- pgv_ctx(w)
    objective <- function(th) psychnets:::.pgv_state(th, ctx)$f +
      rf$lambda * sum(abs(th[pen]))
    expect_lte(objective(pp$estimate), objective(rf$estimate) + 1e-10)
  }
})

test_that("EBIC estimates carry no p-values; refit = TRUE supplies them", {
  skip_on_cran()
  long <- panel_sim()
  fit <- fit_long(long, select = "ebic")
  edges <- as.data.frame(fit)
  expect_true(all(is.na(edges$p)) && all(is.na(edges$significant)))
  expect_true(all(is.na(summary(fit)$significant)))
  rf <- fit_long(long, select = "ebic", refit = TRUE)
  er <- as.data.frame(rf)
  expect_identical(subset(er, select = c(network, from, to)),
                   subset(edges, select = c(network, from, to)))
  expect_false(anyNA(er$p))
  expect_lt(attr(rf, "kkt"), 1e-8)            # an unpenalized ML optimum
  expect_true(attr(rf, "settings")$refit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit), fit)
  expect_snapshot(print(fit), transform = function(x)
    sub("residual\\): .*$", "residual): <certificate>", x))
})

test_that("the EBIC engine stays inside the model on a hard small-n design", {
  skip_on_cran()
  # n = 45 over four waves sends the unpenalized fit to a boundary; the
  # penalized fit must still be a valid model with a certified KKT residual.
  long <- panel_sim(n = 45, p = 4, waves = 4, seed = 21)
  fit <- quiet_boundary(fit_long(long, select = "ebic"))
  pp <- as.data.frame(fit, what = "parameters")
  w <- panel_sim_wide(n = 45, p = 4, waves = 4, seed = 21)
  expect_false(is.null(psychnets:::.pgv_state(pp$estimate, pgv_ctx(w))))
  expect_lt(attr(fit, "kkt"), 1e-6)
  expect_gte(sum(attr(fit, "starts")$best), 1L)
})

# ---- regressions from the Codex adversarial review (2026-10-03) --------------

test_that("the returned EBIC model is the EBIC minimum after final solves", {
  skip_on_cran()
  # A final, certified solve can find a denser optimum than the short path fit
  # and raise that lambda's EBIC; the choice must then move.
  w <- panel_sim_wide(n = 100, p = 3, waves = 4, seed = 8,
                      trait_scale = c(1, 0, 1))
  fit <- quiet_boundary(panel_gvar(w$data, vars = w$design, select = "ebic",
                                   nlambda = 12))
  path <- attr(fit, "settings")$path
  chosen <- which(path$lambda == attr(fit, "settings")$lambda)
  expect_true(path$final[chosen])
  expect_identical(chosen, which.min(path$ebic))
  ft <- as.data.frame(fit, what = "fit")
  ebic_returned <- -2 * ft$logl + ft$npar * log(100) +
    4 * ft$npar * 0.5 * log(12)
  expect_equal(ebic_returned, path$ebic[chosen], tolerance = 1e-8)
})

test_that("final EBIC path certificates describe the rounded coefficients", {
  skip_on_cran()
  w <- panel_sim_wide(n = 200, seed = 1)
  cols <- grep("^V1_", names(w$data))
  w$data[cols] <- w$data[cols] * 1e5
  # Capture the optimizer solutions, then independently certify their rounded
  # estimates by the forward-mode score, rather than the solver's gradient.
  solve <- psychnets:::.pgv_pen_fit
  records <- list()
  testthat::local_mocked_bindings(.pgv_pen_fit = function(ctx, theta0, pen, w,
                                                         lambda, ...) {
    fit <- solve(ctx, theta0, pen, w, lambda, ...)
    if (!is.null(fit)) {
      theta <- fit$theta
      theta[pen & abs(theta * w) < 1e-5] <- 0
      st <- psychnets:::.pgv_state(theta, ctx)
      jac <- psychnets:::.pgv_jacobian(st, ctx, seq_along(theta))
      g <- psychnets:::.pgv_score(st, ctx, jac)
      residual <- abs(g)
      nonzero <- pen & theta != 0
      zero <- pen & theta == 0
      residual[nonzero] <- abs(g[nonzero] +
        lambda * w[nonzero] * sign(theta[nonzero]))
      residual[zero] <- pmax(0, abs(g[zero]) - lambda * w[zero])
      records[[length(records) + 1L]] <<- c(lambda = lambda,
        objective = fit$objective, score = max(residual))
    }
    fit
  }, .package = "psychnets")
  fit <- panel_gvar(w$data, vars = w$design, select = "ebic", nlambda = 8)
  path <- attr(fit, "settings")$path
  records <- as.data.frame(do.call(rbind, records))
  for (i in which(path$final)) {
    # A final solve uses three starts after the scan; choose the same optimum
    # and compare its independently evaluated, post-threshold certificate.
    trials <- tail(records[records$lambda == path$lambda[i], ], 3L)
    expected <- trials$score[which.min(trials$objective)]
    expect_equal(path$kkt[i], expected, tolerance = 1e-7)
  }
  expect_true(any(path$kkt[path$final] > 0.1))
  chosen <- which.min(path$ebic)
  expect_equal(path$kkt[chosen], attr(fit, "kkt"), tolerance = 1e-8)
})

test_that("only the optimizer's bound is an active constraint", {
  w <- panel_sim_wide()
  fit <- panel_gvar(w$data, vars = w$design)
  ctx <- pgv_ctx(w)                          # unit-scale data: scales ~1
  pp <- as.data.frame(fit, what = "parameters")
  j <- which(pp$matrix == "delta_zeta_between" & pp$row == "V2")
  theta <- pp$estimate
  theta[j] <- 5e-5                           # near zero, but interior
  st <- psychnets:::.pgv_state(theta, ctx)
  expect_false(psychnets:::.pgv_at_bound(st, ctx)[j])
  g <- psychnets:::.pgv_gradient_full(st, ctx)
  inf <- suppressWarnings(psychnets:::.pgv_inference(st, ctx))
  expect_gte(inf$kkt, abs(g[j]))             # its score is not suppressed
  theta[j] <- 1e-8                           # on the bound
  expect_true(psychnets:::.pgv_at_bound(psychnets:::.pgv_state(theta, ctx),
                                        ctx)[j])
})

test_that("identified SEs keep identifiable combinations of unidentified pairs", {
  # Only b2 - b3 is unidentified; b1 is conditioned on b2 + b3, not on nothing.
  h <- matrix(c(2, 1, 1, 1, 1, 1, 1, 1, 1), 3)
  expect_equal(psychnets:::.pgv_identified_se(h), c(1, Inf, Inf),
               tolerance = 1e-10)
  expect_equal(psychnets:::.pgv_identified_se(h)[1],
               sqrt(solve(matrix(c(2, 1, 1, 1), 2))[1, 1]), tolerance = 1e-10)
  # Rank is judged on the correlation scale: one enormous, perfectly identified
  # coordinate must not make the rest look null (a raw relative cut-off of
  # 1e-10 x 1e12 would flag the whole second block).
  big <- matrix(0, 4, 4)
  big[1, 1] <- 1e12
  big[2:4, 2:4] <- h
  expect_equal(psychnets:::.pgv_identified_se(big), c(1e-6, 1, Inf, Inf),
               tolerance = 1e-10)
})

test_that("units of any size fit identically (scale equivariance)", {
  long <- panel_sim(n = 400, seed = 1)
  v <- c("V1", "V2", "V3")
  a <- fit_long(long)
  tiny <- long
  tiny[v] <- long[v] * 1e-6
  huge <- long
  huge$V2 <- long$V2 * 1e5
  b <- fit_long(tiny)
  d <- fit_long(huge)
  for (net in c("contemporaneous", "between")) {
    expect_equal(b[[net]]$weights, a[[net]]$weights, tolerance = 1e-6)
    expect_equal(d[[net]]$weights, a[[net]]$weights, tolerance = 1e-6)
  }
  scale <- c(1, 1e5, 1)
  expect_equal(d$temporal$beta, a$temporal$beta * outer(scale, 1 / scale),
               tolerance = 1e-6)
  expect_equal(as.data.frame(d, what = "fit")$chisq,
               as.data.frame(a, what = "fit")$chisq, tolerance = 1e-6)
  expect_lt(attr(d, "kkt"), 1e-8)
})

test_that("a singular information keeps exact SEs for identified parameters", {
  # Three waves, saturated: the between-person block is not identified. Its
  # SEs are NA; every other parameter keeps a finite standard error.
  long <- panel_sim(n = 500, p = 3, waves = 3, seed = 7)
  expect_warning(
    fit <- withCallingHandlers(fit_long(long),
      psychnet_boundary = function(w) invokeRestart("muffleWarning")),
    class = "psychnet_singular_information")
  pp <- as.data.frame(fit, what = "parameters")
  between <- pp$matrix %in% c("omega_zeta_between", "delta_zeta_between")
  expect_true(all(is.finite(pp$se[!between])))
  expect_true(any(is.na(pp$se[between])))
})

# ---- standardize = "z" / "z_per_wave" (psychonetrics panelvar(standardize =)) ----

test_that("standardize equals fitting the hand-standardized data", {
  w <- panel_sim_wide()
  w$data$V2_w3 <- w$data$V2_w3 + 2          # drift in a mean ...
  w$data$V1_w4 <- w$data$V1_w4 * 1.5        # ... and in a variance
  cols <- as.vector(w$design)
  per_wave <- w$data
  per_wave[cols] <- lapply(w$data[cols], function(x) (x - mean(x)) / stats::sd(x))
  pooled <- w$data
  for (v in rownames(w$design)) {
    vals <- unlist(w$data[w$design[v, ]])
    pooled[w$design[v, ]] <- lapply(w$data[w$design[v, ]],
                                    function(x) (x - mean(vals)) / stats::sd(vals))
  }
  est <- function(fit) as.data.frame(fit, what = "parameters")$estimate
  expect_equal(est(panel_gvar(w$data, vars = w$design, standardize = "z_per_wave")),
               est(panel_gvar(per_wave, vars = w$design)), tolerance = 1e-10)
  expect_equal(est(panel_gvar(w$data, vars = w$design, standardize = "z")),
               est(panel_gvar(pooled, vars = w$design)), tolerance = 1e-10)
  # Per-wave z-scoring removes the drift: every wave mean is exactly zero, so
  # the fitted stationary means are zero too.
  fz <- panel_gvar(w$data, vars = w$design, standardize = "z_per_wave")
  mu <- subset(as.data.frame(fz, what = "parameters"), matrix == "mu")$estimate
  expect_lt(max(abs(mu)), 1e-8)
  expect_identical(attr(fz, "settings")$standardize, "z_per_wave")
})

test_that("standardize refuses a variable constant within a wave", {
  long <- panel_sim(n = 200)
  long$V1[long$wave == 2] <- 3
  expect_error(fit_long(long, standardize = "z_per_wave"),
               class = "psychnet_bad_input")
})

test_that("srl_panel is a complete, balanced long panel on a 0-100 scale", {
  vars <- c("efficacy", "planning", "effort", "enjoyment", "anxiety")
  expect_named(srl_panel, c("student", "week", vars))
  expect_identical(nrow(srl_panel), 400L)
  expect_identical(length(unique(srl_panel$student)), 100L)
  expect_true(all(table(srl_panel$student) == 4L))
  expect_setequal(unique(srl_panel$week), 1:4)
  expect_false(anyNA(srl_panel))
  expect_true(all(unlist(srl_panel[vars]) >= 0 & unlist(srl_panel[vars]) <= 100))
  fit <- panel_gvar(srl_panel, vars = vars, id = "student", wave = "week")
  expect_lt(attr(fit, "kkt"), 1e-8)
  expect_false(as.data.frame(fit, what = "diagnostics")$boundary)
})
