# Regenerate the panel_gvar() fixtures.
#
# Two files:
#   golden-panelgvar.rds        psychnets' own frozen output (numerical-stability
#                               policy; checked at 1e-10 by test-golden.R).
#   psychonetrics-panelgvar.rds the reference: psychonetrics::panelgvar() fits
#                               of the same data (saturated + pruned), so the
#                               equivalence test runs without psychonetrics
#                               installed. Regenerate when psychonetrics changes.
#
# Regenerate golden-panelgvar.rds ONLY when a change of numbers is intentional,
# and record it in NEWS.md.
#
# Run from the package root with the package loaded:
#   Rscript -e 'devtools::load_all("."); source("tests/testthat/fixtures/make-panelgvar.R")'

source(file.path("tests", "testthat", "helper-panelgvar.R"))
dir <- file.path("tests", "testthat", "fixtures")

# The reference datasets: the default design, a 5 x 5 design, a small sample
# where BH adjustment prunes differently from no adjustment, and two designs on
# which stepup() re-adds a parameter pruning removed.
cases <- list(
  default = list(sim = list(n = 400, p = 3, waves = 4, seed = 1),
                 alpha = 0.01, adjust = "none"),
  large   = list(sim = list(n = 600, p = 5, waves = 5, seed = 8),
                 alpha = 0.01, adjust = "none"),
  bh      = list(sim = list(n = 150, p = 4, waves = 4, seed = 10),
                 alpha = 0.05, adjust = "BH"),
  stepup_small = list(sim = list(n = 45, p = 4, waves = 4, seed = 21),
                      alpha = 0.01, adjust = "none"),
  stepup_waves = list(sim = list(n = 300, p = 4, waves = 5, seed = 12),
                      alpha = 0.01, adjust = "none")
)

own <- lapply(cases, function(cs) {
  w <- do.call(panel_sim_wide, cs$sim)
  sat <- panel_gvar(w$data, vars = w$design, select = "none")
  pru <- panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                    adjust = cs$adjust, stepup = FALSE)
  up <- panel_gvar(w$data, vars = w$design, select = "prune", alpha = cs$alpha,
                   adjust = cs$adjust)
  list(saturated = as.data.frame(sat, what = "parameters")$estimate,
       pruned = as.data.frame(pru, what = "parameters")$estimate,
       stepup = as.data.frame(up, what = "parameters")$estimate)
})
# EBIC lasso (select = "ebic"): frozen estimates and the selected lambda for
# three designs (psychonetrics PML references below).
ebic_cases <- c("default", "bh", "stepup_waves")
own$ebic <- lapply(cases[ebic_cases], function(cs) {
  w <- do.call(panel_sim_wide, cs$sim)
  fit <- panel_gvar(w$data, vars = w$design, select = "ebic")
  list(estimate = as.data.frame(fit, what = "parameters")$estimate,
       lambda = attr(fit, "settings")$lambda)
})
saveRDS(own, file.path(dir, "golden-panelgvar.rds"))

if (requireNamespace("psychonetrics", quietly = TRUE)) {
  ref <- lapply(cases, function(cs) {
    w <- do.call(panel_sim_wide, cs$sim)
    sat <- psychonetrics::runmodel(
      psychonetrics::panelgvar(w$data, vars = w$design, verbose = FALSE),
      verbose = FALSE)
    pru <- psychonetrics::prune(sat, alpha = cs$alpha, adjust = cs$adjust,
                                verbose = FALSE)
    up <- psychonetrics::stepup(pru, alpha = cs$alpha, verbose = FALSE)
    fit_tab <- psychonetrics::fit(sat)
    list(estimate = sat@parameters$est, se = sat@parameters$se,
         p = sat@parameters$p, objective = sat@objective,
         fit = stats::setNames(fit_tab$Value, fit_tab$Measure),
         pruned_estimate = pru@parameters$est,
         pruned_fixed = pru@parameters$fixed,
         pruned_objective = pru@objective,
         pruned_mi = pru@parameters$mi,
         stepup_estimate = up@parameters$est,
         stepup_fixed = up@parameters$fixed)
  })
  ref$pml <- lapply(cases[ebic_cases], function(cs) {
    w <- do.call(panel_sim_wide, cs$sim)
    m <- psychonetrics::runmodel(
      psychonetrics::panelgvar(w$data, vars = w$design, estimator = "PML",
                               verbose = FALSE), verbose = FALSE)
    list(estimate = m@parameters$est,
         lambda = m@optim$lambda_search$best_lambda,
         lambda_max = m@optim$lambda_search$lambda_max)
  })
  attr(ref, "psychonetrics_version") <-
    as.character(utils::packageVersion("psychonetrics"))
  saveRDS(ref, file.path(dir, "psychonetrics-panelgvar.rds"))
}
