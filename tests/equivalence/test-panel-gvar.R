# panel_gvar() and its accessors: the panel graphical VAR (Epskamp 2020) by maximum likelihood, against
# psychonetrics::panelgvar() |> runmodel(), run live on simulated panels (a stationary random-intercept
# VAR(1)). Saturated estimates, standard errors and fit indices; the supports chosen by pruning and by
# pruning followed by step-up search.

designs <- list(default = list(n = 400, p = 3, waves = 4, seed = 1, alpha = 0.01, adjust = "none"),
                bh = list(n = 150, p = 4, waves = 4, seed = 10, alpha = 0.05, adjust = "BH"),
                large = list(n = 600, p = 5, waves = 5, seed = 8, alpha = 0.01, adjust = "none"))

#' psychnets' parameter table with the mean's empty column named as in ref_psychonetrics_parameters()
own_parameters <- function(fit) {
  table <- as.data.frame(fit, what = "parameters")
  transform(table, col = ifelse(is.na(col), "(mean)", col))
}

saturated <- function(panel) {
  psychonetrics::runmodel(psychonetrics::panelgvar(panel$data, vars = panel$design, verbose = FALSE), verbose = FALSE)
}

invisible(lapply(names(designs), \(name) {
  design <- designs[[name]]
  test_that(sprintf("saturated estimates, standard errors and fit equal psychonetrics, design %s", name), {
    skip_without_reference("psychonetrics")
    panel <- ref_panel(n = design$n, p = design$p, waves = design$waves, seed = design$seed)
    reference <- saturated(panel)
    ours <- panel_gvar(panel$data, vars = panel$design, select = "none")
    expected <- ref_psychonetrics_parameters(reference@parameters, rownames(panel$design))
    # psychonetrics stops nlminb with a gradient of about 1e-3 (measured 6.8e-4 on the default design), which
    # moves estimates and standard errors by up to about 5e-3; psychnets' own fit is certified to 1e-8.
    expect_reference(own_parameters(ours), expected[c("matrix", "row", "col", "estimate", "se")],
                     reference = "psychonetrics", covers = c("panel_gvar()", "as.data.frame.psychnet_panelgvar()"),
                     quantity = sprintf("estimate and standard error, design %s", name), tolerance = 5e-3,
                     key = c("matrix", "row", "col"))
    indices <- psychonetrics::fit(reference)
    value <- \(measure) indices$Value[indices$Measure == measure]
    # chi-square, log likelihood and BIC to 1e-3 at that optimizer stop; CFI, TLI and RMSEA to 1e-5 (the
    # package's own test).
    expect_reference(as.data.frame(ours, what = "fit"),
                     data.frame(df = value("df"), logl = value("logl"), chisq = value("chisq"),
                                baseline_chisq = value("baseline.chisq"), bic = value("bic")),
                     reference = "psychonetrics", covers = "as.data.frame.psychnet_panelgvar()",
                     quantity = sprintf("fit statistics, design %s", name), tolerance = 1e-3, key = "df")
    expect_reference(as.data.frame(ours, what = "fit"),
                     data.frame(df = value("df"), cfi = value("cfi"), tli = value("tli"), rmsea = value("rmsea")),
                     reference = "psychonetrics", covers = "as.data.frame.psychnet_panelgvar()",
                     quantity = sprintf("fit indices, design %s", name), tolerance = 1e-5, key = "df")
  })

  test_that(sprintf("the pruned and step-up supports equal psychonetrics::prune() and stepup(), design %s", name), {
    skip_without_reference("psychonetrics")
    panel <- ref_panel(n = design$n, p = design$p, waves = design$waves, seed = design$seed)
    vars <- rownames(panel$design)
    pruned_reference <- psychonetrics::prune(saturated(panel), alpha = design$alpha, adjust = design$adjust,
                                             verbose = FALSE)
    stepped_reference <- psychonetrics::stepup(pruned_reference, alpha = design$alpha, verbose = FALSE)
    pruned <- panel_gvar(panel$data, vars = panel$design, select = "prune", alpha = design$alpha,
                         adjust = design$adjust, stepup = FALSE)
    stepped <- panel_gvar(panel$data, vars = panel$design, select = "prune", alpha = design$alpha,
                          adjust = design$adjust)
    support <- \(table) table[c("matrix", "row", "col", "fixed")]
    # TRUE / FALSE: equal or not.
    expect_reference(support(own_parameters(pruned)),
                     support(ref_psychonetrics_parameters(pruned_reference@parameters, vars)),
                     reference = "psychonetrics", covers = "panel_gvar()",
                     quantity = sprintf("parameter fixed to zero by pruning, design %s", name), tolerance = 0,
                     key = c("matrix", "row", "col"))
    expect_reference(support(own_parameters(stepped)),
                     support(ref_psychonetrics_parameters(stepped_reference@parameters, vars)),
                     reference = "psychonetrics", covers = "panel_gvar()",
                     quantity = sprintf("parameter fixed to zero after step-up, design %s", name), tolerance = 0,
                     key = c("matrix", "row", "col"))
  })
}))

test_that("the networks of the saturated model equal psychonetrics' estimates as edges", {
  skip_without_reference("psychonetrics")
  panel <- ref_panel()
  reference <- ref_psychonetrics_parameters(saturated(panel)@parameters, rownames(panel$design))
  temporal <- subset(reference, matrix == "beta")
  # psychnets' temporal edge `from -> to` is beta[to, from]: the effect of from at t - 1 on to at t.
  expected <- data.frame(network = "temporal", from = temporal$col, to = temporal$row, weight = temporal$estimate)
  ours <- subset(as.data.frame(panel_gvar(panel$data, vars = panel$design, select = "none"), include_zero = TRUE),
                 network == "temporal")
  expect_reference(ours, expected, reference = "psychonetrics", covers = "as.data.frame.psychnet_panelgvar()",
                   quantity = "temporal edge (beta)", tolerance = 5e-3, key = c("network", "from", "to"))
  ours_summary <- subset(summary(panel_gvar(panel$data, vars = panel$design, select = "none")), network == "temporal")
  # The summary's mean absolute weight leaves out the autoregressive (self-loop) effects.
  cross_lagged <- temporal$estimate[temporal$row != temporal$col]
  expect_reference(ours_summary, data.frame(network = "temporal", mean_abs_weight = mean(abs(cross_lagged))),
                   reference = "psychonetrics", covers = "summary.psychnet_panelgvar()",
                   quantity = "mean absolute cross-lagged weight", tolerance = 5e-3, key = "network")
})
