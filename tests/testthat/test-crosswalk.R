# net_crosswalk: tidy argument map (psychnet as a substitute for the references).

test_that("net_crosswalk returns a tidy one-row-per-argument data.frame", {
  cw <- net_crosswalk()
  expect_s3_class(cw, "data.frame")
  expect_named(cw, c("reference", "psychnet", "ref_arg", "psychnet_arg",
                     "status", "note"))
  expect_gt(nrow(cw), 50L)
  expect_true(all(cw$status %in% c("identical", "renamed", "default differs",
                                   "semantics differ", "reference only",
                                   "psychnet only")))
  # every row maps a real argument on at least one side (no "-"/"-" rows)
  expect_false(any(cw$ref_arg == "-" & cw$psychnet_arg == "-"))
})

test_that("reference= filters to a single estimator and they cover all six", {
  for (r in c("EBICglasso", "cor_auto", "ggmModSelect", "IsingFit", "mgm",
              "panelgvar")) {
    tab <- net_crosswalk(r)
    expect_gt(nrow(tab), 0L)
    expect_equal(unique(sub("^.*::", "", tab$reference)), r)
  }
  expect_error(net_crosswalk("graphicalVAR"))   # time-series VAR: out of scope
})

test_that("crosswalk psychnet_arg names match the real formals", {
  # every non-'-' psychnet_arg must be an actual argument of its verb
  cw <- net_crosswalk()
  verbs <- list(ebic_glasso = ebic_glasso, cor_auto = cor_auto,
                ggm_modselect = ggm_modselect, ising_fit = ising_fit,
                mgm_fit = mgm_fit, panel_gvar = panel_gvar)
  for (v in names(verbs)) {
    rows <- cw[cw$psychnet == v & cw$psychnet_arg != "-", ]
    expect_true(all(rows$psychnet_arg %in% names(formals(verbs[[v]]))),
                info = v)
  }
})
