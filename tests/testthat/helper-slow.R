# Resampling verbs (net_boot, net_compare, net_stability, casedrop_reliability,
# network_reliability, difference_test) re-estimate the whole network once per
# draw, so a single test_that() block can cost several seconds. Run in full
# locally and in CI; skipped on CRAN to keep R CMD check inside the 10-minute
# budget CRAN enforces.
#
# Run the complete suite the way CI does with:
#   NOT_CRAN=true Rscript -e 'testthat::test_dir("tests/testthat")'
skip_slow <- function() {
  testthat::skip_on_cran()
}
