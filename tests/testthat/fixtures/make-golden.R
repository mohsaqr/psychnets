# Regenerate the golden regression fixtures used by test-golden.R.
#
# These freeze the *current* numeric output of the estimators that downstream
# consumers pin against (Nestimate holds frozen baselines at 1e-12, and it has a
# CRAN reverse dependency), so a refactor that moves a weight is caught here
# before it is caught there. Regenerate ONLY when a change of numbers is
# intentional, and record it in NEWS.md per the numerical-stability policy.
#
# Run from the package root with the package loaded:
#   Rscript -e 'devtools::load_all("."); source("tests/testthat/fixtures/make-golden.R")'
#
# The data generators live in tests/testthat/helper-golden.R so this script and
# the tests that check its output can never drift apart.

source(file.path("tests", "testthat", "helper-golden.R"))

dir <- file.path("tests", "testthat", "fixtures")
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

saveRDS(
  list(
    # The glmnet paths are the ones Nestimate delegates to for bit-exact
    # mgm::mgm / IsingFit parity; the base paths are psychnets' own default.
    mgm_glmnet    = mgm_fit(golden_data_mixed(), native = FALSE)$weights,
    mgm_base      = mgm_fit(golden_data_mixed(), native = TRUE)$weights,
    ising_glmnet  = ising_fit(golden_data_binary(), native = FALSE)$weights,
    ising_base    = ising_fit(golden_data_binary(), native = TRUE)$weights,
    glasso_base   = ebic_glasso(cor_matrix = golden_cor(), n = 300)$weights,
    glasso_lambda = ebic_glasso(cor_matrix = golden_cor(), n = 300)$lambda
  ),
  file.path(dir, "golden-estimators.rds")
)

message("wrote ", file.path(dir, "golden-estimators.rds"))
