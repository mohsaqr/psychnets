# Regenerate the network comparison test fixtures used by test-nct.R.
#
# nct-baseline.rds  frozen net_compare() output on complete data (regression:
#                   input-handling changes must not move a complete-data number).
# The comparison with NetworkComparisonTest::NCT() runs live in tests/equivalence/test-net-compare.R.
#
# Regenerate ONLY when a change of numbers is intentional (record it in NEWS.md).
# Run from the package root with the package loaded:
#   Rscript -e 'devtools::load_all("."); source("tests/testthat/fixtures/make-nct.R")'

source(file.path("tests", "testthat", "helper-nct.R"))
dir <- file.path("tests", "testthat", "fixtures")

run_cfg <- function(i, fn) {
  cfg <- nct_configs[i, ]
  g <- nct_groups(cfg$seed)
  set.seed(cfg$seed)
  fn(g, cfg)
}

baseline <- lapply(seq_len(nrow(nct_configs)), run_cfg, fn = function(g, cfg)
  nct_numbers(net_compare(g$a, g$b, iter = nct_iter, paired = cfg$paired,
                          abs = cfg$abs, weighted = cfg$weighted)))
saveRDS(baseline, file.path(dir, "nct-baseline.rds"))

message("wrote nct fixtures to ", dir)
