# Regenerate the network comparison test fixtures used by test-nct.R.
#
# nct-baseline.rds  frozen net_compare() output on complete data (regression:
#                   input-handling changes must not move a complete-data number).
# nct-reference.rds NetworkComparisonTest::NCT() on the same data and seed,
#                   so the ungated suite checks parity without a Suggests entry.
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

if (requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  reference <- lapply(seq_len(nrow(nct_configs)), run_cfg, fn = function(g, cfg) {
    ref <- suppressWarnings(suppressMessages(NetworkComparisonTest::NCT(
      g$a, g$b, gamma = 0.5, it = nct_iter, paired = cfg$paired, abs = cfg$abs,
      weighted = cfg$weighted, test.edges = TRUE, edges = "all",
      progressbar = FALSE)))
    # NCT field names: glstrinv = global strength (S), nwinv = structure (M)
    list(strength_perm = ref$glstrinv.perm, strength_p = ref$glstrinv.pval,
         strength_by_group = ref$glstrinv.sep,
         structure_perm = ref$nwinv.perm, structure_p = ref$nwinv.pval,
         E_p = ref$einv.pvals[["p-value"]], nw1 = unname(ref$nw1),
         nw2 = unname(ref$nw2),
         version = as.character(utils::packageVersion("NetworkComparisonTest")))
  })
  saveRDS(reference, file.path(dir, "nct-reference.rds"))
}
message("wrote nct fixtures to ", dir)
