# Writes the input files the equivalence tests share (tests/equivalence/data). These are inputs, not
# reference values: every test computes its reference from them in the run.
#
#   bfi10.csv          ten items of psych::bfi (A1-A5, C1-C5), complete rows; six ordered answers
#   bfi10_gender.csv   the same items with gender (1 = male, 2 = female), complete rows
#   bfi6_binary.csv    A1-A5 and C1 of those rows as binary answers: 1 = agree (4, 5 or 6)
#   network.tsv        a signed weighted network of the ten items: qgraph::EBICglasso() of their
#                      Pearson correlations (gamma = 0.5), rounded to six decimals
#
# Run from the package root on a machine with psych and qgraph (plain Rscript; not a test):
#   Rscript tests/equivalence/generators/make-data.R
stopifnot("psych is needed" = requireNamespace("psych", quietly = TRUE),
          "qgraph is needed" = requireNamespace("qgraph", quietly = TRUE))
out <- file.path("tests", "equivalence", "data")
items <- c(paste0("A", seq_len(5L)), paste0("C", seq_len(5L)))
bfi10 <- stats::na.omit(psych::bfi[items])
utils::write.csv(bfi10, file.path(out, "bfi10.csv"), row.names = FALSE)
with_gender <- stats::na.omit(psych::bfi[c(items, "gender")])
utils::write.csv(with_gender, file.path(out, "bfi10_gender.csv"), row.names = FALSE)
binary <- as.data.frame(lapply(bfi10[c(paste0("A", seq_len(5L)), "C1")], \(item) as.integer(item >= 4L)))
utils::write.csv(binary, file.path(out, "bfi6_binary.csv"), row.names = FALSE)
# qgraph warns that a dense network was selected; that is about interpreting small edges, not this file.
network <- withCallingHandlers(
  qgraph::EBICglasso(stats::cor(bfi10), n = nrow(bfi10), gamma = 0.5, verbose = FALSE),
  warning = \(w) if (grepl("dense regularized network", conditionMessage(w))) invokeRestart("muffleWarning"))
utils::write.table(round(network, 6L), file.path(out, "network.tsv"), sep = "\t", quote = FALSE, col.names = NA)
message(sprintf("psych %s, qgraph %s: wrote %s", utils::packageVersion("psych"), utils::packageVersion("qgraph"),
                toString(list.files(out))))
