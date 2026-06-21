# =====================================================================
# Case study: reproducing a published psychological-network analysis.
#
# Source: Learning Analytics Methods (LAMethods), Book 1, Chapter 19,
#   "Psychological Networks" (lamethods.org/book1).
# Data:   university COVID well-being survey, Finland and Austria, from
#   https://github.com/sonsoleslp/labook-data (11_universityCovid/data.sav).
#
# The chapter estimates a 6-construct Gaussian graphical model with
#   bootnet::estimateNetwork(default = "EBICglasso", corMethod = "cor_auto",
#                            tuning = 0.5)
# and reports: 15 edges, mean edge weight ~0.14, mean node predictability
# R^2 ~0.333, Competence and Motivation highest in expected influence,
# and a centrality-stability coefficient ~0.95.
#
# This script reproduces that pipeline and compares psychnet against
# bootnet/qgraph on the same cor_auto correlation. Unlike run_validation.R,
# it downloads one external file (the .sav), so it needs an internet
# connection on first run; the file is cached in tempdir().
#
# Reproduce:  Rscript validation/case_study_lamethods.R   (from package root)
# =====================================================================

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(".", quiet = TRUE) else library(psychnet)
  for (p in c("haven", "qgraph", "bootnet", "mgm")) suppressWarnings(library(p, character.only = TRUE))
})

url  <- "https://github.com/sonsoleslp/labook-data/raw/main/11_universityCovid/data.sav"
dest <- file.path(tempdir(), "labook_universityCovid.sav")
if (!file.exists(dest)) utils::download.file(url, dest, mode = "wb", quiet = TRUE)
d <- haven::read_sav(dest)
d <- d[stats::complete.cases(d), ]                      # chapter: import(...) |> drop_na()

# Six constructs = row mean of each scale's three recoded items.
mk <- function(pre) rowMeans(d[, paste0(pre, 1:3, ".rec")])
ad <- data.frame(Competence = mk("comp"), Autonomy = mk("auto"),
                 Motivation = mk("lm"),   Emotion  = mk("pa"),
                 Relatedness = mk("gp"),  SRL = mk("sr"))
n  <- nrow(ad)
cat(sprintf("n after listwise deletion = %d (chapter: 6071)\n\n", n))

S <- qgraph::cor_auto(ad, verbose = FALSE)              # chapter's corMethod
ut <- upper.tri(S)

# --- network: psychnet vs bootnet/qgraph, both on the same cor_auto ---
boot <- bootnet::estimateNetwork(ad, default = "EBICglasso",
                                 corMethod = "cor_auto", tuning = 0.5)
B <- boot$graph
fit <- estimate_network(NULL, "EBICglasso", cor_matrix = S, n = n, gamma = 0.5)
P <- fit$graph

cat("Network (EBICglasso):\n")
cat(sprintf("  bootnet/qgraph : %d edges, mean|w| %.4f, max|w| %.4f\n",
            sum(B[ut] != 0), mean(abs(B[ut])), max(abs(B[ut]))))
cat(sprintf("  psychnet       : %d edges, mean|w| %.4f, max|w| %.4f\n",
            sum(P[ut] != 0), mean(abs(P[ut])), max(abs(P[ut]))))
cat(sprintf("  agreement      : max edge delta %.6f, structure %.3f\n",
            max(abs(P[ut] - B[ut])), mean((abs(P[ut]) > 1e-6) == (abs(B[ut]) > 1e-6))))

# --- centrality: expected influence + strength ---
cp <- psychnet::centrality(P); cq <- qgraph::centrality(P)
cat(sprintf("\nCentrality (psychnet vs qgraph): expected influence max|diff| %.6f, strength max|diff| %.6f\n",
            max(abs(cp$expected_influence - cq$OutExpectedInfluence)),
            max(abs(cp$strength - cq$OutDegree))))
top <- cp$node[order(-cp$expected_influence)][1:2]
cat(sprintf("  highest expected influence: %s (chapter: Competence, Motivation)\n",
            paste(top, collapse = ", ")))

# --- predictability: closed-form GGM R^2 ---
pr <- psychnet::predictability(fit)
cat(sprintf("\nPredictability: mean R^2 %.3f (chapter mgm: 0.333)\n", mean(pr$predictability)))

# --- centrality stability (CS coefficient), case-dropping resampling ---
cs <- tryCatch(psychnet::centrality_stability(ad, method = "EBICglasso", iter = 200L)$cs,
               error = function(e) c(strength = NA, expected_influence = NA))
cat(sprintf("\nCentrality stability (CS coefficient): strength %.2f, expected influence %.2f (chapter: ~0.95)\n",
            cs[["strength"]], cs[["expected_influence"]]))

# --- write a one-row result record ---
res <- data.frame(
  case = "LAMethods Ch.19 (university COVID)", source = "lamethods.org/book1",
  p = ncol(P), n = n, edges = sum(P[ut] != 0),
  mean_weight = round(mean(abs(P[ut])), 4),
  max_edge_delta_vs_bootnet = round(max(abs(P[ut] - B[ut])), 6),
  struct_vs_bootnet = round(mean((abs(P[ut]) > 1e-6) == (abs(B[ut]) > 1e-6)), 4),
  centrality_max_diff = round(max(abs(cp$expected_influence - cq$OutExpectedInfluence)), 6),
  mean_R2 = round(mean(pr$predictability), 3),
  cs_strength = round(cs[["strength"]], 2), cs_expinf = round(cs[["expected_influence"]], 2),
  stringsAsFactors = FALSE)
utils::write.csv(res, "validation/results_lamethods.csv", row.names = FALSE)
cat("\nWrote validation/results_lamethods.csv\n")
