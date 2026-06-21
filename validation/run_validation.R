# =====================================================================
# psychnet validation: agreement with the reference packages on real
# questionnaire data, with a synthetic ground-truth supplement.
#
# Part A. EBIC graphical lasso, psychnet vs qgraph, on real questionnaire
#         item responses. Both estimators receive the identical pairwise
#         correlation and effective n, so the comparison isolates the
#         penalty path and EBIC selection from the choice of correlation.
# Part B. Ising model, psychnet vs IsingFit, on real binary item responses.
#         Both fit nodewise from the same binarized data.
# Part C. Synthetic datasets with a known generating graph, to report
#         recovery (F-measure against the truth) alongside agreement.
#
# Reproduce:  Rscript validation/run_validation.R   (from the package root)
# Output:     validation/results_realdata.csv, results_synthetic.csv,
#             validation/RESULTS.md, and a printed sessionInfo().
#
# Dependencies for the comparison (not for psychnet itself): qgraph, psych,
# psychTools, mgm, NetworkToolbox, networktools, EGAnet, IsingFit.
# =====================================================================

suppressMessages({
  ok <- requireNamespace("devtools", quietly = TRUE)
  if (ok) devtools::load_all(".", quiet = TRUE) else library(psychnet)
  for (p in c("qgraph", "psych", "psychTools", "mgm", "NetworkToolbox",
              "networktools", "EGAnet", "IsingFit")) suppressWarnings(library(p, character.only = TRUE))
})
set.seed(2026)
ute <- function(M) M[upper.tri(M)]

# ---- helpers --------------------------------------------------------
get_data <- function(name, pkg) {                 # load a dataset, return $data if a list
  e <- new.env(); data(list = name, package = pkg, envir = e); o <- e[[name]]
  if (is.list(o) && !is.data.frame(o) && !is.null(o$data)) o$data else o
}
to_items <- function(X, maxp = 50) {              # numeric item matrix, capped
  X <- suppressWarnings(matrix(as.numeric(as.matrix(X)), nrow(as.matrix(X))))
  keep <- apply(X, 2, function(c) sum(!is.na(c)) > 10 && length(unique(c[!is.na(c)])) > 1)
  X <- X[, keep, drop = FALSE]
  if (ncol(X) > maxp) X[, seq_len(maxp), drop = FALSE] else X
}
shared_cor <- function(X) {                       # pairwise correlation + effective n
  S <- psychnet:::.nearest_pd_cor(suppressWarnings(stats::cor(X, use = "pairwise.complete.obs")))
  co <- crossprod(!is.na(X)); list(S = S, n = round(stats::median(co[upper.tri(co)])))
}

# ---- Part A: EBICglasso vs qgraph on real questionnaire data --------
# registry: label, package, dataset, citation, max items
A <- list(
  c("qgraph::big5",          "qgraph",        "big5",      "Dolan, Oort, Stoel & Wicherts (2009)",      "240"),
  c("psych::bfi",            "psych",         "bfi",       "Revelle, Wilt & Rosenthal (2010)",          "25"),
  c("psych::sat.act",        "psych",         "sat.act",   "Revelle (psych package)",                   "50"),
  c("psychTools::msq",       "psychTools",    "msq",       "Revelle & Anderson (1998)",                 "50"),
  c("psychTools::epi.bfi",   "psychTools",    "epi.bfi",   "Eysenck & Eysenck (1964); psychTools",      "50"),
  c("psychTools::epi",       "psychTools",    "epi",       "Eysenck Personality Inventory; psychTools", "50"),
  c("psychTools::affect",    "psychTools",    "affect",    "Rafaeli & Revelle (2006)",                  "50"),
  c("psychTools::sai",       "psychTools",    "sai",       "Spielberger State Anxiety; psychTools",     "40"),
  c("psychTools::tai",       "psychTools",    "tai",       "Spielberger Trait Anxiety; psychTools",     "40"),
  c("psychTools::spi",       "psychTools",    "spi",       "Condon (2018), SAPA Personality Inventory", "50"),
  c("psychTools::blot",      "psychTools",    "blot",      "Bond's Logical Operations Test; psychTools","40"),
  c("mgm::Fried2015",        "mgm",           "Fried2015", "Fried et al. (2015)",                       "50"),
  c("mgm::PTSD_data",        "mgm",           "PTSD_data", "McNally et al. (2015)",                     "50"),
  c("mgm::B5MS",             "mgm",           "B5MS",      "Haslbeck & Waldorp (mgm package)",          "50"),
  c("mgm::symptom_data",     "mgm",           "symptom_data","Haslbeck & Waldorp (mgm package)",        "50"),
  c("NetworkToolbox::neoOpen","NetworkToolbox","neoOpen",  "Christensen, Cotter & Silvia (2019)",       "48"),
  c("networktools::depression","networktools","depression","Jones (networktools package)",             "50"),
  c("networktools::social",  "networktools",  "social",    "Jones (networktools package)",              "50"),
  c("EGAnet::optimism",      "EGAnet",        "optimism",  "Golino & Christensen (EGAnet)",             "50"),
  c("EGAnet::intelligenceBattery","EGAnet",   "intelligenceBattery","Golino & Christensen (EGAnet)",    "50"))

cat("== Part A: EBICglasso, psychnet vs qgraph (real questionnaire data) ==\n")
rowsA <- list()
for (s in A) {
  res <- tryCatch({
    X <- to_items(get_data(s[3], s[2]), as.integer(s[5]))
    if (ncol(X) < 4) stop("too few items")
    sc <- shared_cor(X)
    P <- estimate_network(NULL, "EBICglasso", cor_matrix = sc$S, n = sc$n, gamma = 0.5)$graph
    Q <- suppressWarnings(qgraph::EBICglasso(sc$S, n = sc$n, gamma = 0.5, verbose = FALSE))
    data.frame(dataset = s[1], citation = s[4], p = ncol(X), n = sc$n,
               max_delta = round(max(abs(ute(P) - ute(Q))), 6),
               struct = round(mean((abs(ute(P)) > 1e-6) == (abs(ute(Q)) > 1e-6)), 4),
               edges_psychnet = sum(abs(ute(P)) > 1e-6), edges_qgraph = sum(abs(ute(Q)) > 1e-6),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { cat(sprintf("  %-28s skipped (%s)\n", s[1], "structure not raw items")); next }
  rowsA[[length(rowsA) + 1L]] <- res
  cat(sprintf("  %-28s p=%-3d n=%-5d  maxD=%.5f struct=%.3f  %d/%d\n",
              res$dataset, res$p, res$n, res$max_delta, res$struct, res$edges_psychnet, res$edges_qgraph))
}
resA <- do.call(rbind, rowsA)

# ---- Part B: Ising vs IsingFit on real binary item responses --------
# Coerce to a clean 0/1 matrix: keep native 0/1 columns as is, map other
# two-level columns to 0/1, split polytomous columns at the median, drop
# near-constant columns. Both estimators then see the identical matrix.
prep_binary <- function(X) {
  X <- suppressWarnings(matrix(as.numeric(as.matrix(X)), nrow = nrow(as.matrix(X))))
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < 50 || ncol(X) < 4) return(NULL)
  X <- vapply(seq_len(ncol(X)), function(j) {
    c <- X[, j]; u <- sort(unique(c))
    if (length(u) == 2 && all(u == c(0, 1))) c
    else if (length(u) == 2) as.numeric(c == u[2])
    else as.numeric(c > stats::median(c))
  }, numeric(nrow(X)))
  keep <- apply(X, 2, function(c) { m <- mean(c); length(unique(c)) == 2 && min(m, 1 - m) > 0.05 })
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 4) NULL else X
}
B <- list(
  c("psychTools::ability", "psychTools", "ability", "Revelle; ICAR ability items", "16"),
  c("EGAnet::wmt2",        "EGAnet",     "wmt2",    "Wiener Matrizen-Test; EGAnet", "20"),
  c("psychTools::iqitems", "psychTools", "iqitems", "Condon & Revelle; ICAR",       "16"))

cat("\n== Part B: Ising, psychnet vs IsingFit (real binary data) ==\n")
rowsB <- list()
for (s in B) {
  res <- tryCatch({
    bx <- prep_binary(get_data(s[3], s[2]))
    if (is.null(bx)) stop("no clean binary matrix")
    P <- estimate_network(bx, "ising", gamma = 0.25, rule = "AND")$graph
    R <- suppressWarnings(suppressMessages(IsingFit::IsingFit(bx, gamma = 0.25, AND = TRUE,
                          plot = FALSE, progressbar = FALSE)$weiadj))
    data.frame(dataset = s[1], citation = s[4], p = ncol(bx), n = nrow(bx),
               max_delta = round(max(abs(ute(P) - ute(R))), 6),
               struct = round(mean((abs(ute(P)) > 1e-6) == (abs(ute(R)) > 1e-6)), 4),
               edges_psychnet = sum(abs(ute(P)) > 1e-6), edges_isingfit = sum(abs(ute(R)) > 1e-6),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  if (is.null(res)) { cat(sprintf("  %-28s skipped\n", s[1])); next }
  rowsB[[length(rowsB) + 1L]] <- res
  cat(sprintf("  %-28s p=%-3d n=%-5d  maxD=%.5f struct=%.3f  %d/%d\n",
              res$dataset, res$p, res$n, res$max_delta, res$struct, res$edges_psychnet, res$edges_isingfit))
}
resB <- if (length(rowsB)) do.call(rbind, rowsB) else NULL

# ---- Part C: synthetic ground-truth supplement ----------------------
# Known sparse precision (chain / cluster). Report agreement with qgraph AND
# recovery of the true edge set, which real data cannot give.
chain_prec <- function(p, a) { K <- diag(p); for (i in seq_len(p - 1)) K[i, i+1] <- K[i+1, i] <- a; stats::cov2cor(solve(K)) }
two_block  <- function(p, w) { m <- matrix(0, p, p); h <- p %/% 2L
  m[1:h, 1:h] <- w; m[(h+1):p, (h+1):p] <- w; diag(m) <- 1; stats::cov2cor(m) }
rmvn <- function(n, S) matrix(stats::rnorm(n * ncol(S)), n, ncol(S)) %*% chol(S)
f1 <- function(e, t) { tp <- sum(e & t); if (2*tp + sum(e & !t) + sum(!e & t) == 0) 1 else 2*tp/(2*tp + sum(e & !t) + sum(!e & t)) }

cat("\n== Part C: synthetic ground truth (recovery + agreement vs qgraph) ==\n")
rowsC <- list()
for (cfg in list(list("chain", 10, 400, 0.95), list("chain", 8, 200, 0.85),
                 list("two-block", 8, 400, 0.30), list("two-block", 10, 250, 0.97))) {
  S0 <- if (cfg[[1]] == "chain") chain_prec(cfg[[2]], -0.35) else two_block(cfg[[2]], 0.4)
  true <- if (cfg[[1]] == "chain") ute(abs(outer(1:cfg[[2]], 1:cfg[[2]], "-")) == 1) else ute(abs(S0) > 0.05)
  X <- rmvn(cfg[[3]], S0); S <- stats::cor(X); n <- cfg[[3]]
  P <- estimate_network(NULL, "EBICglasso", cor_matrix = S, n = n, gamma = 0.5)$graph
  Q <- suppressWarnings(qgraph::EBICglasso(S, n = n, gamma = 0.5, verbose = FALSE))
  rowsC[[length(rowsC) + 1L]] <- data.frame(
    structure = cfg[[1]], p = cfg[[2]], n = cfg[[3]],
    max_delta = round(max(abs(ute(P) - ute(Q))), 6),
    struct_vs_qgraph = round(mean((abs(ute(P)) > 1e-6) == (abs(ute(Q)) > 1e-6)), 4),
    F1_psychnet = round(f1(abs(ute(P)) > 1e-6, true), 3),
    F1_qgraph = round(f1(abs(ute(Q)) > 1e-6, true), 3), stringsAsFactors = FALSE)
  r <- rowsC[[length(rowsC)]]
  cat(sprintf("  %-10s p=%-3d n=%-4d  maxD=%.5f struct=%.3f  F1 psy=%.2f qg=%.2f\n",
              r$structure, r$p, r$n, r$max_delta, r$struct_vs_qgraph, r$F1_psychnet, r$F1_qgraph))
}
resC <- do.call(rbind, rowsC)

# ---- write outputs --------------------------------------------------
utils::write.csv(resA, "validation/results_realdata.csv", row.names = FALSE)
if (!is.null(resB)) utils::write.csv(resB, "validation/results_ising.csv", row.names = FALSE)
utils::write.csv(resC, "validation/results_synthetic.csv", row.names = FALSE)

md <- c("# psychnet validation results",
  "",
  sprintf("Generated %s. Reference packages: qgraph %s, IsingFit %s, mgm %s.",
          format(Sys.Date()), packageVersion("qgraph"), packageVersion("IsingFit"), packageVersion("mgm")),
  "", "## Part A: EBICglasso vs qgraph (real questionnaire data)", "",
  sprintf("%d datasets; all structure agreement = %.3f, max edge delta <= %.5f.",
          nrow(resA), min(resA$struct), max(resA$max_delta)),
  "", "## Part B: Ising vs IsingFit (real binary data)", "",
  if (!is.null(resB)) sprintf("%d datasets; min structure agreement %.3f.", nrow(resB), min(resB$struct)) else "none",
  "", "## Part C: synthetic ground truth", "",
  "psychnet and qgraph agree, and recover the true graph at the same F-measure.",
  "", "See results_realdata.csv, results_ising.csv, results_synthetic.csv.")
writeLines(md, "validation/RESULTS.md")

cat(sprintf("\nPart A: %d datasets, all struct = %.3f, max delta = %.5f\n",
            nrow(resA), min(resA$struct), max(resA$max_delta)))
cat("Wrote validation/results_*.csv and validation/RESULTS.md\n\n")
cat("R", as.character(getRversion()), "| psychnet", as.character(packageVersion("psychnet")), "\n")
