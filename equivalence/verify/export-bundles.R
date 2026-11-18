# Exports psychnets comparison cases as bundles for Prova Verify (format:
# plan/24-verify.md in the prova repository), with Prova's R writer.
#
# Every live side of a case is an R script, stored in the bundle and run there
# in its own R process: its console is the side's raw output and the table it
# leaves in `result` holds the side's values. Two cases set psychnets beside a
# reference that is not installed here (NetworkComparisonTest, psychonetrics);
# their reference side is RETAINED: the numbers the package's own test fixtures
# stored when the reference was run (tests/testthat/fixtures/make-nct.R and
# make-panelgvar.R).
#
# Run from the package root:
#   PROVA_VERIFY_R=<prova>/reporters/verify.R \
#     Rscript equivalence/verify/export-bundles.R [case ...]
source(Sys.getenv("PROVA_VERIFY_R", file.path("..", "prova", "reporters", "verify.R")))

.verify_base <- file.path("equivalence", "verify")
.verify_generator <- "psychnets equivalence/verify/export-bundles.R"
.verify_fixtures <- file.path("tests", "testthat", "fixtures")
.psychnets_version <- unname(read.dcf("DESCRIPTION")[, "Version"])

# ---- data --------------------------------------------------------------------
# All data files are written to one temporary folder and copied into the bundles
# that name them. Numbers are written by write.csv() (15 significant digits).

.verify_data_dir <- function() {
  folder <- file.path(tempdir(), "psychnets-verify-data")
  dir.create(folder, showWarnings = FALSE)
  folder
}

.write_csv <- function(table, name) {
  file <- file.path(.verify_data_dir(), name)
  utils::write.csv(table, file, row.names = FALSE)
  file
}

#' The bundled SRL_Claude data (300 rows, five scale scores)
.data_srl <- function() {
  scope <- new.env()
  load(file.path("data", "SRL_Claude.rda"), envir = scope)
  .write_csv(scope$SRL_Claude, "srl.csv")
}

#' Ten items of psych::bfi (A1-A5, C1-C5), complete rows; six answer categories
.data_bfi10 <- function() {
  items <- c(paste0("A", seq_len(5L)), paste0("C", seq_len(5L)))
  .write_csv(stats::na.omit(psych::bfi[items]), "bfi10.csv")
}

#' Six of those items (A1-A5, C1) as binary answers: 1 = agree (4, 5 or 6)
.data_bfi6_binary <- function() {
  items <- c(paste0("A", seq_len(5L)), "C1")
  answers <- stats::na.omit(psych::bfi[items])
  .write_csv(as.data.frame(lapply(answers, function(item) as.integer(item >= 4L))), "bfi6_binary.csv")
}

#' A signed weighted network of ten nodes: the EBIC graphical lasso network of
#' the ten bfi items as estimated by qgraph, rounded to six decimals
.data_network <- function() {
  items <- utils::read.csv(.data_bfi10())
  # qgraph warns that a dense network was selected; that warning is about the
  # interpretation of small edges and says nothing about this file.
  network <- withCallingHandlers(
    qgraph::EBICglasso(stats::cor(items), n = nrow(items), gamma = 0.5, verbose = FALSE),
    warning = function(w) if (grepl("dense regularized network", conditionMessage(w))) invokeRestart("muffleWarning"))
  file <- file.path(.verify_data_dir(), "network.tsv")
  utils::write.table(round(network, 6L), file, sep = "\t", quote = FALSE, col.names = NA)
  file
}

#' The two groups of the package's NCT fixture, configuration 1 (seed 1)
.data_nct <- function() {
  scope <- new.env()
  sys.source(file.path("tests", "testthat", "helper-nct.R"), envir = scope)
  groups <- scope$nct_groups(1)
  c("Group 1" = .write_csv(as.data.frame(groups$a), "group_a.csv"),
    "Group 2" = .write_csv(as.data.frame(groups$b), "group_b.csv"))
}

#' The panel of the package's panel GVAR fixture, "default" design, long form
.data_panel <- function() {
  scope <- new.env()
  sys.source(file.path("tests", "testthat", "helper-panelgvar.R"), envir = scope)
  .write_csv(scope$panel_sim(n = 400, p = 3, waves = 4, seed = 1), "panel.csv")
}

# ---- sides -------------------------------------------------------------------

.psychnets <- function(script, key = NULL) {
  verify_side("psychnets", package = "psychnets", version = .psychnets_version, primary = TRUE,
              load = ".", key = key, script = c("library(psychnets)", script))
}
.reference <- function(package, script, key = NULL, id = tolower(package), attach = TRUE) {
  verify_side(id, name = package, package = package, authority = "open-source", key = key,
              script = c(if (attach) sprintf("library(%s)", package), script))
}
.base_r <- function(script, key = NULL) {
  verify_side("stats", name = "base R (stats)", package = "stats", authority = "independent", key = key,
              script = script)
}

# The reference scripts name an edge as psychnets' as.data.frame() does: the
# upper triangle, column by column, "from -- to".
.edge_names <- c(
  "edge_names <- function(m) {",
  "  upper <- upper.tri(m)",
  "  paste(rownames(m)[row(m)[upper]], \"--\", colnames(m)[col(m)[upper]])",
  "}")
.read_network <- "network <- as.matrix(utils::read.delim(\"../data/network.tsv\", row.names = 1))"

# ---- cases whose sides are all live ------------------------------------------
# Each tolerance was fixed before the case was first exported; the comment says
# why that value. `sides` is a function, so that a side of a reference package
# that is not installed is never described.

.verify_cases <- function() {
  list(
    correlation_networks = list(
      title = "Correlation and partial correlation networks",
      description = paste(
        "Pearson correlations and partial correlations of the five scales of the bundled SRL_Claude data",
        "(n = 300). Base R computes the partial correlations from the inverse of the correlation matrix;",
        "qgraph computes both networks with qgraph(graph = \"cor\" / \"pcor\"). No threshold on any side."),
      covers = c("cor_network()", "pcor_network()"),
      requires = "qgraph", data = function() c("SRL_Claude" = .data_srl()),
      # Closed-form quantities (a correlation, one matrix inverse of a 5 x 5 matrix): only rounding
      # of the last bits may differ.
      tolerance = 1e-10,
      sides = function() list(
        .psychnets(c(
          "srl <- utils::read.csv(\"../data/srl.csv\")",
          "correlations <- as.data.frame(cor_network(data = srl, cor_method = \"pearson\"))",
          "partials <- as.data.frame(pcor_network(data = srl, cor_method = \"pearson\"))",
          "correlations",
          "partials",
          "result <- rbind(",
          "  with(correlations, data.frame(quantity = paste(\"cor:\", from, \"--\", to), kind = \"correlation\", value = weight)),",
          "  with(partials, data.frame(quantity = paste(\"pcor:\", from, \"--\", to), kind = \"partial correlation\", value = weight)))")),
        .base_r(c(
          "srl <- utils::read.csv(\"../data/srl.csv\")",
          .edge_names,
          "correlations <- cor(srl, method = \"pearson\")",
          "partials <- -cov2cor(solve(correlations))",
          "result <- rbind(",
          "  data.frame(quantity = paste(\"cor:\", edge_names(correlations)), kind = \"correlation\",",
          "             value = correlations[upper.tri(correlations)]),",
          "  data.frame(quantity = paste(\"pcor:\", edge_names(partials)), kind = \"partial correlation\",",
          "             value = partials[upper.tri(partials)]))",
          "result")),
        .reference("qgraph", c(
          "srl <- utils::read.csv(\"../data/srl.csv\")",
          .edge_names,
          "correlations <- getWmat(qgraph(cor(srl), graph = \"cor\", DoNotPlot = TRUE))",
          "partials <- getWmat(qgraph(cor(srl), graph = \"pcor\", DoNotPlot = TRUE))",
          "result <- rbind(",
          "  data.frame(quantity = paste(\"cor:\", edge_names(correlations)), kind = \"correlation\",",
          "             value = correlations[upper.tri(correlations)]),",
          "  data.frame(quantity = paste(\"pcor:\", edge_names(partials)), kind = \"partial correlation\",",
          "             value = partials[upper.tri(partials)]))",
          "result")))),

    ebic_glasso = list(
      title = "EBIC graphical lasso network",
      description = paste(
        "The EBIC graphical lasso network of ten items of psych::bfi (A1-A5, C1-C5; 2,632 complete rows),",
        "from Pearson correlations. Both sides: gamma = 0.5, 100 penalties spaced on the log scale from",
        "max|r| down to 0.01 x max|r|, diagonal not penalized, no thresholding, the regularized network",
        "returned (no unregularized refit). psychnets uses its own solver (native = TRUE) and returns the",
        "optimum at the selected penalty; qgraph uses the glasso package, here with its convergence",
        "threshold tightened from the default 1e-4 to 1e-10 (thr, passed through EBICglasso's dots) so that",
        "the two solvers can be read against each other. All 45 node pairs are compared, the zeros included,",
        "and the selected penalty."),
      covers = "ebic_glasso()",
      requires = "qgraph", data = function() c("bfi items" = .data_bfi10()),
      # Two independent solvers of one strictly convex problem. psychnets certifies its optimum to a
      # KKT residual below 1e-7 and glasso is asked for thr = 1e-10, so partial correlations should
      # agree far below 1e-6. (A probe at glasso's DEFAULT thr = 1e-4 differed by 1.4e-5; that run is
      # not this case.) The same value covers the selected penalty, which both sides compute in closed
      # form from the same grid.
      tolerance = 1e-6,
      sides = function() list(
        .psychnets(c(
          "items <- utils::read.csv(\"../data/bfi10.csv\")",
          "fit <- ebic_glasso(data = items, gamma = 0.5, nlambda = 100L, lambda_min_ratio = 0.01,",
          "                   threshold = 0, cor_method = \"pearson\", penalize_diagonal = FALSE,",
          "                   native = TRUE, refit = TRUE)",
          "fit",
          "edges <- as.data.frame(fit, include_zero = TRUE)",
          "result <- rbind(",
          "  with(edges, data.frame(quantity = paste(\"weight:\", from, \"--\", to), kind = \"partial correlation\", value = weight)),",
          "  data.frame(quantity = \"selected lambda\", kind = \"penalty\", value = fit$lambda))",
          "result")),
        .reference("qgraph", c(
          "items <- utils::read.csv(\"../data/bfi10.csv\")",
          .edge_names,
          "fit <- EBICglasso(cor(items), n = nrow(items), gamma = 0.5, nlambda = 100, lambda.min.ratio = 0.01,",
          "                  penalize.diagonal = FALSE, refit = FALSE, threshold = FALSE,",
          "                  returnAllResults = TRUE, thr = 1e-10, maxit = 1e5)",
          "network <- fit$optnet",
          "result <- rbind(",
          "  data.frame(quantity = paste(\"weight:\", edge_names(network)), kind = \"partial correlation\",",
          "             value = network[upper.tri(network)]),",
          "  data.frame(quantity = \"selected lambda\", kind = \"penalty\", value = fit$lambda[which.min(fit$ebic)]))",
          "result")))),

    polychoric_correlations = list(
      title = "Polychoric correlations of ordinal items",
      description = paste(
        "Polychoric correlations of the ten psych::bfi items (six ordered categories each), as used by",
        "cor_method = \"auto\". psych::polychoric() and qgraph::cor_auto() (which calls lavaan) are two",
        "independent two-step estimators; psych runs with its defaults, which include a continuity",
        "correction of 0.5 for empty cells."),
      covers = c("cor_auto()", "cor_network()"),
      requires = c("psych", "qgraph", "lavaan"), data = function() c("bfi items" = .data_bfi10()),
      # The package's own test (tests/testthat/test-cor_auto.R) holds cor_auto() to 1e-3 of both
      # references: three two-step estimators with different one-dimensional optimizers and different
      # treatment of empty cells. That existing tolerance is used unchanged.
      tolerance = 1e-3,
      sides = function() list(
        .psychnets(c(
          "items <- utils::read.csv(\"../data/bfi10.csv\")",
          "cor_auto(data = items)",
          "result <- with(as.data.frame(cor_network(data = items, cor_method = \"auto\")),",
          "               data.frame(quantity = paste(from, \"--\", to), kind = \"polychoric correlation\", value = weight))",
          "result")),
        .reference("psych", c(
          "items <- utils::read.csv(\"../data/bfi10.csv\")",
          .edge_names,
          "correlations <- polychoric(items)$rho",
          "result <- data.frame(quantity = edge_names(correlations), kind = \"polychoric correlation\",",
          "                     value = correlations[upper.tri(correlations)])",
          "result")),
        .reference("qgraph", c(
          "items <- utils::read.csv(\"../data/bfi10.csv\")",
          .edge_names,
          "correlations <- cor_auto(items, verbose = TRUE)",
          "result <- data.frame(quantity = edge_names(correlations), kind = \"polychoric correlation\",",
          "                     value = correlations[upper.tri(correlations)])",
          "result")))),

    centralities = list(
      title = "Node centralities of a signed weighted network",
      description = paste(
        "Strength, expected influence, closeness and betweenness of a ten-node signed network (the EBIC",
        "graphical lasso network of the bfi items, rounded to six decimals). Shortest paths run over",
        "1 / |weight|. psychnets reports closeness and betweenness normalized as igraph does",
        "(normalized = TRUE); qgraph::centrality_auto() reports closeness as 1 / sum of distances and",
        "betweenness as a count of paths, so the qgraph script multiplies closeness by (p - 1) and divides",
        "betweenness by (p - 1)(p - 2) / 2, as psychnets' documentation states the relation."),
      covers = "net_centralities()",
      requires = c("qgraph", "igraph"), data = function() c(Network = .data_network()),
      # Deterministic sums and shortest paths on the same matrix; the package's tests use 1e-8, and
      # nothing iterative is involved, so only rounding may differ.
      tolerance = 1e-10,
      sides = function() list(
        .psychnets(key = "node", c(
          .read_network,
          "result <- net_centralities(network, measures = c(\"strength\", \"expected_influence\", \"closeness\", \"betweenness\"))",
          "result")),
        .reference("qgraph", key = "node", c(
          .read_network,
          "central <- centrality_auto(network)$node.centrality",
          "p <- nrow(network)",
          "result <- data.frame(node = rownames(central), strength = central$Strength,",
          "                     expected_influence = central$ExpectedInfluence,",
          "                     closeness = central$Closeness * (p - 1),",
          "                     betweenness = central$Betweenness / ((p - 1) * (p - 2) / 2))",
          "result")),
        .reference("igraph", key = "node", c(
          .read_network,
          "graph <- graph_from_adjacency_matrix(abs(network), mode = \"undirected\", weighted = TRUE)",
          "lengths <- 1 / E(graph)$weight",
          "# expected influence is the signed row sum; igraph has no function of that name",
          "result <- data.frame(node = V(graph)$name, strength = strength(graph),",
          "                     expected_influence = rowSums(network),",
          "                     closeness = closeness(graph, weights = lengths, normalized = TRUE),",
          "                     betweenness = betweenness(graph, weights = lengths, normalized = TRUE))",
          "result")))),

    clustering = list(
      title = "Local clustering coefficients of a signed weighted network",
      description = paste(
        "The seven local clustering coefficients (Watts-Strogatz, Zhang, Onnela, each unsigned and",
        "signed, and Barrat) of the same ten-node signed network, against qgraph::clustcoef_auto()."),
      covers = "net_clustering()",
      requires = "qgraph", data = function() c(Network = .data_network()),
      # Closed-form ratios of sums of products of weights on the same matrix.
      tolerance = 1e-10,
      sides = function() list(
        .psychnets(key = "node", c(
          .read_network,
          "result <- net_clustering(network)",
          "result")),
        .reference("qgraph", key = "node", c(
          .read_network,
          "coefficients <- clustcoef_auto(network)",
          "result <- data.frame(node = rownames(coefficients), coefficients)",
          "result")))),

    ising_unregularized = list(
      title = "Unregularized Ising model against nodewise logistic regression",
      description = paste(
        "An Ising model of six binary items (psych::bfi A1-A5 and C1, 1 = agree; 2,632 rows) without a",
        "penalty and without pruning (alpha = NULL). The reference is the definition of the estimator",
        "written with stats::glm(): each item regressed on the other five by logistic regression, an",
        "edge being the mean of its two coefficients (AND rule on an unpruned model) and a threshold the",
        "intercept."),
      covers = "ising_sampler()",
      requires = "psych", data = function() c("Binary bfi items" = .data_bfi6_binary()),
      # Two IRLS solvers of six concave problems. glm() is run to epsilon = 1e-12 and psychnets
      # certifies its score to ~1e-8, so coefficients of order 1 should agree well below 1e-6.
      tolerance = 1e-6,
      sides = function() list(
        .psychnets(c(
          "answers <- utils::read.csv(\"../data/bfi6_binary.csv\")",
          "fit <- ising_sampler(data = answers, rule = \"AND\", alpha = NULL)",
          "fit",
          "result <- rbind(",
          "  with(as.data.frame(fit), data.frame(quantity = paste(\"weight:\", from, \"--\", to), kind = \"edge weight\", value = weight)),",
          "  data.frame(quantity = paste(\"threshold:\", names(fit$thresholds)), kind = \"threshold\", value = unname(fit$thresholds)))",
          "result")),
        .base_r(c(
          "answers <- as.matrix(utils::read.csv(\"../data/bfi6_binary.csv\"))",
          .edge_names,
          "items <- colnames(answers)",
          "# one column per regressed item: its intercept, then its coefficient on every item (0 on itself)",
          "nodewise <- sapply(items, function(item) {",
          "  others <- setdiff(items, item)",
          "  fit <- glm(answers[, item] ~ answers[, others], family = binomial,",
          "             control = glm.control(epsilon = 1e-12, maxit = 100))",
          "  replace(setNames(numeric(length(items) + 1), c(\"threshold\", items)), c(\"threshold\", others), coef(fit))",
          "})",
          "nodewise",
          "weights <- (nodewise[items, ] + t(nodewise[items, ])) / 2",
          "result <- rbind(",
          "  data.frame(quantity = paste(\"weight:\", edge_names(weights)), kind = \"edge weight\", value = weights[upper.tri(weights)]),",
          "  data.frame(quantity = paste(\"threshold:\", items), kind = \"threshold\", value = nodewise[\"threshold\", ]))",
          "result")))),

    # NOT RUN in the session that wrote this file (IsingFit was not installed): the case is skipped
    # with a classed message until it is. Settings follow tests/testthat/test-equivalence.R.
    ising_regularized = list(
      title = "EBIC-regularized Ising model",
      description = paste(
        "The eLasso Ising model of six binary bfi items against IsingFit::IsingFit(): gamma = 0.25, AND",
        "rule, 100 penalties per node. psychnets uses its own solver (native = TRUE)."),
      covers = "ising_fit()",
      requires = c("psych", "IsingFit"), data = function() c("Binary bfi items" = .data_bfi6_binary()),
      # psychnets' own solver against glmnet's path (IsingFit): independent coordinate-descent solvers
      # at glmnet's default convergence threshold (1e-7 of the null deviance).
      tolerance = 1e-4,
      sides = function() list(
        .psychnets(c(
          "answers <- utils::read.csv(\"../data/bfi6_binary.csv\")",
          "fit <- ising_fit(data = answers, gamma = 0.25, rule = \"AND\", nlambda = 100L, native = TRUE)",
          "fit",
          "result <- with(as.data.frame(fit, include_zero = TRUE),",
          "               data.frame(quantity = paste(\"weight:\", from, \"--\", to), kind = \"edge weight\", value = weight))",
          "result")),
        .reference("IsingFit", c(
          "answers <- utils::read.csv(\"../data/bfi6_binary.csv\")",
          .edge_names,
          "fit <- IsingFit(answers, gamma = 0.25, AND = TRUE, plot = FALSE, progressbar = FALSE)",
          "result <- data.frame(quantity = paste(\"weight:\", edge_names(fit$weiadj)), kind = \"edge weight\",",
          "                     value = fit$weiadj[upper.tri(fit$weiadj)])",
          "result"))))
  )
}

#' Export one case whose sides are all live
.export_live <- function(name, case) {
  verify_export(bundle = file.path(.verify_base, "bundles", name), title = case$title,
                description = case$description, sides = case$sides(), data = case$data(),
                tolerance = case$tolerance, covers = case$covers, generator = .verify_generator)
}

# ---- cases with a retained reference -----------------------------------------

#' Export a case: psychnets live, the reference kept from a test fixture
#'
#' @param name The case's name.
#' @param title,description,covers The case's words.
#' @param side The psychnets side (a script that leaves a long `result`).
#' @param reference The retained side (`verify_side(source = "retained")`).
#' @param data Named data files, as `verify_export()`.
#' @param retained A function of the psychnets `result` returning the reference's
#'   values: a data frame with `quantity`, `kind`, `value`.
#' @param maker The fixture script that ran the reference; stored as its code.
#' @param tolerances A named numeric vector: the tolerance of each `kind`.
#' @return One row, as `verify_export()`.
.export_retained <- function(name, title, description, covers, side, reference, data, retained, maker,
                             tolerances) {
  bundle <- file.path(.verify_base, "bundles", name)
  unlink(bundle, recursive = TRUE)
  dir.create(bundle, recursive = TRUE)
  shared <- unname(Map(function(file, label) verify_copy(file, bundle, file.path("data", basename(file)), label = label),
                       unname(data), names(data)))
  run <- verify_run(bundle, side, reads = vapply(shared, function(item) item$path, character(1)))
  own <- verify_values(run$result)
  kept <- retained(own)
  stopifnot("every kind needs a tolerance" = all(c(own$kind, kept$kind) %in% names(tolerances)))
  kept_file <- file.path(bundle, reference$id, "values.csv")
  dir.create(dirname(kept_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(transform(kept, value = format(value, digits = 17L)), kept_file, row.names = FALSE, quote = FALSE)
  comparisons <- do.call(c, lapply(names(tolerances), function(kind) {
    pair <- stats::setNames(list(own[own$kind == kind, ], kept[kept$kind == kind, ]), c(side$id, reference$id))
    if (nrow(pair[[1L]]) + nrow(pair[[2L]]) == 0L) return(list())
    verify_compare(pair, tolerance = tolerances[[kind]], baseline = reference$id)
  }))
  sections <- list(
    list(id = "code", title = "Code", items = list(
      run$code,
      verify_copy(maker, bundle, file.path(reference$id, basename(maker)), side = reference$id,
                  label = "Fixture script that ran the reference", type = "code", language = "r"))),
    list(id = "output", title = "Raw output", items = list(
      run$output,
      verify_item(bundle, file.path(reference$id, "values.csv"), side = reference$id,
                  label = "Values kept in the fixture", type = "table"))),
    list(id = "data", title = "Data", items = shared))
  verify_write(bundle, title = title, description = description, sides = list(side, reference),
               sections = sections, comparisons = comparisons, covers = covers, generator = .verify_generator)
  states <- vapply(comparisons, function(row) row$status, character(1))
  data.frame(bundle = bundle, sides = 2L, compared = length(states), agree = sum(states == "agree"),
             disagree = sum(states == "disagree"), missing = sum(states == "missing"))
}

.export_network_comparison <- function(name) {
  fixture <- readRDS(file.path(.verify_fixtures, "nct-reference.rds"))[[1L]]
  nodes <- paste0("V", seq_len(5L))
  upper <- upper.tri(fixture$nw1)
  edges <- paste(nodes[row(fixture$nw1)[upper]], "--", nodes[col(fixture$nw1)[upper]])
  .export_retained(
    name,
    title = "Network comparison test (permutation)",
    description = paste(
      "The network comparison test of two independent groups (n = 120 each, five variables; the",
      "package's fixture configuration 1): EBIC glasso networks at gamma = 0.5, absolute weighted",
      "differences, 15 permutations after set.seed(1). The reference is NetworkComparisonTest::NCT() on",
      "the same data and seed, kept in tests/testthat/fixtures/nct-reference.rds (the reference package",
      "is not installed here, so its side is retained, not rerun). With the same seed both draw the same",
      "permutations, so the p-values are comparable one to one. The fixture holds edge p-values in NCT's",
      "edge order (upper triangle, column by column), which is the order psychnets reports."),
    covers = "net_compare()",
    side = .psychnets(c(
      "group_a <- utils::read.csv(\"../data/group_a.csv\")",
      "group_b <- utils::read.csv(\"../data/group_b.csv\")",
      "set.seed(1)",
      "fit <- net_compare(data1 = group_a, data2 = group_b, iter = 15L, gamma = 0.5, paired = FALSE,",
      "                   abs = TRUE, weighted = TRUE, p_adjust = \"none\")",
      "fit",
      "tests <- summary(fit)",
      "edges <- as.data.frame(fit)",
      "result <- rbind(",
      "  with(tests, data.frame(quantity = paste(\"p:\", test), kind = \"p value\", value = p_value)),",
      "  with(edges, data.frame(quantity = paste(\"p: edge\", from, \"--\", to), kind = \"p value\", value = p_value)),",
      "  with(edges, data.frame(quantity = paste(\"weight, group 1:\", from, \"--\", to), kind = \"edge weight\", value = weight_1)),",
      "  with(edges, data.frame(quantity = paste(\"weight, group 2:\", from, \"--\", to), kind = \"edge weight\", value = weight_2)),",
      "  data.frame(quantity = paste(\"global strength, group\", seq_len(2L)), kind = \"global strength\", value = fit$strength$by_group))",
      "result")),
    reference = verify_side("nct", name = "NetworkComparisonTest", version = fixture$version,
                            authority = "open-source", source = "retained",
                            role = "kept in tests/testthat/fixtures/nct-reference.rds",
                            command = "NetworkComparisonTest::NCT(a, b, gamma = 0.5, it = 15, paired = FALSE, abs = TRUE, weighted = TRUE, test.edges = TRUE, edges = \"all\")"),
    data = .data_nct(),
    retained = function(own) rbind(
      data.frame(quantity = c("p: network structure", "p: global strength"), kind = "p value",
                 value = c(fixture$structure_p, fixture$strength_p)),
      data.frame(quantity = paste("p: edge", edges), kind = "p value", value = fixture$E_p),
      data.frame(quantity = paste("weight, group 1:", edges), kind = "edge weight", value = fixture$nw1[upper]),
      data.frame(quantity = paste("weight, group 2:", edges), kind = "edge weight", value = fixture$nw2[upper]),
      data.frame(quantity = paste("global strength, group", seq_len(2L)), kind = "global strength",
                 value = fixture$strength_by_group)),
    maker = file.path(.verify_fixtures, "make-nct.R"),
    # From the package's own test (tests/testthat/test-nct.R): NCT estimates with qgraph, whose glasso
    # stops at thr = 1e-4, so weights and strengths agree to 1e-3. A permutation p-value moves in steps
    # of 1 / 16, so the same 1e-3 demands the identical p-value.
    tolerances = c("p value" = 1e-3, "edge weight" = 1e-3, "global strength" = 1e-3))
}

.export_panel_gvar <- function(name) {
  fixtures <- readRDS(file.path(.verify_fixtures, "psychonetrics-panelgvar.rds"))
  fixture <- fixtures$default
  indices <- c(logl = "logl", chisq = "chisq", baseline_chisq = "baseline.chisq", bic = "bic",
               cfi = "cfi", tli = "tli", rmsea = "rmsea")
  .export_retained(
    name,
    title = "Panel graphical VAR, saturated model",
    description = paste(
      "The saturated panel graphical VAR (temporal, contemporaneous and between-person networks) of a",
      "simulated panel: 400 persons, four waves, three variables (the package's fixture design",
      "\"default\"), by maximum likelihood without pruning. The reference is",
      "psychonetrics::panelgvar() |> runmodel() on the same data, kept in",
      "tests/testthat/fixtures/psychonetrics-panelgvar.rds (psychonetrics is not installed here, so its",
      "side is retained, not rerun). The fixture holds the 24 estimates and standard errors as bare",
      "vectors in psychonetrics' parameter order; they are named here by position from psychnets'",
      "parameter table, the pairing the package's own test uses. psychnets is called with the long",
      "panel; the fixture was made from the same panel in wide form."),
    covers = "panel_gvar()",
    side = .psychnets(c(
      "panel <- utils::read.csv(\"../data/panel.csv\")",
      "fit <- panel_gvar(data = panel, vars = c(\"V1\", \"V2\", \"V3\"), id = \"id\", wave = \"wave\", select = \"none\")",
      "fit",
      "parameters <- as.data.frame(fit, what = \"parameters\")",
      "parameters",
      "indices <- utils::stack(as.data.frame(fit, what = \"fit\"),",
      "                        select = c(logl, chisq, baseline_chisq, bic, cfi, tli, rmsea))",
      "result <- rbind(",
      "  with(parameters, data.frame(quantity = paste(\"estimate:\", matrix, row, col), kind = \"estimate\", value = estimate)),",
      "  with(parameters, data.frame(quantity = paste(\"se:\", matrix, row, col), kind = \"standard error\", value = se)),",
      "  with(indices, data.frame(quantity = paste(\"fit:\", ind), kind = ifelse(ind %in% c(\"cfi\", \"tli\", \"rmsea\"), \"fit index\", \"fit statistic\"), value = values)))",
      "result")),
    reference = verify_side("psychonetrics", name = "psychonetrics",
                            version = attr(fixtures, "psychonetrics_version"),
                            authority = "open-source", source = "retained",
                            role = "kept in tests/testthat/fixtures/psychonetrics-panelgvar.rds",
                            command = "psychonetrics::runmodel(psychonetrics::panelgvar(data, vars = design))"),
    data = c("Panel (long)" = .data_panel()),
    retained = function(own) {
      estimates <- subset(own, kind == "estimate")
      errors <- subset(own, kind == "standard error")
      stopifnot("the fixture and the parameter table must have the same length" =
                  nrow(estimates) == length(fixture$estimate) && nrow(errors) == length(fixture$se))
      rbind(
        data.frame(quantity = estimates$quantity, kind = "estimate", value = fixture$estimate),
        data.frame(quantity = errors$quantity, kind = "standard error", value = fixture$se),
        data.frame(quantity = paste("fit:", names(indices)),
                   kind = ifelse(names(indices) %in% c("cfi", "tli", "rmsea"), "fit index", "fit statistic"),
                   value = unname(fixture$fit[indices])))
    },
    maker = file.path(.verify_fixtures, "make-panelgvar.R"),
    # From the package's own test (tests/testthat/test-panelgvar.R): psychonetrics stops nlminb with a
    # gradient of about 1e-3, which moves estimates and standard errors by up to 5e-3; chi-square, the
    # log likelihood and BIC are held to 1e-3 and CFI, TLI and RMSEA to 1e-5. (The test holds standard
    # errors to a RATIO within 5e-3; the absolute 5e-3 here is looser for errors below 1.)
    tolerances = c("estimate" = 5e-3, "standard error" = 5e-3, "fit statistic" = 1e-3, "fit index" = 1e-5))
}

# ---- run ---------------------------------------------------------------------

.verify_retained <- list(
  network_comparison = list(export = .export_network_comparison, requires = character()),
  panel_gvar = list(export = .export_panel_gvar, requires = character()))

#' Export one case, or say why it was skipped
#'
#' @param name The case's name.
#' @return One row: the bundle, its sides and how its compared values came out;
#'   `NULL` (with a message of class `psychnets_verify_skipped`) when a
#'   reference package of the case is not installed.
export_case <- function(name) {
  live <- .verify_cases()
  stopifnot("`name` must be one exported case" =
              length(name) == 1L && name %in% c(names(live), names(.verify_retained)))
  requires <- if (name %in% names(live)) live[[name]]$requires else .verify_retained[[name]]$requires
  absent <- requires[!vapply(requires, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    message(structure(
      class = c("psychnets_verify_skipped", "message", "condition"),
      list(message = sprintf("Case `%s` skipped: %s not installed.\n", name, toString(absent)), call = NULL)))
    return(NULL)
  }
  if (name %in% names(live)) .export_live(name, live[[name]]) else .verify_retained[[name]]$export(name)
}

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- c(names(.verify_cases()), names(.verify_retained))
print(do.call(rbind, lapply(requested, export_case)), row.names = FALSE)
