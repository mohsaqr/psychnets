# tmfg_network(), logo_network() and tmfg_certificate(): the triangulated maximally filtered graph and the
# LoGo model of the ten bfi items' Pearson correlations, against NetworkToolbox::TMFG(), EGAnet::TMFG() and
# NetworkToolbox::LoGo(); the TMFG structure count against igraph and a base-R check of the defining
# properties.

items <- utils::read.csv(reference_data("data/bfi10.csv"))
correlations <- stats::cor(items)

#' The upper triangle of a reference's weighted adjacency matrix, diagonal dropped
weights_of <- \(m) {
  m <- as.matrix(m)
  dimnames(m) <- dimnames(correlations)
  ref_edges(m)
}

test_that("the TMFG equals NetworkToolbox::TMFG() and EGAnet::TMFG() on the same correlations", {
  skip_without_reference(c("NetworkToolbox", "EGAnet"))
  # The graph is chosen greedily from the same matrix and its weights are the correlations: identical up to
  # the last bit.
  expect_reference(as.data.frame(tmfg_network(data = items), include_zero = TRUE),
                   weights_of(NetworkToolbox::TMFG(correlations)$A), reference = "NetworkToolbox",
                   covers = "tmfg_network()", quantity = "TMFG edge weight", tolerance = 1e-12, key = c("from", "to"))
  expect_reference(as.data.frame(tmfg_network(data = items), include_zero = TRUE),
                   weights_of(EGAnet::TMFG(correlations, n = nrow(items))), reference = "EGAnet",
                   covers = "tmfg_network()", quantity = "TMFG edge weight", tolerance = 1e-12, key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(items, method = "tmfg"), include_zero = TRUE),
                   weights_of(NetworkToolbox::TMFG(correlations)$A), reference = "NetworkToolbox",
                   covers = "psychnet()", quantity = "TMFG edge weight, psychnet(method = \"tmfg\")",
                   tolerance = 1e-12, key = c("from", "to"))
})

test_that("the LoGo partial correlations equal NetworkToolbox::LoGo()", {
  skip_without_reference("NetworkToolbox")
  # Closed form from inverses of 4 x 4 and 3 x 3 blocks of one matrix.
  expect_reference(as.data.frame(logo_network(data = items), include_zero = TRUE),
                   weights_of(NetworkToolbox::LoGo(correlations, partial = TRUE)), reference = "NetworkToolbox",
                   covers = "logo_network()", quantity = "LoGo partial correlation", tolerance = 1e-10,
                   key = c("from", "to"))
  expect_reference(as.data.frame(psychnet(items, method = "logo"), include_zero = TRUE),
                   weights_of(NetworkToolbox::LoGo(correlations, partial = TRUE)), reference = "NetworkToolbox",
                   covers = "psychnet()", quantity = "LoGo partial correlation, psychnet(method = \"logo\")",
                   tolerance = 1e-10, key = c("from", "to"))
})

test_that("the LoGo precision equals the clique-separator formula of Barfuss et al. (2016)", {
  skip_without_reference("NetworkToolbox")
  tmfg <- NetworkToolbox::TMFG(correlations)
  pad <- \(block) {
    full <- matrix(0, ncol(correlations), ncol(correlations))
    full[block, block] <- solve(correlations[block, block])
    full
  }
  # Theta = sum over the 4-cliques of [S_cc^-1]^0 minus the sum over the 3-separators of [S_tt^-1]^0
  theta <- Reduce(`+`, lapply(seq_len(nrow(tmfg$cliques)), \(k) pad(tmfg$cliques[k, ]))) -
    Reduce(`+`, lapply(seq_len(nrow(tmfg$separators)), \(k) pad(tmfg$separators[k, ])))
  dimnames(theta) <- dimnames(correlations)
  expect_reference(as.data.frame(logo_network(data = items), include_zero = TRUE), ref_edges(-stats::cov2cor(theta)),
                   reference = "independent: the LoGo clique-separator formula on NetworkToolbox's TMFG cliques, base R",
                   covers = "logo_network()", quantity = "LoGo partial correlation", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("tmfg_certificate() counts no violated TMFG property, as igraph and a base-R check find", {
  skip_without_reference("igraph")
  fit <- tmfg_network(data = items)
  adjacency <- abs(as.matrix(fit$weights)) > 0
  graph <- igraph::graph_from_adjacency_matrix(adjacency * 1, mode = "undirected", diag = FALSE)
  p <- ncol(adjacency)
  triangles <- igraph::triangles(graph) |> as.integer() |> matrix(ncol = 3L, byrow = TRUE)
  # A stacked triangulation: no triangle has more than two common neighbours, and two common neighbours of a
  # triangle are never adjacent to each other (Massara et al. 2016).
  stacked <- all(vapply(seq_len(nrow(triangles)), \(k) {
    common <- which(colSums(adjacency[triangles[k, ], , drop = FALSE]) == 3L)
    length(common) <= 2L && !(length(common) == 2L && adjacency[common[[1L]], common[[2L]]])
  }, logical(1)))
  violations <- sum(igraph::ecount(graph) != 3L * (p - 2L), !igraph::is_chordal(graph)$chordal,
                    !igraph::is_connected(graph), !stacked)
  expect_reference(tmfg_certificate(fit), violations, reference = "igraph", covers = "tmfg_certificate()",
                   quantity = "number of violated TMFG properties", tolerance = 0,
                   note = "Edge count 3(p - 2), chordality and connectedness by igraph; the stacked-triangulation test in base R.")
  expect_reference(certificate(fit)$certificate, violations, reference = "igraph", covers = "certificate()",
                   quantity = "certificate of tmfg_network()", tolerance = 0)
})

# The graphs tmfg_certificate() is given below: TMFGs built by tmfg_network() on several correlation matrices,
# and graphs that are not TMFGs. tmfg_certificate() reads the graph a network stores as `adjacency`; the only way
# to give it a graph that tmfg_network() would never build is to put that graph in a network of the same size.
#' A TMFG network of p nodes whose stored graph is replaced by `adjacency`
as_network_of <- \(adjacency) {
  p <- ncol(adjacency)
  fit <- tmfg_network(cor_matrix = 0.5^abs(outer(seq_len(p), seq_len(p), "-")))
  dimnames(adjacency) <- dimnames(fit$adjacency)
  fit$adjacency <- adjacency
  fit
}
#' The complete graph on p nodes, as a logical adjacency matrix
complete_graph <- \(p) {
  adjacency <- matrix(TRUE, p, p)
  diag(adjacency) <- FALSE
  adjacency
}
#' A logical adjacency matrix of p nodes with the given edges (a two-column matrix of node numbers)
graph_of <- \(p, edges) {
  adjacency <- matrix(FALSE, p, p)
  adjacency[edges] <- TRUE
  adjacency[edges[, 2:1, drop = FALSE]] <- TRUE
  adjacency
}
#' The correlations of n draws of a p-variable three-factor model, seeded
factor_correlations <- \(p, n, seed) {
  set.seed(seed)
  loadings <- matrix(stats::runif(p * 3L, -0.8, 0.8), p, 3L)
  draws <- matrix(stats::rnorm(n * 3L), n, 3L) %*% t(loadings) + matrix(stats::rnorm(n * p), n, p)
  colnames(draws) <- sprintf("X%02d", seq_len(p))
  stats::cor(draws)
}

#' The icosahedron: 12 nodes, 30 edges (3(p - 2)), maximal planar and not chordal; its nodes are the cyclic
#' permutations of (0, +-1, +-golden ratio), two nodes adjacent at distance 2
icosahedron <- \() {
  golden <- (1 + sqrt(5)) / 2
  base <- as.matrix(expand.grid(0, c(-1, 1), c(-golden, golden)))
  points <- rbind(base, base[, c(3L, 1L, 2L)], base[, c(2L, 3L, 1L)])
  distance <- as.matrix(stats::dist(points))
  abs(distance - 2) < 1e-9
}

bfi_tmfg <- tmfg_network(data = items)
bfi_graph <- bfi_tmfg$adjacency
bfi_edges <- which(upper.tri(bfi_graph) & bfi_graph, arr.ind = TRUE)
bfi_gaps <- which(upper.tri(bfi_graph) & !bfi_graph, arr.ind = TRUE)
#' The bfi TMFG with edges removed and added (rows of bfi_edges / bfi_gaps)
bfi_changed <- \(remove = integer(), add = integer()) {
  adjacency <- bfi_graph
  adjacency[bfi_edges[remove, , drop = FALSE]] <- FALSE
  adjacency[bfi_edges[remove, 2:1, drop = FALSE]] <- FALSE
  adjacency[bfi_gaps[add, , drop = FALSE]] <- TRUE
  adjacency[bfi_gaps[add, 2:1, drop = FALSE]] <- TRUE
  as_network_of(adjacency)
}

certificate_networks <- list(
  "TMFG, p = 4" = tmfg_network(cor_matrix = 0.5^abs(outer(1:4, 1:4, "-"))),
  "TMFG, p = 6" = tmfg_network(cor_matrix = 0.3 + 0.7 * diag(6) + 0.2 * (abs(outer(1:6, 1:6, "-")) == 1)),
  "TMFG, bfi10 (p = 10)" = bfi_tmfg,
  "TMFG, three-factor model, p = 20" = tmfg_network(cor_matrix = factor_correlations(20L, 300L, seed = 11L)),
  "TMFG, three-factor model, p = 45" = tmfg_network(cor_matrix = factor_correlations(45L, 500L, seed = 12L)),
  "bfi TMFG, first edge removed" = bfi_changed(remove = 1L),
  "bfi TMFG, last edge removed" = bfi_changed(remove = nrow(bfi_edges)),
  "bfi TMFG, first non-edge added" = bfi_changed(add = 1L),
  "bfi TMFG, last non-edge added" = bfi_changed(add = nrow(bfi_gaps)),
  "bfi TMFG, one edge moved to a non-edge" = bfi_changed(remove = 1L, add = 1L),
  "K5" = as_network_of(complete_graph(5L)),
  "K3,3" = as_network_of(graph_of(6L, as.matrix(expand.grid(1:3, 4:6)))),
  "K5 with a sixth node on one edge" = as_network_of(graph_of(6L, rbind(which(upper.tri(complete_graph(5L)), arr.ind = TRUE),
                                                                         c(1L, 6L), c(2L, 6L)))),
  "triangle with three common neighbours" = as_network_of(graph_of(6L, rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L),
                                                                              as.matrix(expand.grid(1:3, 4:6))))),
  "two disjoint K4" = as_network_of(graph_of(8L, rbind(which(upper.tri(complete_graph(4L)), arr.ind = TRUE),
                                                       which(upper.tri(complete_graph(4L)), arr.ind = TRUE) + 4L))),
  "K7 and two isolated nodes" = as_network_of(graph_of(9L, which(upper.tri(complete_graph(7L)), arr.ind = TRUE))),
  "6-cycle" = as_network_of(graph_of(6L, cbind(1:6, c(2:6, 1L)))),
  "octahedron" = as_network_of(complete_graph(6L) & !graph_of(6L, rbind(c(1L, 2L), c(3L, 4L), c(5L, 6L)))),
  "icosahedron" = as_network_of(icosahedron())
)

# The Python side: for each graph, its edge count, chordality, connectedness and planarity in NetworkX.
networkx_properties_script <- '
import json, sys
import networkx as nx
graphs = json.load(open(sys.argv[1]))
out = {"networkx": nx.__version__}
for name, spec in graphs.items():
    G = nx.Graph()
    G.add_nodes_from(range(1, spec["p"] + 1))
    G.add_edges_from(tuple(edge) for edge in spec["edges"])
    out[name] = {"edges": G.number_of_edges(), "chordal": nx.is_chordal(G),
                 "connected": nx.is_connected(G), "planar": nx.check_planarity(G)[0]}
json.dump(out, open(sys.argv[2], "w"))
'
#' Edge count, chordality, connectedness and planarity of each network's graph, by NetworkX
networkx_properties <- \(networks) {
  graphs <- lapply(networks, \(fit) {
    adjacency <- unname(fit$adjacency)
    list(p = ncol(adjacency), edges = unname(which(upper.tri(adjacency) & adjacency, arr.ind = TRUE)))
  })
  files <- c(input = tempfile(fileext = ".json"), script = tempfile(fileext = ".py"), output = tempfile(fileext = ".json"))
  on.exit(unlink(files), add = TRUE, after = FALSE)
  jsonlite::write_json(graphs, files[["input"]], auto_unbox = TRUE, digits = NA)
  writeLines(networkx_properties_script, files[["script"]])
  status <- system2(verify_python(), c("-I", files[["script"]], files[["input"]], files[["output"]]))
  stopifnot("NetworkX did not run" = identical(status, 0L))
  found <- jsonlite::read_json(files[["output"]])
  stopifnot("NetworkX is not the version verify.json states (3.7)" = identical(found$networkx, "3.7"))
  found
}
#' Skip, recorded, when the interpreter of verify_python() cannot import NetworkX (reticulate may be bound to
#' another Python, so the helper's own module check is not asked)
skip_without_networkx <- \() {
  python <- verify_python()
  found <- nzchar(python) &&
    identical(system2(python, c("-c", shQuote("import networkx")), stdout = FALSE, stderr = FALSE), 0L)
  if (!found) skip_without_reference("networkx")
}
#' The violated TMFG properties of a graph of p nodes, from its edge count, chordality, connectedness, planarity
violated_properties <- \(p, edges, chordal, connected, planar) {
  c("edge count 3(p - 2)", "chordal", "connected", "planar")[c(edges != 3L * (p - 2L), !chordal, !connected, !planar)]
}

test_that("tmfg_certificate() counts the TMFG properties NetworkX finds violated, graph by graph", {
  skip_without_networkx()
  found <- networkx_properties(certificate_networks)
  expected <- \(name) {
    fit <- certificate_networks[[name]]
    p <- ncol(fit$adjacency)
    graph <- found[[name]]
    violated <- violated_properties(p, graph$edges, graph$chordal, graph$connected, graph$planar)
    list(fit = fit, count = length(violated), quantity = sprintf("number of violated TMFG properties: %s", name),
         note = sprintf("NetworkX %s: %d edges for p = %d; violated: %s.", found$networkx, graph$edges, p,
                        if (length(violated)) toString(violated) else "none"))
  }
  # A count of properties: an integer, compared exactly.
  invisible(lapply(lapply(names(certificate_networks), expected), \(graph) {
    expect_reference(tmfg_certificate(graph$fit), graph$count, reference = "networkx", covers = "tmfg_certificate()",
                     quantity = graph$quantity, tolerance = 0, note = graph$note)
  }))
})

test_that("tmfg_certificate() counts the TMFG properties igraph finds violated, graph by graph", {
  skip_without_reference("igraph")
  skip_without_networkx()
  planar <- vapply(networkx_properties(certificate_networks)[names(certificate_networks)], \(graph) graph$planar, logical(1))
  expected <- \(name) {
    fit <- certificate_networks[[name]]
    graph <- igraph::graph_from_adjacency_matrix(fit$adjacency * 1, mode = "undirected", diag = FALSE)
    p <- ncol(fit$adjacency)
    violated <- violated_properties(p, igraph::ecount(graph), igraph::is_chordal(graph)$chordal,
                                    igraph::is_connected(graph), planar[[name]])
    list(fit = fit, count = length(violated), quantity = sprintf("number of violated TMFG properties: %s", name),
         note = sprintf("Edge count, chordality and connectedness by igraph, planarity by NetworkX's check_planarity (igraph has no planarity test); violated: %s.",
                        if (length(violated)) toString(violated) else "none"))
  }
  invisible(lapply(lapply(names(certificate_networks), expected), \(graph) {
    expect_reference(tmfg_certificate(graph$fit), graph$count, reference = "igraph", covers = "tmfg_certificate()",
                     quantity = graph$quantity, tolerance = 0, note = graph$note)
  }))
})
