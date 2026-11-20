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
