# Actor/event networks: frequency conversion + the cross-sectional vs
# within/between branching.

test_that("event_frequencies counts actions per occasion", {
  ev <- data.frame(
    Actor   = rep(c("a", "b"), each = 6),
    Session = rep(rep(1:2, each = 3), 2),
    Action  = c("read", "quiz", "read", "quiz", "read", "note",
                "note", "note", "read", "read", "quiz", "quiz"))
  f <- event_frequencies(ev, session = "Session")
  expect_s3_class(f, "data.frame")
  expect_identical(nrow(f), 4L)                       # 2 actors x 2 sessions
  expect_true(all(c("actor", "note", "quiz", "read") %in% names(f)))
  fa1 <- f[f$actor == "a", ]                          # actor a, first session = row 1
  expect_equal(fa1$read[1], 2L); expect_equal(fa1$quiz[1], 1L)
  expect_equal(fa1$note[1], 0L)
  # without a session, one occasion per actor
  expect_identical(nrow(event_frequencies(ev)), 2L)
})

test_that("psychnet eventdata: only an actor -> one ordinary GGM", {
  set.seed(1)
  ev <- data.frame(
    Actor  = rep(paste0("s", 1:40), each = 25),
    Action = sample(c("read", "quiz", "note", "watch"), 1000, replace = TRUE))
  net <- psychnet(ev, actor = "Actor", action = "Action")
  expect_s3_class(net, "psychnet")
  expect_false(inherits(net, "psychnet_multilevel"))   # not within/between
})

test_that("psychnet eventdata: sessions -> within + between, separating structure", {
  set.seed(7)
  J <- 80; K <- 6; N <- J * K
  actor <- rep(seq_len(J), each = K)
  bf <- stats::rnorm(J)[actor]; wf <- stats::rnorm(N)
  am <- stats::rnorm(J)[actor]
  V1 <- am + 0.9 * bf + stats::rnorm(N, 0, .3)
  V2 <- am + 0.9 * bf + stats::rnorm(N, 0, .3)        # share BETWEEN factor
  V3 <- 0.9 * wf + stats::rnorm(N, 0, .3)
  V4 <- 0.9 * wf + stats::rnorm(N, 0, .3)             # share WITHIN factor
  X <- cbind(V1, V2, V3, V4); colnames(X) <- paste0("V", 1:4)

  ml <- psychnet(data.frame(id = actor, X), id = "id", standardize = FALSE)
  expect_s3_class(ml, "psychnet_multilevel")
  expect_s3_class(ml, "netobject_group")          # plots as a cograph group
  expect_equal(attr(ml, "n_actors"), J)
  ed <- as.data.frame(ml)
  has <- function(level, a, b) any(ed$level == level &
    ((ed$from == a & ed$to == b) | (ed$from == b & ed$to == a)))
  expect_true(has("within", "V3", "V4"))
  expect_true(has("between", "V1", "V2"))
  expect_lt(certificate(ml$within)$certificate, 1e-6)
  expect_lt(certificate(ml$between)$certificate, 1e-6)
})

test_that("psychnet eventdata: standardize -> a single de-clustered network", {
  set.seed(7)
  J <- 60; K <- 5; N <- J * K
  actor <- rep(seq_len(J), each = K)
  X <- matrix(stats::rnorm(J, 0, 2)[actor] + stats::rnorm(N), N, 3)
  X <- cbind(X, X[, 1] + stats::rnorm(N, 0, .3))      # a within-correlated pair
  colnames(X) <- paste0("V", 1:4)
  net <- psychnet(data.frame(id = actor, X), id = "id", standardize = TRUE)
  expect_s3_class(net, "psychnet")
  expect_false(inherits(net, "psychnet_multilevel"))
})

test_that("compute_sessions splits actors into sessions from time gaps", {
  set.seed(3)
  base <- as.POSIXct("2026-01-01", tz = "UTC")
  mk <- function(a) {
    # three bursts, within-burst gaps 60s, between-burst gaps 3600s (> 900)
    times <- c(base + cumsum(c(0, rep(60, 4))),
               base + 3600 + cumsum(c(0, rep(60, 4))),
               base + 7200 + cumsum(c(0, rep(60, 4))))
    data.frame(Actor = a, Time = times,
               Action = sample(c("read", "quiz", "note", "watch"), 15, TRUE))
  }
  ev <- do.call(rbind, lapply(paste0("s", 1:10), mk))
  # compute_sessions = FALSE -> time ignored -> one occasion per actor
  f0 <- event_frequencies(ev, time = "Time", compute_sessions = FALSE)
  expect_identical(nrow(f0), 10L)
  # default (compute_sessions = TRUE): gap > 900 -> 3 sessions per actor
  f1 <- event_frequencies(ev, time = "Time")
  expect_identical(nrow(f1), 30L)                     # 10 actors x 3 sessions
  expect_identical(length(unique(f1$actor)), 10L)
  # standardize = FALSE to get the within/between pair
  ml <- psychnet(ev, actor = "Actor", action = "Action", time = "Time",
                 standardize = FALSE)
  expect_s3_class(ml, "psychnet_multilevel")
  expect_equal(attr(ml, "n_actors"), 10)
})

test_that("psychnet keeps threshold/gamma/labels in their positional slots", {
  set.seed(4)
  X <- matrix(stats::rnorm(200 * 5), 200, 5)
  # third positional argument is still `threshold` (not `vars`)
  a <- psychnet(X, "glasso", 0.1)
  b <- psychnet(X, method = "glasso", threshold = 0.1)
  expect_equal(a$weights, b$weights)
  if (nrow(a$edges)) expect_true(all(abs(a$edges$weight) >= 0.1))
})

test_that("event data fits TMFG (no spurious n argument) and honors labels", {
  set.seed(5)
  ev <- data.frame(Actor = rep(paste0("s", 1:40), each = 15),
                   Action = sample(c("a", "b", "c", "d"), 600, replace = TRUE))
  expect_s3_class(psychnet(ev, actor = "Actor", method = "tmfg"), "psychnet")
  n <- psychnet(ev, actor = "Actor", labels = c("Alpha", "Beta", "Gamma", "Delta"))
  expect_identical(n$nodes$label, c("Alpha", "Beta", "Gamma", "Delta"))
  expect_error(psychnet(ev, actor = "Actor", labels = c("only", "three", "here")),
               "must match the number of features")
})

test_that("psychnet errors on an event log without an actor or id", {
  d <- data.frame(a = stats::rnorm(20), b = stats::rnorm(20))
  expect_error(psychnet(d, source = "eventdata"), "Supply `actor`")
})

test_that("psychnet vars selects columns tidily (name range, index, vector)", {
  set.seed(11)
  X <- matrix(stats::rnorm(150 * 5), 150, 5)
  colnames(X) <- c("joy", "fear", "calm", "anger", "trust")
  n_range <- psychnet(X, method = "pcor", vars = joy:calm)
  expect_identical(n_range$nodes$label, c("joy", "fear", "calm"))
  i_range <- psychnet(X, method = "pcor", vars = 1:3)
  expect_identical(i_range$nodes$label, c("joy", "fear", "calm"))
  vec_chr <- psychnet(X, method = "pcor", vars = c("anger", "trust"))
  expect_identical(vec_chr$nodes$label, c("anger", "trust"))
  vec_sym <- psychnet(X, method = "pcor", vars = c(anger, trust))
  expect_identical(vec_sym$nodes$label, c("anger", "trust"))
})

test_that("event-data GGMs honor cor_method (single and group paths)", {
  set.seed(1)
  df <- data.frame(id = 1:120, matrix(stats::rnorm(120 * 5), 120, 5))
  names(df)[-1] <- paste0("F", 1:5)
  pe <- psychnet(df, id = "id", method = "pcor", cor_method = "pearson")
  sp <- psychnet(df, id = "id", method = "pcor", cor_method = "spearman")
  expect_gt(max(abs(pe$weights - sp$weights)), 1e-6)        # cor_method changes it

  g <- data.frame(id = 1:200, grp = rep(c("A", "B"), 100),
                  matrix(stats::rnorm(200 * 4), 200, 4))
  names(g)[-(1:2)] <- paste0("F", 1:4)
  gpe <- psychnet(g, id = "id", group = "grp", method = "pcor", cor_method = "pearson")
  gsp <- psychnet(g, id = "id", group = "grp", method = "pcor", cor_method = "spearman")
  expect_gt(max(abs(gpe[["A"]]$weights - gsp[["A"]]$weights)), 1e-6)  # dots forwarded
})
