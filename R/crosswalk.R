# Argument crosswalk: psychnet as a drop-in substitute for the reference
# packages it reimplements (qgraph, IsingFit, mgm). A static, curated map --
# one row per argument -- so a user coming from those packages can find the
# psychnet equivalent (and see what psychnet adds or fixes by default).

# Small row builder; vectors recycle so each call is one pairing block.
#' @noRd
.cw <- function(reference, psychnet, ref_arg, psychnet_arg, status, note) {
  data.frame(reference = reference, psychnet = psychnet,
             ref_arg = ref_arg, psychnet_arg = psychnet_arg,
             status = status, note = note, stringsAsFactors = FALSE)
}

#' @noRd
.crosswalk_table <- function() {
  rbind(
    # ---- qgraph::EBICglasso  <->  ebic_glasso ----------------------------
    .cw("qgraph::EBICglasso", "ebic_glasso",
      c("S", "n", "gamma", "nlambda", "lambda.min.ratio", "threshold",
        "penalize.diagonal", "refit", "checkPD", "penalizeMatrix",
        "countDiagonal", "returnAllResults", "verbose",
        "-", "-", "-", "-", "-"),
      c("cor_matrix", "n", "gamma", "nlambda", "lambda_min_ratio", "threshold",
        "-", "-", "-", "-", "-", "-", "-",
        "data", "cor_method", "na_method", "native", "labels"),
      c("renamed", "identical", "identical", "identical", "renamed",
        "semantics differ", "reference only", "reference only",
        "reference only", "reference only", "reference only",
        "reference only", "reference only",
        "psychnet only", "psychnet only", "psychnet only", "psychnet only",
        "psychnet only"),
      c("correlation matrix in; psychnet also accepts raw data", "",
        "EBIC tuning, default 0.5", "default 100", "same default 0.01",
        "qgraph: logical sig-rule; psychnet: numeric weight cutoff",
        "psychnet always FALSE (W_ii = S_ii)",
        "psychnet always refits the selected lambda (two-tier)",
        "psychnet always validates/PD-checks cor_matrix",
        "per-edge penalty matrix; psychnet uses scalar rho",
        "EBIC df toggle; psychnet never counts the diagonal",
        "psychnet always stores $kkt, $lambda, ...", "psychnet is silent",
        "estimate straight from raw data", "pearson/spearman/kendall/auto",
        "pairwise/listwise missing data", "native solver (TRUE) / glasso Fortran (FALSE)",
        "node labels")),

    # ---- qgraph::cor_auto  <->  cor_auto ---------------------------------
    .cw("qgraph::cor_auto", "cor_auto",
      c("data", "ordinalLevelMax", "missing", "detectOrdinal", "forcePD",
        "npn.SKEPTIC", "select", "verbose"),
      c("data", "ordinal_max_levels", "na_method", "-", "-", "-", "-", "-"),
      c("identical", "renamed", "renamed", "reference only", "reference only",
        "reference only", "reference only", "reference only"),
      c("", "same default 7", "same default pairwise",
        "psychnet always auto-detects ordinal",
        "psychnet always nearest-PD projects",
        "nonparanormal SKEPTIC; see huge_network(npn=)",
        "variable subset selection", "psychnet is silent")),

    # ---- qgraph::ggmModSelect  <->  ggm_modselect ------------------------
    .cw("qgraph::ggmModSelect", "ggm_modselect",
      c("S", "n", "gamma", "stepwise", "start", "considerPerStep", "nCores",
        "checkPD", "criterion", "verbose",
        "-", "-", "-", "-", "-", "-", "-"),
      c("cor_matrix", "n", "gamma", "stepwise", "-", "-", "-", "-", "-", "-",
        "nlambda", "lambda_min_ratio", "threshold", "cor_method", "na_method",
        "native", "labels"),
      c("renamed", "identical", "identical", "identical",
        "reference only", "reference only", "reference only", "reference only",
        "reference only", "reference only",
        "psychnet only", "psychnet only", "psychnet only", "psychnet only",
        "psychnet only", "psychnet only", "psychnet only"),
      c("also accepts raw data", "", "aligned: default 0 (BIC) in both", "default TRUE",
        "psychnet seeds candidates from the glasso path",
        "single-edge add/drop in psychnet", "psychnet is serial", "always on",
        "psychnet is EBIC-only", "psychnet is silent",
        "candidate-path length", "smallest penalty fraction",
        "partial-correlation weight cutoff", "pearson/spearman/kendall/auto",
        "pairwise/listwise", "base / glasso", "node labels")),

    # ---- IsingFit::IsingFit  <->  ising_fit ------------------------------
    .cw("IsingFit::IsingFit", "ising_fit",
      c("x", "gamma", "AND", "min_sum", "family", "plot", "progressbar",
        "lowerbound.lambda", "-", "-", "-", "-", "-"),
      c("data", "gamma", "rule", "min_sum", "-", "-", "-", "-",
        "nlambda", "lambda_min_ratio", "weights", "na_method", "labels"),
      c("renamed", "identical", "renamed", "default differs", "reference only",
        "reference only", "reference only", "reference only",
        "psychnet only", "psychnet only", "psychnet only", "psychnet only",
        "psychnet only"),
      c("binary 0/1 data", "default 0.25",
        "AND=TRUE/FALSE -> rule='AND'/'OR'",
        "-Inf vs NULL (both keep all rows)",
        "psychnet is binomial-only", "psychnet never plots",
        "psychnet has no progress bar",
        "psychnet uses its own full lambda path",
        "lambda-path length", "smallest penalty fraction",
        "observation weights", "pairwise/listwise", "node labels")),

    # ---- mgm::mgm  <->  mgm_fit ------------------------------------------
    .cw("mgm::mgm", "mgm_fit",
      c("data", "type", "lambdaGam", "ruleReg", "threshold", "weights",
        "level", "regularize", "lambdaSeq", "lambdaSel", "lambdaFolds",
        "alphaSeq", "alphaSel", "alphaFolds", "alphaGam", "k", "moderators",
        "method", "binarySign", "scale", "overparameterize", "thresholdCat",
        "signInfo", "verbatim", "pbar", "warnings", "saveModels", "saveData",
        "-", "-", "-", "-"),
      c("data", "types", "gamma", "rule", "threshold", "weights",
        rep("-", 22),
        "nlambda", "lambda_min_ratio", "na_method", "labels"),
      c("identical", "renamed", "renamed", "renamed", "identical", "identical",
        rep("reference only", 22),
        "psychnet only", "psychnet only", "psychnet only", "psychnet only"),
      c("", "node types; psychnet 'g'/'c' (gaussian/binary)",
        "EBIC tuning, default 0.25", "'and'/'or' -> 'AND'/'OR'",
        "'LW'/'HW'/'none' thresholding rule", "observation weights",
        "explicit #levels; psychnet infers binary vs gaussian",
        "psychnet always lasso-regularized", "psychnet builds its own path",
        "psychnet is EBIC-only (no CV)", "no CV folds in psychnet",
        "elastic-net alpha; psychnet is pure lasso", "no alpha selection",
        "no alpha CV", "no alpha EBIC", "max interaction order; psychnet pairwise",
        "no moderation in psychnet", "estimation method (EBIC fixed)",
        "binary sign handling", "psychnet scales gaussian nodes internally",
        "no overparameterized coding", "categorical threshold",
        "sign-info messages", "verbose toggle", "progress bar",
        "warning toggle", "model saving", "data saving",
        "lambda-path length", "smallest penalty fraction",
        "pairwise/listwise", "node labels"))
  )
}

#' Argument crosswalk: psychnet as a substitute for qgraph / IsingFit / mgm
#'
#' A tidy, one-row-per-argument map from each reference package's estimator to
#' its `psychnet` equivalent, so users migrating from `qgraph::EBICglasso`,
#' `qgraph::cor_auto`, `qgraph::ggmModSelect`, `IsingFit::IsingFit`, or
#' `mgm::mgm` can find the matching argument and see what `psychnet` changes by
#' default. Cross-sectional estimators only (temporal models are out of scope).
#'
#' @param reference Which reference function to show: `"all"` (default),
#'   `"EBICglasso"`, `"cor_auto"`, `"ggmModSelect"`, `"IsingFit"`, or `"mgm"`.
#' @return A tidy `data.frame`, one row per argument, with columns `reference`
#'   (the `pkg::fn` being substituted), `psychnet` (the psychnet verb),
#'   `ref_arg`, `psychnet_arg` (`"-"` when there is no counterpart), `status`
#'   (`identical` / `renamed` / `default differs` / `semantics differ` /
#'   `reference only` / `psychnet only`), and a short `note`.
#' @examples
#' net_crosswalk("EBICglasso")
#' net_crosswalk("IsingFit")
#' @export
net_crosswalk <- function(reference = c("all", "EBICglasso", "cor_auto",
                                        "ggmModSelect", "IsingFit", "mgm")) {
  reference <- match.arg(reference)
  cw <- .crosswalk_table()
  if (reference != "all") {
    cw <- cw[sub("^.*::", "", cw$reference) == reference, , drop = FALSE]
  }
  rownames(cw) <- NULL
  cw
}
