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
        "-", "-", "-", "-", "-", "-"),
      c("cor_matrix", "n", "gamma", "nlambda", "lambda_min_ratio", "threshold",
        "penalize_diagonal", "refit", "-", "-", "-", "-", "-",
        "data", "cor_method", "na_method", "native", "labels", "lambda_path"),
      c("renamed", "identical", "identical", "identical", "renamed",
        "semantics differ", "renamed", "semantics differ",
        "reference only", "reference only", "reference only",
        "reference only", "reference only",
        "psychnet only", "psychnet only", "psychnet only", "psychnet only",
        "psychnet only", "psychnet only"),
      c("correlation matrix in; psychnet also accepts raw data", "",
        "EBIC tuning, default 0.5", "default 100", "same default 0.01",
        "qgraph: logical sig-rule; psychnet: numeric weight cutoff >= 0 (a logical is refused)",
        "same default FALSE (W_ii = S_ii)",
        "qgraph refit = TRUE (unpenalized refit of the selected graph) is psychnet refit = \"unregularized\"; psychnet refit = TRUE (default) re-solves the penalized optimum to machine precision",
        "psychnet always validates/PD-checks cor_matrix",
        "per-edge penalty matrix; psychnet uses scalar rho",
        "EBIC df toggle; psychnet never counts the diagonal",
        "psychnet always stores $kkt, $lambda, ...", "psychnet is silent",
        "estimate straight from raw data", "pearson/spearman/kendall/auto",
        "pairwise/listwise missing data", "native solver (TRUE) / glasso Fortran (FALSE)",
        "node labels",
        "fixed penalty grid, for resampling callers")),

    # ---- qgraph::cor_auto  <->  cor_auto ---------------------------------
    .cw("qgraph::cor_auto", "cor_auto",
      c("data", "ordinalLevelMax", "missing", "detectOrdinal", "forcePD",
        "npn.SKEPTIC", "select", "verbose"),
      c("data", "ordinal_max_levels", "na_method", "-", "-", "-", "-", "-"),
      c("semantics differ", "renamed", "renamed", "reference only",
        "reference only", "reference only", "reference only", "reference only"),
      c(paste("ordered factors are ordinal (unordered factors are dropped with a",
              "message); the polyserial is the Olsson ad-hoc estimator, not",
              "qgraph's two-step ML (differences up to ~0.006)"),
        "same default 7", "same default pairwise",
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
        "both select on glmnet's default lambda path",
        "lambda-path length (glmnet default 100)",
        "smallest penalty fraction (glmnet default: 1e-4, 0.01 if n < predictors)",
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
        "psychnet always lasso-regularized", "both use glmnet's default path",
        "psychnet is EBIC-only (no CV)", "no CV folds in psychnet",
        "elastic-net alpha; psychnet is pure lasso", "no alpha selection",
        "no alpha CV", "no alpha EBIC", "max interaction order; psychnet pairwise",
        "single moderator supported via moderators=", "estimation method (EBIC fixed)",
        "binary sign handling", "psychnet scales gaussian nodes internally",
        "no overparameterized coding", "categorical threshold",
        "sign-info messages", "verbose toggle", "progress bar",
        "warning toggle", "model saving", "data saving",
        "lambda-path length (glmnet default 100)",
        "smallest penalty fraction (glmnet default: 1e-4, 0.01 if n < predictors)",
        "pairwise/listwise", "node labels")),

    # ---- psychonetrics::panelgvar (+ runmodel / prune)  <->  panel_gvar ---
    .cw("psychonetrics::panelgvar", "panel_gvar",
      c("data", "vars", "idvar", "beepvar", "missing", "alpha", "adjust",
        "recursive", "criterion", "estimator", "nLambda", "lambda_min_ratio",
        "standardize", "datatype", "within_latent", "between_latent",
        "stationary_means", "covtype", "optimizer", "-", "-", "-", "-", "-"),
      c("data", "vars", "id", "wave", "na_method", "alpha", "adjust",
        "recursive", "criterion", "select", "nlambda", "lambda_min_ratio",
        "standardize", rep("-", 6), "select", "stepup", "gamma", "refit",
        "labels"),
      c("identical", "identical", "renamed", "renamed", "semantics differ",
        "identical", "identical", "identical", "identical",
        "semantics differ", "renamed", "identical", "semantics differ",
        rep("reference only", 6), "psychnet only", "default differs",
        "psychnet only", "psychnet only", "psychnet only"),
      c("wide (one row per person) or long data",
        "design matrix (variables x waves); psychnet also takes a character vector for long data",
        "person id column (long data)", "wave column (long data); waves ordered by sort()",
        "psychonetrics: FIML when data have NA; psychnet: listwise only",
        "prune() argument, same default 0.01", "prune() argument, same default \"none\"",
        "prune() argument, same default FALSE",
        "stepup() argument, same default \"bic\"",
        "estimator = \"PML\" (penalized, auto lambda) is select = \"ebic\"; FIML is not supported",
        "find_penalized_lambda() argument, same default 50",
        "find_penalized_lambda() argument, same default 1e-4",
        "\"none\"/\"z\"/\"z_per_wave\" as psychonetrics; \"quantile\" is not offered",
        "psychnet infers wide vs long from `vars`",
        "psychnet fixes \"ggm\" (the panelgvar default)",
        "psychnet fixes \"ggm\" (the panelgvar default)",
        "psychnet: stationary means (the psychonetrics default)",
        "psychnet: ML covariance (denominator n), as psychonetrics uses for raw data",
        "psychnet: nlminb + Fisher-scoring polish, ML score stored as a certificate",
        "\"none\" (default) = runmodel() full model with every edge's p-value; \"prune\" = prune() |> stepup(); \"ebic\" = PML with EBIC lambda",
        "stepup() after prune(); psychnet runs it by default (FALSE = prune() alone)",
        "EBIC gamma for select = \"ebic\" (psychonetrics encodes it in criterion = \"ebic.5\")",
        "unpenalized refit of the EBIC-selected network",
        "node labels"))
  )
}

#' Map the arguments of established network packages to psychnets
#'
#' @description
#' `net_crosswalk()` maps the arguments of `qgraph`, `IsingFit`, `mgm` and
#' `psychonetrics` functions to their psychnets counterparts, as a table with
#' one row per argument. It shows which function fits the same model, which
#' argument corresponds to which, and where the defaults or the meaning of an
#' argument differ, for analyses moved to this package from those packages.
#'
#' The table covers `qgraph::EBICglasso()` ([ebic_glasso()]),
#' `qgraph::cor_auto()` ([cor_auto()]), `qgraph::ggmModSelect()`
#' ([ggm_modselect()]), `IsingFit::IsingFit()` ([ising_fit()]), `mgm::mgm()`
#' ([mgm_fit()]) and `psychonetrics::panelgvar()` ([panel_gvar()]); `reference`
#' selects one of them. Each row records whether the argument is identical,
#' renamed, has a different default or a different meaning, or exists in only
#' one of the two functions, with a short note.
#'
#' The result is a data frame with one row per argument; nothing is fitted.
#'
#' Intensive time-series models (`graphicalVAR`, `mlVAR`) are not covered by
#' this package.
#'
#' @param reference The reference function to show: `"all"` (default),
#'   `"EBICglasso"`, `"cor_auto"`, `"ggmModSelect"`, `"IsingFit"`, `"mgm"` or
#'   `"panelgvar"`.
#' @return A data frame with one row per argument and columns `reference` (the
#'   reference function, as `pkg::fn`), `psychnet` (the corresponding function
#'   in this package), `ref_arg` and `psychnet_arg` (the argument names, `"-"`
#'   when there is no counterpart), `status` (`"identical"`, `"renamed"`,
#'   `"default differs"`, `"semantics differ"`, `"reference only"` or
#'   `"psychnet only"`) and `note`.
#' @section Conditions:
#'   `psychnet_bad_input` when `reference` is not one of the listed values.
#' @seealso [psychnet()], which also accepts the `qgraph` and `bootnet` method
#'   names.
#' @examples
#' net_crosswalk("EBICglasso")
#' net_crosswalk("IsingFit")
#' net_crosswalk("panelgvar")
#' @export
net_crosswalk <- function(reference = c("all", "EBICglasso", "cor_auto",
                                        "ggmModSelect", "IsingFit", "mgm",
                                        "panelgvar")) {
  reference <- .psn_match_arg(reference)
  cw <- .crosswalk_table()
  if (reference != "all") {
    cw <- cw[sub("^.*::", "", cw$reference) == reference, , drop = FALSE]
  }
  rownames(cw) <- NULL
  cw
}
