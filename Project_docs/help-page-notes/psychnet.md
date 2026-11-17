# Implementation notes removed from `?psychnet` (2026-11-17)

Still true of R/psychnet.R, R/groups.R and R/multilevel.R; kept for
maintainers.

- **Method aliases** (`.resolve_method()`): matching lower-cases the name and
  strips every non-letter, then looks it up in: cor/correlation -> cor;
  pcor/partial/partialcor -> pcor; glasso/ebicglasso/ebic -> glasso;
  huge/npn/nonparanormal -> huge; ggm/ggms/ggmmodselect/modselect/stepwise ->
  ggm; tmfg; logo; relimp/relativeimportance/lmg -> relimp; ising/isingfit ->
  ising; isingsampler -> ising_sampler; mgm. The stored `method` field is the
  short name.
- **Routing of shared arguments**: `threshold` goes only to
  `.threshold_methods` (cor, pcor, glasso, huge, ggm, logo), and to mgm when
  it is character; `gamma` only to `.gamma_methods` (glasso, huge, ggm, ising,
  mgm) and only when non-NULL. Method defaults: gamma 0.5 for glasso/huge, 0
  for ggm (as `qgraph::ggmModSelect`), 0.25 for ising/mgm.
  `.warn_ignored_args()` warns only for a non-default value (threshold != 0,
  gamma non-NULL). In group mode the warning is raised once at the front door
  and muffled in the per-level refits.
- **within**: implemented with `scale_within()` (mlVAR's `scaleWithin`) on all
  rows before grouping, so person means and SDs use every occasion. With
  `group`, each member stores the person ids so `net_compare(paired = TRUE)`
  pairs by person.
- **Event data**: `.net_eventdata()` allows only glasso, pcor, ggm, huge, tmfg
  and logo (the old help said "the Gaussian graphical methods"; pcor and tmfg
  are included, cor and relimp are not). The log is reduced with
  `event_frequencies()`; with several occasions per actor it is decomposed into
  within/between networks (`standardize = FALSE`) or person-centred into one
  network (`standardize = TRUE`).
- **source = "auto"** picks event data when `actor` or `id` is given (the old
  help mentioned only `actor`).
