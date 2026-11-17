# Implementation notes removed from `?tmfg_network` (2026-11-17)

Still true of R/tmfg.R; kept for maintainers.

- **Old framing**: "Equivalent in purpose to `NetworkToolbox::TMFG()` /
  `bootnet`'s `"TMFG"` default, pure base R; correctness is certified
  structurally by tmfg_certificate()."
- **Seed rule**: `order(rowSums(W * (W > mean(W))), decreasing = TRUE)[1:4]`
  on `W = abs(S)`, with the mean taken over the full matrix including the
  diagonal (as in NetworkToolbox) before the diagonal is zeroed. Insertion
  ties go to the first (vertex, face) pair found (strict `>`).
- **n**: the argument exists so a uniform `(cor_matrix =, n =)` call does not
  partial-match `na_method`.
- **CLAUDE.md** lists the TMFG edge weights as "signed (partial cor)"; the
  code stores `S * adjacency`, i.e. marginal correlations, as the old and new
  help say.
