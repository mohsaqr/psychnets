# accessors return tidy data frames

    Code
      print(fit)
    Output
      <psychnet_panelgvar> panel graphical VAR (temporal + contemporaneous + between)
        persons: 400   waves: 4   nodes: 3
        edges (significant/estimated, p < 0.01): temporal 2/6 (+3/3 autoregressive)   contemporaneous 1/3   between 1/3
        selection: none (full model)
        fit: chisq(66) = 70.97, p = 0.316   CFI = 0.989   TLI = 0.986   RMSEA = 0.014
        optimality (ML score): <certificate>
        starts: 3 of 3 reached the best saturated optimum

# EBIC estimates carry no p-values; refit = TRUE supplies them

    Code
      print(fit)
    Output
      <psychnet_panelgvar> panel graphical VAR (temporal + contemporaneous + between)
        persons: 400   waves: 4   nodes: 3
        edges (EBIC-selected): temporal 4 (+3 autoregressive)   contemporaneous 2   between 1
        selection: EBIC lasso (gamma = 0.5, lambda = 0.09805)
        fit: chisq(71) = 80.97, p = 0.196   CFI = 0.978   TLI = 0.975   RMSEA = 0.019
        optimality (lasso KKT residual): <certificate>
        starts: 3 of 3 reached the best penalized optimum

