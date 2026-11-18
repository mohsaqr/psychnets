library(psychnets)
panel <- utils::read.csv("../data/panel.csv")
fit <- panel_gvar(data = panel, vars = c("V1", "V2", "V3"), id = "id", wave = "wave", select = "none")
fit
parameters <- as.data.frame(fit, what = "parameters")
parameters
indices <- utils::stack(as.data.frame(fit, what = "fit"),
                        select = c(logl, chisq, baseline_chisq, bic, cfi, tli, rmsea))
result <- rbind(
  with(parameters, data.frame(quantity = paste("estimate:", matrix, row, col), kind = "estimate", value = estimate)),
  with(parameters, data.frame(quantity = paste("se:", matrix, row, col), kind = "standard error", value = se)),
  with(indices, data.frame(quantity = paste("fit:", ind), kind = ifelse(ind %in% c("cfi", "tli", "rmsea"), "fit index", "fit statistic"), value = values)))
result
