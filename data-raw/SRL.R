# Generates the SRL_<LLM> datasets: a random 300-row sample of the five MSLQ
# construct scores each "valid" (structured) LLM produced in the Computers in
# Human Behavior (2025) study. The two incoherent generators (ChatGPT, LeChat,
# whose items carry no estimable structure) are excluded. Source CSV is
# local-only; this script records how the shipped data was made.
src <- "/Users/mohammedsaqr/Downloads/rstudio-export (18) 2/mohammed_all.csv"
d <- read.csv(src, check.names = FALSE)

# Construct score = the mean of a construct's items (the item suffix names the
# construct). Returns a tidy 5-column data frame: CSU, IV, SE, SR, TA.
score_constructs <- function(df) {
  items <- grep("^Q_", names(df), value = TRUE)
  groups <- split(items, sub("^Q_[0-9]+_", "", items))
  as.data.frame(lapply(groups, function(cols) rowMeans(df[cols])))
}

valid <- c("GPT", "Gemini", "Claude", "Mistral", "LLaMa")

set.seed(2025)
for (llm in valid) {
  rows <- which(d$Source == llm)
  take <- sort(sample(rows, 300L))
  df <- score_constructs(d[take, ])
  rownames(df) <- NULL
  assign(paste0("SRL_", llm), df)
  save(list = paste0("SRL_", llm),
       file = file.path("data", paste0("SRL_", llm, ".rda")), compress = "xz")
}
