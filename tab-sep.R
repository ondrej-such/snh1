library(tidyr)
library(dplyr)
library(xtable)

df <- read.csv("data/separation.csv")

df1 <- df |> group_by(dataset, status, K) |> summarize(count = n() ) |>
  mutate(pct = count / (10 * K * (K-1))) |>
  pivot_wider(id_cols = c("dataset", "K"), 
              names_from = "status", 
              values_from = "pct", 
              values_fill = 0) |>
  mutate(separable = sprintf("%.0f %%   ", 100 * optimal))  |>
  select(dataset, K, separable)

colnames(df1) <- c("Dataset", "$K$", "Percentage linearly separable")

sink("tab-sep.tex")
print(xtable(df1, label = "tab:sep", align = c("l", "l", "r", "r"), digits = 3,
             caption = "The benchmark datasets with indication of the number of classes $K$ and the percentage of binary classification problems that are linearly separable",
),  include.rownames = FALSE, sanitize.colnames.function = identity)
sink()
