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
  mutate(separable = optimal) |>
  select(dataset, K, separable)

sink("tab-sep.tex")
print(xtable(df1, label = "tab:sep", digits = 2,
             caption = "The benchmark datasets with indication of the number of classes $K$ and percentage of binary problems deemed separable",
),  include.rownames = FALSE)
sink()
