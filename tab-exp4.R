library(tidyr)
library(dplyr)
library(xtable)



df <- read.csv("data/exp4.csv")

df1 <- df |> group_by(dataset, method) |> mutate(method = ifelse(method == "wlw2", "$\\boldsymbol{wlw}$", method)) |>
    summarize(c2c1 = mean(c2c1)) |> pivot_wider(names_from = "dataset", values_from = c2c1)
# colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-exp4.tex")
print(xtable(df1, label = "tab:exp4", digits = 3,
             caption = "Comparison of the difference in accuracies between models $\\boldsymbol{c}_1$ and
             $\\boldsymbol{c}_2$"
),  include.rownames = FALSE, sanitize.text.function = identity)
sink()

