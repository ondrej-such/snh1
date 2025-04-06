library(tidyr)
library(dplyr)
library(xtable)



df <- read.csv("data/exp2.csv")

df1 <- df |> group_by(dataset) |> summarize(cv_acc = mean(cv_acc), wlw2 = mean(wlw2))
# colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-multi.tex")
print(xtable(df1, label = "tab:exp2", digits = 2,
             caption = "Comparison of combinations of trees and Wu-Lin-Weng's method"
),  include.rownames = FALSE)
sink()

