library(tidyr)
library(dplyr)
library(xtable)

df <- read.csv("data/multi-acc.csv")

df1 <- df |> group_by(dataset) |> summarize(normal = mean(normal)/1000, 
                                            radial = mean(radial)/1000, 
                                            wlw2 = mean(wlw2)/1000)
colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-multi.tex")
print(xtable(df1, label = "tab:multi", digits = 2,
             caption = "Comparison of 0-1 score for three parameterless coupling methods"
),  include.rownames = FALSE)
sink()
