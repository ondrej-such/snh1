
library(xtable)
library(readr)
library(dplyr)
library(e1071)
library(tidyr)

df <- read.csv("data/glm1-multi.csv")

summary <- df |> group_by(dataset, n, method) |> 
    mutate(correct_n = correct / if_else(n == 300, 500, 1000)) |>
    summarize(macc = mean(correct_n )) 

summary2 <- summary |> 
    pivot_wider(names_from = method, values_from = macc)

addtorow <- list()
addtorow$pos <- list(-1, 0)
addtorow$command <- c(
  "\\hline \\multicolumn{2}{|c|}{Datasets} & \\multicolumn{3}{c|}{Accuracy} \\\\ \\hline\n",
  "\\hline"
)

sink("glm1-multi.tex")
print(xtable(summary2, label = "glm1-multi1", digits = 2,
             caption = "Mean accuracy of multiclass model based on penalized logistic regerssion",
            ), 
      add.to.row = addtorow, include.rownames = FALSE)
sink()

