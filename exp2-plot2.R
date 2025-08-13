library(tidyr)
library(dplyr)
library(ggplot2)

source("theme.R")


df1 <- read.csv("data/exp2.csv") |> 
    mutate(outcome = if_else(cv_acc >= wlw2, "gt", "lt")) |> 
    group_by(dataset, outcome) |>
    summarize(count = n())

df1$outcome = as.factor(df1$outcome)

plt1 <- ggplot(df1, aes(x = dataset, y = count, fill = outcome)) + 
    geom_col() + theme_minimal(base_size = 9) + coord_flip() + pub_theme + 
    	scale_fill_manual(
          values = c("gt" = "seagreen", "lt" = "brown"),
	  labels = c(expression(parametric >= bolditalic(wlw)),
		     expression(parametric < bolditalic(wlw))),
	  name = "accuracy comparison"
	  )

ggsave(
  filename = "exp2-plot2.pdf",
    plot = plt1,
    width = 5.5,   # Width in inches
    height = 2.5,  # Height in inches
    units = "in"   # Specify units as inches
)

