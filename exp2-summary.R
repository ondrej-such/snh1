library(tidyr)
library(dplyr)
library(ggplot2)

source('lda.R')
source("theme.R")


pt <- par_triples(limit = 1000)
df1 <- pt$df |> 
    mutate(outcome = if_else(cv_acc >= wlw2, "parametric >= WLW", "parametric < WLW")) |> 
    group_by(dataset, outcome) |>
    summarize(count = n())

plt1 <- ggplot(df1, aes(x = dataset, y = count, color = outcome, fill = outcome)) + 
    geom_col() + theme_minimal(base_size = 9) + coord_flip() + pub_theme

write.csv(pt$df, file = "data/exp2.csv")


ggsave(
  filename = "exp2-summary.pdf",
    plot = plt1,
    width = 5.5,   # Width in inches
    height = 2.5,  # Height in inches
    units = "in"   # Specify units as inches
)

