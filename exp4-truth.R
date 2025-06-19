library(ggplot2)
library(dplyr)

df <- read.csv("data/exp4.csv") |> mutate(diff = multi-binary) |> mutate(step2 = multi) |> mutate(step1 = binary)
df$class = as.factor(df$changed - 1)

df1 <- df |> filter(dataset ==  "mnist")

plt1 <- ggplot(df1, aes(x = class, y = step2 - step1, color = method)) +  facet_wrap(~ method) + 
geom_col(data = df1, aes(x = class, y = step2 - step1, color = method), fill = "black", width = 0.7) +
theme_minimal() + ylab("count")

ggsave(
    filename = "exp4-truth.pdf",
    plot = plt1,
    width = 5.5,
    height  = 2.5,
    units = "in"
    )
