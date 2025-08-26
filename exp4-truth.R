library(ggplot2)
library(dplyr)
source("theme.R")

df <- read.csv("data/exp4.csv") |>  mutate(step2 = multi) |> mutate(step1 = binary)

df$class = as.factor(df$changed - 1)

df1 <- df |> filter(dataset ==  "mnist") |> group_by(class, method) |> mutate(diff = mean(step2 - step1)) |> mutate(method = ifelse(method == "wlw2", "wlw", method))

plt1 <- ggplot(df1, aes(x = class, y = diff, fill = method, color = method)) +  facet_wrap(~ method) + 
	geom_col(data = df1, aes(x = class, y = diff, color = method), width = 0.7) +
	ylab("accuracy difference") + pub_theme +  guides(fill = guide_legend(keyheight = unit(1.2, "lines"))) +
	theme(legend.position = "none")
# legend.key.height = unit(1.2, "cm"),       # increase box height in legend
# legend.spacing.y = unit(0.6, "cm"))

ggsave(
    filename = "exp4-truth.pdf",
    plot = plt1,
    width = 5.5,
    height  = 2.4,
    units = "in"
    )
