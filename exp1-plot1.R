library(ggplot2)
library(dplyr)
library(tidyr)

df <- read.csv("data/multi-acc.csv") 

plt1 <- ggplot(df, aes(x = log2(tol), 
                       y = acc, 
                       group = method, 
                       color = method )) + 
  geom_line() + facet_wrap(~dataset) + 
  scale_x_continuous(breaks = c(-12, -8, -4)) +
  ylab("0-1 score")

ggsave(
  filename = "graphs/exp1-plot1.pdf",
  plot = plt1,
  width = 5.5,
  height = 3.5,
  units = "in"
)

#print(plt1)
