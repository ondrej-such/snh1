library(ggplot2)
library(dplyr)
library(tidyr)

source("theme.R")

df <- read.csv("data/multi-acc.csv")  
df$method = if_else(df$method == "wlw2", "Wu-Lin-Weng", df$method)

plt1 <- ggplot(df, aes(x = log2(tol), 
                       y = acc, 
                       group = method, 
                       linetype = method,
                       color = method )) + 
  geom_line() + facet_wrap(~dataset) + 
  scale_x_continuous(breaks = c(-12, -8, -4)) +
  ylab("0-1 score")  + pub_theme

ggsave(
  filename = "graphs/exp1-plot1.pdf",
  plot = plt1,
  width = 5.5,
  height = 3.5,
  units = "in"
)

library(lattice)

plt2 <- xyplot(acc ~ log2(tol) | dataset, data = df, type= "l", group = method, auto.key = list(space = "top", columns =
2), par.settings =
list(superpose.line = list(lty = c(1, 2,3 ,4 ))), 
ylab = "0-1 score")

pdf("graphs/exp1-plot2.pdf", width = 5.5, height = 3.5)
print(plt2)
dev.off()

#print(plt1)
