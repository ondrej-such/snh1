library(tidyr)
library(dplyr)
library(ggplot2)

source('lda.R')
source("theme.R")


set.seed(123)
pt <- par_triples(limit = 1000)
write.csv(pt$df, file = "data/exp2.csv")


