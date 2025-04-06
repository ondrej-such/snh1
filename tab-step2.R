library(tidyr)
library(dplyr)
library(xtable)

df <- read.csv("data/triples.csv")

df1 <- df |> mutate(omit = pmax(omit12, omit13, omit23)) |> 
            summarize(`normal >= WLW` = mean(normal >= wlw2),
                      `omit >= WLW` = mean(omit >= wlw2),
                      `WLW > max(normal, omit)` = mean(wlw2 > normal & wlw2 > omit)) |>
      pivot_longer(everything(), names_to = 'outcome', values_to = "percentage")
            
                      
                      
                      
           
# colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-step2.tex")
print(xtable(df1, label = "tab:step2", digits = 2,
             caption = "Comparison of Wu-Lin-Weng's method and  parameterless coupling methods on three class subset"
),  include.rownames = FALSE)
sink()